"""Free-control C2 cubic courses with a continuous curvature certificate.

For each cubic interval its derivative is a quadratic Bezier curve. Bounding
all three derivative controls below by a positive projection gives a speed
lower bound on the WHOLE interval. Its second derivative is linear, so two
endpoint SOC norm bounds bound the whole acceleration. Therefore
    curvature <= acceleration_bound / speed_lower_bound**2.
Intermediate points remain free unless explicit source-owned attachment jets
are supplied. Clinical/packing acceptance is separate.
"""
import time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy import sparse
import clarabel
from .tube_sections import unit


class C2Course:
    def __init__(self, reference, length_mm, controls, minimum_radius_mm, *, speed_lower=.85, trust_radius_mm=8., boundary_condition='natural'):
        self.length=float(length_mm);self.knots=np.linspace(0,self.length,controls)
        if boundary_condition not in ('natural','not-a-knot'):
            raise ValueError('Supported homogeneous course boundary conditions: natural, not-a-knot')
        self.boundary_condition=boundary_condition
        self.basis=CubicSpline(self.knots,np.eye(controls),axis=0,bc_type=boundary_condition)
        self.reference=np.asarray(reference,float);self.minimum_radius=float(minimum_radius_mm)
        self.speed_lower=float(speed_lower);self.acceleration_bound=self.speed_lower**2/self.minimum_radius
        self.trust_radius=float(trust_radius_mm);self.controls=controls
        delta=np.diff(self.reference,axis=0);self.directions=unit(delta)
        self.reference_spline=CubicSpline(self.knots,self.reference,axis=0,bc_type=boundary_condition)

    def solve(self, constraints, *, time_limit=5., step_center=None, step_radius_mm=None, equalities=()):
        """Constrain a course without changing source station parameters.

        Inequalities: (parameter mm, inward normal, minimum dot product mm).
        Equalities: (parameter mm, derivative order 0/1/2, xyz value, optional
        owner label). Derivatives are with respect to the physical parameter,
        so a tangent value includes its declared speed. A root/attachment jet
        can pin position, first and second derivative. Natural end conditions
        already require zero endpoint second derivative; use not-a-knot when
        prescribing a nonzero endpoint curvature jet. This is a local course
        certificate, not global branch separation or finite-wall acceptance.
        """
        began=time.monotonic();k=self.controls;nv=k*3;sample=np.linspace(0,self.length,max(64,4*k))
        equality_rows=[];equality_rhs=[];declared=[]
        for item in equalities:
            if len(item) not in (3,4):raise ValueError('Equality requires parameter, derivative, xyz and optional owner')
            parameter,derivative,value=item[:3];parameter=float(parameter)
            if not np.isfinite(parameter)or not 0<=parameter<=self.length or derivative not in (0,1,2):
                raise ValueError('Finite in-range course parameter and derivative order 0/1/2 required')
            value=np.asarray(value,float)
            if value.shape!=(3,)or not np.isfinite(value).all():raise ValueError('Finite equality xyz value required')
            equality_rows.append(np.kron(self.basis(parameter,int(derivative))[None],np.eye(3)))
            equality_rhs.extend(value)
            declared.append(dict(parameter_mm=parameter,derivative_order=int(derivative),value=value.tolist(),owner=item[3]if len(item)>3 else None))
        B=self.basis(sample);D=self.basis(sample,2);desired=self.reference_spline(sample)
        metric=B.T@B+80*(D.T@D)+np.eye(k)*1e-6
        P=sparse.kron(sparse.csc_matrix(metric),sparse.eye(3),format='csc');q=-(B.T@desired).ravel()
        blocks=[];rhs=[];cones=[]
        def add(A,b,cone):blocks.append(sparse.csc_matrix(A));rhs.append(np.asarray(b,float));cones.append(cone)
        endpoints=np.zeros((6,nv));endpoints[:3,:3]=np.eye(3);endpoints[3:,-3:]=np.eye(3)
        add(endpoints,np.r_[self.reference[0],self.reference[-1]],clarabel.ZeroConeT(6))
        if equality_rows:add(np.vstack(equality_rows),equality_rhs,clarabel.ZeroConeT(len(equality_rhs)))
        speed=[]
        for i,h in enumerate(np.diff(self.knots)):
            q0=self.basis(self.knots[i],1);q1=q0+h/2*self.basis(self.knots[i],2);q2=self.basis(self.knots[i+1],1)
            for coefficient in (q0,q1,q2):speed.append(-np.outer(coefficient,self.directions[i]).ravel())
        add(np.asarray(speed),np.full(len(speed),-self.speed_lower),clarabel.NonnegativeConeT(len(speed)))
        if constraints:
            A=[];b=[]
            for constraint in constraints:
                parameter,normal,value=constraint[:3]
                A.append(-np.outer(self.basis(parameter),normal).ravel());b.append(-value)
            add(np.asarray(A),b,clarabel.NonnegativeConeT(len(A)))
        for parameter in self.knots:
            coefficient=self.basis(parameter,2)
            matrix=np.vstack([np.zeros((1,nv)),-np.kron(coefficient[None],np.eye(3))])
            add(matrix,[self.acceleration_bound,0,0,0],clarabel.SecondOrderConeT(4))
        for i in range(1,k-1):
            matrix=np.zeros((4,nv));matrix[1:,3*i:3*i+3]=-np.eye(3)
            add(matrix,np.r_[self.trust_radius,-self.reference[i]],clarabel.SecondOrderConeT(4))
        if step_center is not None:
            center=np.asarray(step_center,float);radius=float(step_radius_mm)
            if center.shape!=(k,3)or not np.isfinite(center).all()or not np.isfinite(radius)or radius<=0:
                raise ValueError('Finite local-step centers and positive radius required')
            for i in range(1,k-1):
                matrix=np.zeros((4,nv));matrix[1:,3*i:3*i+3]=-np.eye(3)
                add(matrix,np.r_[radius,-center[i]],clarabel.SecondOrderConeT(4))
        settings=clarabel.DefaultSettings();settings.verbose=False;settings.max_iter=90;settings.time_limit=time_limit
        settings.tol_feas=1e-8;settings.tol_gap_abs=1e-8;settings.tol_gap_rel=1e-8
        solver=clarabel.DefaultSolver(sparse.triu(P).tocsc(),q,sparse.vstack(blocks,format='csc'),np.concatenate(rhs),cones,settings)
        answer=solver.solve();status=str(answer.status)
        if status not in ('Solved','AlmostSolved'):
            z=np.asarray(answer.z);A=sparse.vstack(blocks,format='csc');b=np.concatenate(rhs);offset=6+len(equality_rhs)+len(speed)
            native_dual=z[offset:offset+len(constraints)];order=np.argsort(abs(native_dual))[::-1][:12]
            witnesses=[dict(weight=float(native_dual[j]),constraint=constraints[j][3]if len(constraints[j])>3 else dict(parameter_mm=float(constraints[j][0])))for j in order]
            return None,dict(status=status,iterations=int(answer.iterations),dual_b_dot_z=float(b@z),dual_constraint_balance_norm=float(np.linalg.norm(A.T@z)),
                             leading_native_conic_witnesses=witnesses,declared_station_equalities=declared,scope='Certificate of the stated linearised conic problem only; not proof of global anatomical route impossibility',elapsed_seconds=time.monotonic()-began)
        controls=np.asarray(answer.x).reshape(k,3)
        controls[[0,-1]]=self.reference[[0,-1]]
        certificate=self.certificate(controls)
        replay=[dict(item,maximum_absolute_residual=float(np.max(np.abs(self.evaluate(controls,item['parameter_mm'],item['derivative_order'])-item['value']))))for item in declared]
        maximum=max([x['maximum_absolute_residual']for x in replay],default=0.)
        if maximum>1e-7:return None,dict(status='EqualityReplayFailed',solver_status=status,station_equalities=replay,maximum_station_equality_residual=maximum,elapsed_seconds=time.monotonic()-began)
        return controls,dict(status=status,iterations=int(answer.iterations),certificate=certificate,station_equalities=replay,maximum_station_equality_residual=maximum,elapsed_seconds=time.monotonic()-began)

    def evaluate(self,controls,parameters,derivative=0):
        """Evaluate at saved native physical parameters, retaining attachments."""
        return CubicSpline(self.knots,controls,axis=0,bc_type=self.boundary_condition)(parameters,derivative)

    def certificate(self,controls):
        curve=CubicSpline(self.knots,controls,axis=0,bc_type=self.boundary_condition);speed=[]
        for i,h in enumerate(np.diff(self.knots)):
            v0=curve(self.knots[i],1);v1=v0+h/2*curve(self.knots[i],2);v2=curve(self.knots[i+1],1)
            speed.extend(float(v@self.directions[i])for v in (v0,v1,v2))
        acceleration=np.linalg.norm(curve(self.knots,2),axis=1);minimum_speed=min(speed);maximum_acceleration=float(acceleration.max())
        radius=minimum_speed**2/maximum_acceleration if maximum_acceleration else np.inf
        return dict(whole_interval_speed_lower_bound=minimum_speed,whole_interval_acceleration_upper_bound_per_mm=maximum_acceleration,
                    whole_interval_bend_radius_lower_bound_mm=radius,requested_minimum_radius_mm=self.minimum_radius,
                    continuous_curvature_bound_pass=bool(minimum_speed>=self.speed_lower-1e-7 and maximum_acceleration<=self.acceleration_bound+1e-7 and radius>=self.minimum_radius-1e-7),
                    proof='Quadratic derivative Bezier convex hull projection and linear acceleration convex norm bounds on every complete C2 cubic interval')

    def sample(self,controls,*,chord_error_mm=.003,max_step_mm=.25):
        curve=CubicSpline(self.knots,controls,axis=0,bc_type=self.boundary_condition);points=[curve(0.)];parameters=[0.]
        def flatten(bezier,a,b,depth=0):
            chord=bezier[3]-bezier[0];fraction=np.clip((bezier-bezier[0])@chord/max(chord@chord,1e-30),0,1)
            error=np.linalg.norm(bezier-bezier[0]-fraction[:,None]*chord,axis=1).max()
            if error<=chord_error_mm and np.linalg.norm(chord)<=max_step_mm:
                points.append(bezier[3]);parameters.append(b);return
            if depth>=20:raise ValueError('C2 chord subdivision limit reached')
            one=(bezier[:-1]+bezier[1:])/2;two=(one[:-1]+one[1:])/2;middle=(two[0]+two[1])/2;parameter=(a+b)/2
            flatten(np.array([bezier[0],one[0],two[0],middle]),a,parameter,depth+1)
            flatten(np.array([middle,two[1],one[2],bezier[3]]),parameter,b,depth+1)
        for a,b in zip(self.knots[:-1],self.knots[1:]):
            h=b-a;bezier=np.array([curve(a),curve(a)+h/3*curve(a,1),curve(b)-h/3*curve(b,1),curve(b)])
            flatten(bezier,a,b)
        return np.asarray(points),unit(curve(np.asarray(parameters),1)),np.asarray(parameters)
