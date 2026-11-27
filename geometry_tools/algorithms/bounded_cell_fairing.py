"""Sparse, bounded translation fairing of explicit rigid source cells.

The source dependency graph owns motion. Profile cells retain their gauge;
internal exporter ends are not fixed unless the caller declares them fixed.
Curvature constraints use actual adjacent profile-cell centres. Linearized
section halfspaces are proposals; emitted surfaces still require review.
"""
import time
import clarabel
import numpy as np
from scipy import sparse


def fair_rigid_cells(centres, edges, moving_cells, fixed_cells, targets, *,
                     profile_cells=(), halfspaces=(), plane_count=0,
                     maximum_motion_m=.008, curvature_per_m=80.,
                     time_limit_seconds=8., maximum_iterations=100,
                     infeasibility_callback=None, curvature_edges=None,
                     curvature_reference_centres=None, curvature_limits_per_m=None,
                     elastic_curvature_phase=False,surface_guards=None):
    start=time.monotonic();c=np.asarray(centres,float);edges=np.asarray(edges,int)
    if c.ndim!=2 or c.shape[1]!=3 or not np.isfinite(c).all()or edges.ndim!=2 or edges.shape[1]!=2:
        raise ValueError('Finite source centres and explicit graph edges required')
    if not 0<maximum_motion_m<=.050 or not 0<curvature_per_m<=200 or not 0<time_limit_seconds<=20:
        raise ValueError('Explicit finite motion, curvature and solver limits required')
    count=len(c);fixed=set(map(int,fixed_cells));free=np.array(sorted(set(map(int,moving_cells))-fixed),int);lookup={int(i):j for j,i in enumerate(free)};variables=2*len(free)+plane_count
    slack_index=variables if elastic_curvature_phase else None
    if elastic_curvature_phase:variables+=1
    reference=c if curvature_reference_centres is None else np.asarray(curvature_reference_centres,float)
    if reference.shape!=c.shape or not np.isfinite(reference).all():raise ValueError('Finite current geometric curvature reference required')
    limits={}if curvature_limits_per_m is None else {int(k):float(v)for k,v in curvature_limits_per_m.items()}
    if any(not 0<=k<count or not np.isfinite(v)or not 0<v<=200 for k,v in limits.items()):raise ValueError('Finite source-cell curvature quality limits up to200/m required')
    if not len(free):return np.zeros_like(c),dict(status='Solved',feasible=True,freeCells=0,elapsedSeconds=time.monotonic()-start)
    if np.any(edges<0)or np.any(edges>=count):raise ValueError('Source graph index outside declared cells')
    edges=np.unique(np.sort(edges,axis=1),axis=0);edges=edges[edges[:,0]!=edges[:,1]];length=np.linalg.norm(c[edges[:,0]]-c[edges[:,1]],axis=1);weight=1/np.maximum(length,.00001)**2;adj=sparse.coo_matrix((np.r_[weight,weight],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(count,count)).tocsr();degree=np.asarray(adj.sum(1)).ravel();L=sparse.eye(count,format='csr')-sparse.diags(np.divide(1.,degree,out=np.zeros(count),where=degree>0))@adj;Lf=L[:,free];H=(Lf.T@Lf+.002*sparse.eye(len(free))).tocsc();q=np.zeros(variables)
    target_rows=[]
    for key,entry in targets.items():
        if key not in lookup:continue
        delta=np.asarray(entry['translation_m'],float);strength=float(entry.get('weight',1.))
        if delta.shape!=(2,)or not np.isfinite(delta).all()or strength<=0 or np.linalg.norm(delta)>maximum_motion_m:raise ValueError('Bounded finite source-cell fitting target required')
        j=lookup[key];H[j,j]+=strength;q[2*j:2*j+2]-=strength*delta/.001;target_rows.append((int(key),delta,strength))
    P=sparse.block_diag([sparse.kron(H,sparse.eye(2)),sparse.eye(plane_count)*1e-10],format='csc')if plane_count else sparse.kron(H,sparse.eye(2),format='csc')
    if elastic_curvature_phase:P=sparse.block_diag([P,sparse.csc_matrix([[1e8]])],format='csc')
    blocks=[];rhs=[];cones=[];rr=[];cc=[];vv=[];bb=[]
    for row,halfspace in enumerate(halfspaces):
        for key,n in halfspace['coefficients'].items():
            if int(key)not in lookup:continue
            for axis in range(2):rr.append(row);cc.append(2*lookup[int(key)]+axis);vv.append(-float(n[axis]))
        if halfspace.get('plane')is not None:
            index,value=halfspace['plane'];rr.append(row);cc.append(2*len(free)+index);vv.append(-float(value))
        bb.append(-float(halfspace['minimum_m'])/.001)
    if len(bb):blocks.append(sparse.csc_matrix((vv,(rr,cc)),shape=(len(bb),variables)));rhs.append(np.asarray(bb));cones.append(clarabel.NonnegativeConeT(len(bb)))
    for j in range(len(free)):
        blocks.append(sparse.csc_matrix(([-1.,-1.],([1,2],[2*j,2*j+1])),shape=(3,variables)));rhs.append(np.array([maximum_motion_m/.001,0,0]));cones.append(clarabel.SecondOrderConeT(3))
    curve_edges=edges if curvature_edges is None else np.asarray(curvature_edges,int)
    if curve_edges.ndim!=2 or curve_edges.shape[1]!=2 or np.any(curve_edges<0)or np.any(curve_edges>=count):
        raise ValueError('Explicit source curvature-edge indices required')
    neighbour=[set()for _ in range(count)]
    for a,b in curve_edges:neighbour[a].add(int(b));neighbour[b].add(int(a))
    profile=set(map(int,profile_cells));curvature=[]
    for key in profile:
        if len(neighbour[key])!=2:continue
        a,b=sorted(neighbour[key],key=lambda j:c[j,1])
        if a not in profile or b not in profile or not c[a,1]<c[key,1]<c[b,1]:continue
        old_la,old_lb=np.linalg.norm(c[key]-c[a])/.001,np.linalg.norm(c[b]-c[key])/.001
        la,lb=np.linalg.norm(reference[key]-reference[a])/.001,np.linalg.norm(reference[b]-reference[key])/.001
        if min(la,lb,old_la,old_lb)<.1:continue
        oldbase=2*((c[b]-c[key])/.001/old_lb-(c[key]-c[a])/.001/old_la)/(old_la+old_lb)
        coefficients={a:2/(la*(la+lb)),key:-2/(la*lb),b:2/(lb*(la+lb))};base=sum(value*c[node]/.001 for node,value in coefficients.items());limit=max(limits.get(key,curvature_per_m)*.001,float(np.linalg.norm(oldbase))+1e-10)
        indices=[];columns=[];values=[]
        for node,value in coefficients.items():
            if node in lookup:
                indices.extend([1,3]);columns.extend([2*lookup[node],2*lookup[node]+1]);values.extend([-value,-value])
        if not values:continue
        if elastic_curvature_phase:indices.append(0);columns.append(slack_index);values.append(-1.)
        blocks.append(sparse.csc_matrix((values,(indices,columns)),shape=(4,variables)));rhs.append(np.r_[limit,base]);cones.append(clarabel.SecondOrderConeT(4));curvature.append((key,coefficients,base,limit,a,b))
    surface_count=0
    if surface_guards is not None:
        ids=np.asarray(surface_guards['cell_ids'],int);w=np.asarray(surface_guards['gradient_weights'],float);bounds=np.asarray(surface_guards['maximum_norm'],float);surface_count=len(ids)
        if ids.shape!=(surface_count,3)or w.shape!=(surface_count,3,2)or bounds.shape!=(surface_count,)or not np.isfinite(w).all()or not np.isfinite(bounds).all()or np.any(ids<0)or np.any(ids>=count)or np.any(bounds<=0)or np.any(bounds>=1):raise ValueError('Finite source-bound surface gradient cones required')
        lut=np.full(count,-1,int);lut[free]=np.arange(len(free));columns=lut[ids];valid=columns>=0;faceids=np.broadcast_to(np.arange(surface_count)[:,None],ids.shape);rows=[];cols=[];values=[]
        for axis in range(2):
            for derivative in range(2):
                rows.extend((5*faceids[valid]+1+2*axis+derivative).tolist());cols.extend((2*columns[valid]+axis).tolist());values.extend((-.001*w[:,:,derivative][valid]).tolist())
        blocks.append(sparse.csc_matrix((values,(rows,cols)),shape=(5*surface_count,variables)));bnd=np.zeros((surface_count,5));bnd[:,0]=bounds;rhs.append(bnd.ravel());cones.extend(clarabel.SecondOrderConeT(5)for _ in range(surface_count))
    if elastic_curvature_phase:
        blocks.append(sparse.csc_matrix(([-1.],([0],[slack_index])),shape=(1,variables)));rhs.append(np.array([0.]));cones.append(clarabel.NonnegativeConeT(1))
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.time_limit=time_limit_seconds;settings.max_iter=maximum_iterations;settings.tol_gap_abs=settings.tol_gap_rel=settings.tol_feas=1e-9
    A=sparse.vstack(blocks,format='csc');b=np.concatenate(rhs)
    solution=clarabel.DefaultSolver(sparse.triu(P,format='csc'),q,A,b,cones,settings).solve();status=str(solution.status);proof=dict(status=status,freeCells=len(free),fixedCells=len(fixed),softTargetCells=len(target_rows),sectionHalfspaces=len(halfspaces),curvatureConstraints=len(curvature),maximumMotion_m=maximum_motion_m,requestedCurvature_per_m=curvature_per_m,baselineHigherCurvaturePreservedRatherThanSilentlyClaimedCorrect=True,solverTimeSeconds=float(solution.solve_time),solverIterations=int(solution.iterations),timeLimitSeconds=time_limit_seconds,elapsedSeconds=time.monotonic()-start,anatomicalAcceptance=False)
    if status=='PrimalInfeasible' and infeasibility_callback is not None:
        z=np.asarray(solution.z,float);scale=-float(b@z)
        if np.isfinite(z).all() and scale>0:
            z=z/scale;at=0;cone_error=0.
            if len(halfspaces):
                cone_error=max(cone_error,-float(z[:len(halfspaces)].min()));at=len(halfspaces)
            for width,number in [(3,len(free)),(4,len(curvature)),(5,surface_count)]:
                for _ in range(number):
                    ray=z[at:at+width];cone_error=max(cone_error,float(np.linalg.norm(ray[1:])-ray[0]));at+=width
            if elastic_curvature_phase:cone_error=max(cone_error,-float(z[at]))
            certificate=dict(dualRhsDot=float(b@z),dualStationarityMax=float(np.abs(A.T@z).max(initial=0.)),dualConeViolation=max(0.,cone_error))
            certificate['verifiedNumerically']=certificate['dualStationarityMax']<1e-6 and certificate['dualConeViolation']<1e-6 and abs(certificate['dualRhsDot']+1)<1e-6
            proof['infeasibilityCertificate']=certificate
            infeasibility_callback(A,b,z,free,np.array([row[0]for row in curvature],int))
    if status!='Solved'and not(elastic_curvature_phase and status=='AlmostSolved'):
        proof['feasible']=False;return None,proof
    x=np.asarray(solution.x);delta=np.zeros_like(c);delta[free[:,None],np.array([0,2])]=x[:2*len(free)].reshape(-1,2)*.001;maximum=float(np.linalg.norm(delta,axis=1).max(initial=0));worst_curvature=0.;curvature_error=0.
    geometric_error=0.;maximum_geometric=0.;current=c+delta
    for key,coefficients,base,limit,a,b in curvature:
        actual=base+sum(value*delta[node]/.001 for node,value in coefficients.items());worst_curvature=max(worst_curvature,float(np.linalg.norm(actual))/.001);curvature_error=max(curvature_error,float(np.linalg.norm(actual)-limit)/.001)
        first=current[key]-current[a];second=current[b]-current[key];length_a=np.linalg.norm(first);length_b=np.linalg.norm(second)
        if min(length_a,length_b)<=1e-12:geometric_error=np.inf;maximum_geometric=np.inf;continue
        geometric=float(np.linalg.norm(2*(second/length_b-first/length_a)/(length_a+length_b)));maximum_geometric=max(maximum_geometric,geometric);geometric_error=max(geometric_error,geometric-limit/.001)
    residual=0.
    for row in halfspaces:
        value=sum(float(np.dot(n,delta[int(key),[0,2]]))for key,n in row['coefficients'].items())
        if row.get('plane')is not None:index,coefficient=row['plane'];value+=coefficient*x[2*len(free)+index]*.001
        residual=max(residual,row['minimum_m']-value)
    phase_slack=max(0.,float(x[slack_index])/.001)if elastic_curvature_phase else 0.
    surface_error=0.;surface_maximum=0.
    if surface_guards is not None:
        gradient=np.einsum('nka,nkb->nab',delta[ids][:,:,[0,2]],w);norm=np.linalg.norm(gradient,axis=(1,2));surface_error=float((norm-bounds).max(initial=0.));surface_maximum=float(norm.max(initial=0.))
    proof.update(surfaceGradientConstraints=surface_count,maximumSurfaceDisplacementGradient=surface_maximum,maximumSurfaceGradientViolation=surface_error)
    proposal_feasible=maximum<=maximum_motion_m+1e-8 and residual<=1e-8 and curvature_error<=phase_slack+1e-4 and surface_error<=1e-8
    proof.update(feasible=status=='Solved'and proposal_feasible and phase_slack<=1e-4,solverStatusQualified=status=='Solved',diagnosticProposalFeasible=proposal_feasible,elasticCurvaturePhase=bool(elastic_curvature_phase),curvaturePhaseSlack_per_m=phase_slack,maximumActualMotion_m=maximum,maximumHalfspaceViolation_m=max(0.,residual),maximumCheckedCurvature_per_m=worst_curvature,maximumCurvatureBoundViolation_per_m=max(0.,curvature_error),maximumActualGeometricCurvature_per_m=maximum_geometric,maximumGeometricCurvatureViolation_per_m=max(0.,geometric_error),geometricCurvatureQualified=geometric_error<=1e-4,calibreSpecificQualityLimits=bool(limits),currentArcLinearization=curvature_reference_centres is not None,softTargetRmsError_m=float(np.sqrt(np.mean([np.sum((delta[key,[0,2]]-desired)**2)for key,desired,w in target_rows])))if target_rows else 0.,sourceArcParameterBoundNotTrueCurvature=curvature_reference_centres is None,yCoordinatesLiteral=True,fixedCellsLiteral=True)
    return delta if proposal_feasible else None,proof
