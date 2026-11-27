"""Bounded grid untangling of an explicitly selected shared control patch.

Clinical curves, external joins and all unselected anatomy remain immutable.
This creates a candidate only: complete surface/occupied-volume checks and
anatomical/render review are mandatory after solving.
"""
import itertools,time
import numpy as np
from .faithful_float32 import separating_axis_options
from .tube_sections import TubeGeometryError
from ..distance import triangle_pairs


def untangle(source,encoded,faces,pairs,corner_words,*,origin=(0,0,0),
             mutable_controls=(),original_controls=(),additional_controls=(),
             original_steps=1,new_steps=4,maximum_displacement_m=5e-8,time_limit=10.,
             feasibility_only=False,maximum_rounds=1,point_halfspaces=(),separation_axis_limit=4):
    """Find a certified bounded embedding, optionally without minimizing cost.

    Feasibility mode changes only the optimization objective. Every source,
    orientation, separation, unique-control and displacement check still runs.
    """
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    p=np.asarray(source,float);q=np.asarray(encoded,float);f=np.asarray(faces,int)
    orientation_only=pairs is None
    pairs=np.empty((0,2,3),int)if orientation_only else np.asarray(pairs,int)
    w=np.asarray(corner_words,float);o=np.asarray(origin,float)
    if p.shape!=q.shape or p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all()or not np.isfinite(q).all():raise TubeGeometryError('Finite corresponding control banks required')
    if pairs.ndim!=3 or pairs.shape[1:]!=(2,3)or(not len(pairs)and not orientation_only)or w.shape!=(len(f),3,3):raise TubeGeometryError('Diagnosed pairs and complete corner fields required')
    if orientation_only and not len(additional_controls):raise TubeGeometryError('Orientation-only solving requires explicitly diagnosed controls')
    if time_limit<=0 or maximum_displacement_m<=0:raise TubeGeometryError('Finite positive solve bounds required')
    if not isinstance(feasibility_only,bool):raise TubeGeometryError('feasibility_only must be an explicit boolean')
    if not isinstance(separation_axis_limit,int)or isinstance(separation_axis_limit,bool)or not 1<=separation_axis_limit<=4:raise TubeGeometryError('One to four explicit source separating-axis choices required')
    if not isinstance(maximum_rounds,int)or isinstance(maximum_rounds,bool)or not 1<=maximum_rounds<=32:
        raise TubeGeometryError('A finite refinement budget from one to thirty-two rounds is required')
    relation=np.empty(0,dtype=str)if orientation_only else triangle_pairs(p[pairs[:,0]],p[pairs[:,1]])['relation']
    if np.isin(relation,['crossing','coplanar-overlap']).any():raise TubeGeometryError('Repair the source intersections before grid encoding')
    start=time.monotonic();allowed=set(map(int,mutable_controls));original=set(map(int,original_controls))
    support=[]
    for item in point_halfspaces:
        if not isinstance(item,dict):raise TubeGeometryError('Explicit point halfspace records required')
        control=item.get('control');normal=np.asarray(item.get('normal'),float);point=np.asarray(item.get('point'),float);minimum=float(item.get('minimum',0.))
        if not isinstance(control,(int,np.integer))or isinstance(control,bool)or not 0<=control<len(p):raise TubeGeometryError('Halfspace control outside source bank')
        if normal.shape!=(3,)or point.shape!=(3,)or not np.isfinite(np.r_[normal,point,minimum]).all():raise TubeGeometryError('Finite halfspace normal, point and minimum required')
        length=float(np.linalg.norm(normal))
        if length==0:raise TubeGeometryError('A halfspace normal must be nonzero')
        support.append((int(control),normal/length,point,minimum/length))
    selected=sorted(set(map(int,np.r_[pairs.ravel(),np.asarray(additional_controls,int).ravel()]))&allowed)
    variables=[];lo=q.copy();step=np.zeros_like(q);upper=[];preferred=[]
    for control in selected:
        radius=original_steps if control in original else new_steps
        for axis in range(3):
            local=np.float32(q[control,axis]-o[axis]);unit=float(abs(np.spacing(local)))
            if unit==0 or not np.isfinite(unit):continue
            radius_here=min(int(radius),int(maximum_displacement_m/unit))
            if radius_here<1:continue
            sequence=float(local)+np.arange(-radius_here,radius_here+1)*unit
            valid=sequence.astype(np.float32).astype(float)==sequence
            center=radius_here;left=center;right=center
            while left>0 and valid[left-1]:left-=1
            while right+1<len(valid)and valid[right+1]:right+=1
            if right==left:continue
            low=float(sequence[left])+o[axis]
            lo[control,axis]=low;step[control,axis]=unit
            variables.append((control,axis));upper.append(right-left);preferred.append((p[control,axis]-low)/unit)
    if not variables:raise TubeGeometryError('No selected representable grid choices')
    count=len(variables);columns={key:i for i,key in enumerate(variables)}
    costs=[0.]*count;integrality=[1]*count;lower_bounds=[0.]*count;upper_bounds=list(map(float,upper))
    rows=[];limits=[]
    def add(coeff,minimum=-np.inf,maximum=np.inf):
        coeff={i:v for i,v in coeff.items()if v!=0}
        if not coeff:
            if minimum>0 or maximum<0:raise TubeGeometryError('Immutable controls violate a required constraint')
            return
        scale=max(abs(v)for v in coeff.values());rows.append({i:v/scale for i,v in coeff.items()});limits.append((minimum/scale,maximum/scale))
    max_unit=float(step.max())
    for i,((control,axis),target)in enumerate(zip(variables,preferred)):
        auxiliary=len(costs);costs.append(step[control,axis]/max_unit*(10 if control in original else 1));integrality.append(0);lower_bounds.append(0.);upper_bounds.append(np.inf)
        add({auxiliary:1.,i:-1.},minimum=-target);add({auxiliary:1.,i:1.},minimum=target)
    for control,normal,point,minimum in support:
        coefficients={}
        for axis in range(3):
            column=columns.get((control,axis))
            if column is not None:coefficients[column]=normal[axis]*step[control,axis]
        add(coefficients,minimum=minimum-float((lo[control]-point)@normal))
    for a,b in pairs:
        selectors=[]
        for axis,gap in separating_axis_options(p[a],p[b],limit=separation_axis_limit):
            selector=len(costs);selectors.append(selector);costs.append(0.);integrality.append(1);lower_bounds.append(0.);upper_bounds.append(1.)
            for i,j in itertools.product(a,b):
                coeff={}
                for control,sign in ((int(i),-1),(int(j),1)):
                    for k in range(3):
                        column=columns.get((control,k))
                        if column is not None:coeff[column]=coeff.get(column,0.)+sign*axis[k]*step[control,k]
                rhs=gap*.25-float((lo[j]-lo[i])@axis)
                smallest=sum(min(0.,value*upper_bounds[column])for column,value in coeff.items())
                big=max(0.,rhs-smallest)
                coeff[selector]=-big;add(coeff,minimum=rhs-big)
        add({i:1. for i in selectors},minimum=1.)
    affected=np.flatnonzero(np.isin(f,selected).any(1));tri=p[f[affected]]
    normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);length=np.linalg.norm(normal,axis=1)
    if np.any(length==0):raise TubeGeometryError('Degenerate source facet')
    for row,face in enumerate(f[affected]):
        directions=[normal[row]/length[row]]
        for word in w[affected[row]]:
            size=np.linalg.norm(word)
            if size==0:raise TubeGeometryError('Zero source corner field')
            directions.append(word/size)
        vertices=tri[row]
        for direction in directions:
            area=float(normal[row]@direction)
            if area<=0:raise TubeGeometryError('Inward source corner field')
            gradient=np.array([np.cross(vertices[1]-vertices[2],direction),np.cross(vertices[2]-vertices[0],direction),np.cross(direction,vertices[1]-vertices[0])])
            constant=area;coeff={}
            for corner,control in enumerate(face):
                for k in range(3):
                    column=columns.get((int(control),k))
                    if column is not None:coeff[column]=gradient[corner,k]*step[control,k]
                    constant+=gradient[corner,k]*(lo[control,k]-p[control,k])
            add(coeff,minimum=area*1e-4-constant)
    active=np.unique(f)
    rejected=[]
    def forbid(controls,solution):
        ids=sorted({columns[(int(control),axis)] for control in controls for axis in range(3)if(int(control),axis)in columns})
        if not ids:raise TubeGeometryError('An invalid embedding has no mutable incident grid control')
        selectors=[]
        for column in ids:
            value=int(round(solution[column]));maximum=int(upper_bounds[column])
            if value>0:
                selector=len(costs);selectors.append(selector);costs.append(0.);integrality.append(1);lower_bounds.append(0.);upper_bounds.append(1.)
                big=maximum-value+1
                add({column:1.,selector:float(big)},maximum=float(value-1+big))
            if value<maximum:
                selector=len(costs);selectors.append(selector);costs.append(0.);integrality.append(1);lower_bounds.append(0.);upper_bounds.append(1.)
                big=value+1
                add({column:1.,selector:-float(big)},minimum=float(value+1-big))
        add({selector:1. for selector in selectors},minimum=1.)
    for attempt in range(maximum_rounds):
        rr=[];cc=[];entries=[]
        for i,row in enumerate(rows):
            for j,value in row.items():rr.append(i);cc.append(j);entries.append(value)
        remaining=float(time_limit)-(time.monotonic()-start)
        if remaining<=0:raise TimeoutError('Bounded grid refinement reached its deadline')
        matrix=coo_matrix((entries,(rr,cc)),shape=(len(rows),len(costs))).tocsc()
        objective=np.zeros(len(costs)) if feasibility_only else np.array(costs)
        solved=milp(objective,integrality=np.array(integrality),bounds=Bounds(lower_bounds,upper_bounds),
                    constraints=LinearConstraint(matrix,[x[0]for x in limits],[x[1]for x in limits]),options={'time_limit':remaining,'mip_rel_gap':1e-9})
        if not solved.success:
            error=TimeoutError if solved.status==1 else TubeGeometryError
            raise error('Bounded grid solve did not complete: '+str(solved.message))
        result=q.copy()
        for i,(control,axis)in enumerate(variables):result[control,axis]=lo[control,axis]+step[control,axis]*int(round(solved.x[i]))
        if not np.array_equal((result-o).astype(np.float32).astype(float)+o,result):raise TubeGeometryError('An output coordinate is not representable')
        actual=result[f[affected]];n=np.cross(actual[:,1]-actual[:,0],actual[:,2]-actual[:,0])
        bad=np.flatnonzero((np.einsum('ij,ij->i',n,normal)<=0)|np.any(np.einsum('fci,fi->fc',w[affected],n)<=0,axis=1))
        remaining_pairs=np.empty(0,dtype=str)if orientation_only else triangle_pairs(result[pairs[:,0]],result[pairs[:,1]])['relation']
        wrong=np.flatnonzero(np.isin(remaining_pairs,['crossing','coplanar-overlap']))
        _,inverse,counted=np.unique(result[active],axis=0,return_inverse=True,return_counts=True)
        duplicates=[active[inverse==i]for i in np.flatnonzero(counted>1)]
        bad_support=[i for i,(control,direction,point,minimum)in enumerate(support)if float((result[control]-point)@direction)<minimum]
        if not len(bad)and not len(wrong)and not duplicates and not bad_support:break
        rejected.append(dict(round=attempt+1,inward_or_collapsed_faces=len(bad),crossing_pairs=len(wrong),merged_control_groups=len(duplicates),violated_point_halfspaces=len(bad_support)))
        if attempt+1==maximum_rounds:raise TubeGeometryError('Bounded grid exact-validation refinement budget exhausted: '+str(rejected[-1]))
        for index in bad:forbid(f[affected[index]],solved.x)
        for index in wrong:forbid(pairs[index].ravel(),solved.x)
        for controls in duplicates:forbid(controls,solved.x)
        for index in bad_support:forbid([support[index][0]],solved.x)
    changed=np.flatnonzero(np.any(result!=q,axis=1));assert set(map(int,changed))<=allowed
    displacement=np.linalg.norm(result-q,axis=1)
    if np.max(np.abs(result-q))>maximum_displacement_m+1e-16:raise TubeGeometryError('Candidate exceeded its coordinate displacement limit')
    return result,dict(changed_controls=changed.tolist(),changed_original_controls=[int(i)for i in changed if int(i)in original],
            maximum_control_displacement_m=float(displacement.max()),maximum_coordinate_displacement_m=float(np.abs(result-q).max()),
            declared_coordinate_displacement_bound_m=maximum_displacement_m,source_position_reference_literal=True,
            diagnosed_pairs_checked=len(pairs),orientation_only=orientation_only,remaining_diagnosed_crossings=0,incident_faces_checked=len(affected),
            separation_axis_limit=separation_axis_limit,
            objective_mode='feasibility' if feasibility_only else 'minimum_displacement',
            solver_rounds=attempt+1,rejected_proposals=rejected,maximum_refinement_rounds=maximum_rounds,
            point_halfspace_checks=[dict(control=control,actual_signed_gap_m=float((result[control]-point)@direction),required_signed_gap_m=minimum,passed=True)for control,direction,point,minimum in support],
            variables=len(costs),constraints=len(rows),elapsed_seconds=time.monotonic()-start,full_mesh_acceptance=False)
