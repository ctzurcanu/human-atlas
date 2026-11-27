"""Conforming regions of a supplied piecewise-linear multiclass field.

Each region is the maximum of one explicitly supplied score. No anatomical
identity or boundary is inferred here. Clipping uses one canonical position
bank: a cut vertex is defined by its minimal source simplex and tied labels.
Reciprocal polygons use the same triangulation. Internal calculation faces
cancel, leaving selectable region surfaces rather than overlapping envelopes.
The caller must independently audit the emitted meshes and anatomical fields.
"""
from collections import defaultdict
import time
import numpy as np


def partition_labelled_tetrahedra(positions,tetrahedra,scores,*,seconds=20.,maximum_tetrahedra=120000,
                                 identical_score_priority=None,allowed_tetrahedron_labels=None):
    """Partition supplied fields; an explicit priority resolves volume ties.

    ``identical_score_priority`` lists every label once, highest priority
    first. It applies only when two affine fields are literally identical on
    all four controls of a calculation tetrahedron. Ordinary interface ties
    remain shared. No score or anatomical boundary is perturbed.
    """
    p=np.asarray(positions,float);t=np.asarray(tetrahedra,int);s=np.asarray(scores,float)
    if p.ndim!=2 or p.shape[1]!=3 or t.ndim!=2 or t.shape[1]!=4 or s.ndim!=2 or len(s)!=len(p):
        raise ValueError('Indexed tetrahedra and one score vector per position required')
    if not np.isfinite(p).all()or not np.isfinite(s).all()or not 0<seconds<=300 or len(t)>maximum_tetrahedra:
        raise ValueError('Finite data and explicit work/deadline budgets required')
    if t.min()<0 or t.max()>=len(p)or s.shape[1]<2:
        raise ValueError('Invalid calculation indices or fewer than two labels')
    priority=None
    if identical_score_priority is not None:
        order=list(identical_score_priority)
        if any(not isinstance(i,(int,np.integer))or isinstance(i,(bool,np.bool_))for i in order)or sorted(order)!=list(range(s.shape[1])):
            raise ValueError('Identical-score priority must list every integer label exactly once')
        priority={int(label):rank for rank,label in enumerate(order)}
    determinant=np.linalg.det(np.transpose(p[t[:,1:]]-p[t[:,:1]],(0,2,1)))
    if np.any(determinant<=0):raise ValueError('All calculation tetrahedra must have positive orientation')
    restricted=None;restriction_proofs=[]
    if allowed_tetrahedron_labels is not None:
        if len(allowed_tetrahedron_labels)!=len(t):raise ValueError('One explicit allowed-label set per tetrahedron required')
        restricted=[tuple(map(int,row))for row in allowed_tetrahedron_labels];groups=defaultdict(list)
        for i,row in enumerate(restricted):
            if not row or len(set(row))!=len(row)or min(row)<0 or max(row)>=s.shape[1]:raise ValueError('Invalid allowed-label set')
            groups[row].append(i)
        for labels,ids in groups.items():
            omitted=sorted(set(range(s.shape[1]))-set(labels))
            for begin in range(0,len(ids),10000):
                q=s[t[np.asarray(ids[begin:begin+10000],int)]]
                if omitted:
                    # An affine score dominating each omitted field at all
                    # four controls also dominates it throughout the simplex.
                    # Strictness at an interior control excludes volume ties.
                    delta=q[:,:,labels[0],None]-q[:,:,omitted]
                    if np.any(delta<0)or np.any(np.max(delta,axis=1)<=0):
                        raise ValueError('An omitted label is not strictly volume-dominated by the supplied domain witness')
            restriction_proofs.append(dict(labels=list(labels),tetrahedra=len(ids),omitted_fields_affine_dominated=True))
    started=time.monotonic();deadline=started+seconds
    points=list(p.copy());weights=[{i:1.}for i in range(len(p))]
    ties=[frozenset()for _ in p];cache={};centroids={};maximum_replay=0.;maximum_weight_violation=0.
    faces_by_label=[{}for _ in range(s.shape[1])];raw_counts=np.zeros(s.shape[1],int)
    identical_ownership=[]

    def value(index,a,b):
        if a in ties[index]and b in ties[index]:return 0.
        return sum(w*(s[k,a]-s[k,b])for k,w in weights[index].items())

    def cut_vertex(a,b,owner,other):
        nonlocal maximum_replay,maximum_weight_violation
        da,db=value(a,owner,other),value(b,owner,other)
        if da==0:return a
        if db==0:return b
        if da*db>=0:raise ValueError('Requested a noncrossing calculation edge')
        support=tuple(sorted(set(weights[a])|set(weights[b])))
        # A categorical Dirichlet node can have several exactly equal losing
        # scores. Its plane incidences must be derived from the actual field,
        # even though the original node itself has no stored cut history.
        ta=frozenset(j for j in range(s.shape[1])if value(a,owner,j)==0.)
        tb=frozenset(j for j in range(s.shape[1])if value(b,owner,j)==0.)
        labels=tuple(sorted((ta&tb)|{owner,other}))
        key=(support,labels)
        if key in cache:return cache[key]
        nodes=np.array(support,int);label=np.array(labels,int)
        A=np.r_[np.ones((1,len(nodes))),(s[nodes[:,None],label[1:]]-s[nodes,label[0]][:,None]).T]
        rhs=np.r_[1.,np.zeros(len(label)-1)]
        if A.shape[0]<A.shape[1]:
            raise ValueError('Underdetermined tied field requires an explicit caller-side refinement: '+str(key))
        if A.shape[0]==A.shape[1]:w=np.linalg.solve(A,rhs)
        else:
            w,_,rank,_=np.linalg.lstsq(A,rhs,rcond=None)
            if rank<len(nodes):raise ValueError('Rank-deficient tied field: '+str(key))
        residual=np.max(abs(A@w-rhs));maximum_replay=max(maximum_replay,float(residual))
        maximum_weight_violation=max(maximum_weight_violation,float(max(0.,-w.min(),w.max()-1)))
        if residual>2e-11 or w.min()<-2e-11 or w.max()>1+2e-11:
            raise ValueError('Canonical intersection leaves its source simplex')
        # These are arithmetic gates, not a geometric weld or anatomical tolerance.
        q=w@p[nodes];index=len(points);points.append(q);weights.append(dict(zip(nodes.tolist(),w.tolist())))
        ties.append(frozenset(labels));cache[key]=index
        return index

    def polygon_order(ids,normal):
        ids=np.array(sorted(set(ids)),int)
        if len(ids)<3:raise ValueError('Collapsed region polygon')
        q=np.array([points[i]for i in ids]);center=q.mean(0);normal=np.asarray(normal,float)
        normal=normal/np.linalg.norm(normal);axis=np.eye(3)[np.argmin(abs(normal))]
        u=np.cross(normal,axis);u/=np.linalg.norm(u);v=np.cross(normal,u)
        angle=np.arctan2((q-center)@v,(q-center)@u);ids=ids[np.argsort(angle,kind='stable')]
        return np.roll(ids,-int(ids.argmin())).tolist()

    def add_polygon(label,polygon):
        if len(polygon)==3:triangles=[polygon]
        else:
            key=tuple(sorted(polygon))
            if key not in centroids:
                centroids[key]=len(points);points.append(np.array([points[i]for i in polygon]).mean(0));weights.append({});ties.append(frozenset())
            c=centroids[key];triangles=[[c,a,b]for a,b in zip(polygon,polygon[1:]+polygon[:1])]
        output=faces_by_label[label]
        for f in triangles:
            raw_counts[label]+=1;k=tuple(sorted(f))
            if k not in output:output[k]=f
            else:
                old=output.pop(k)
                if not any(np.array_equal(np.roll(old,i),np.array(f)[[0,2,1]])for i in range(3)):
                    raise ValueError('Calculation seam has duplicate rather than reciprocal faces')

    for number,tet in enumerate(t):
        if number%64==0 and time.monotonic()>deadline:raise TimeoutError('Explicit labelled-region construction deadline')
        local_scores=s[tet];active=[]
        candidates=restricted[number]if restricted is not None else range(s.shape[1])
        for owner in candidates:
            # A score dominated everywhere in this tetrahedron has no region.
            dominated=any(np.all(local_scores[:,j]>=local_scores[:,owner])and np.any(local_scores[:,j]>local_scores[:,owner])for j in candidates if j!=owner)
            if not dominated:active.append(owner)
        groups=[];unseen=set(active)
        while unseen:
            first=min(unseen)
            group=[j for j in sorted(unseen)if np.array_equal(local_scores[:,first],local_scores[:,j])]
            unseen.difference_update(group)
            if len(group)>1:groups.append(group)
        for group in groups:
            if priority is None:
                raise ValueError('Identical affine scores can duplicate positive volume in tetrahedron '+str(number)+
                                 '; supply explicit identical_score_priority for labels '+str(group))
            owner=min(group,key=priority.__getitem__)
            active=[label for label in active if label not in group or label==owner]
            identical_ownership.append(dict(tetrahedron=number,identical_labels=group,owner=owner))
        for owner in active:
            polygons=[tet[[1,2,3]].tolist(),tet[[0,3,2]].tolist(),tet[[0,1,3]].tolist(),tet[[0,2,1]].tolist()]
            for other in active:
                if owner==other:continue
                all_ids={i for face in polygons for i in face};d={i:value(i,owner,other)for i in all_ids}
                if all(x>=0 for x in d.values()):continue
                if all(x<=0 for x in d.values()):polygons=[];break
                new=[];cap=set()
                for face in polygons:
                    polygon=[]
                    for a,b in zip(face,face[1:]+face[:1]):
                        if d[a]>=0:polygon.append(a)
                        if d[a]*d[b]<0:
                            index=cut_vertex(a,b,owner,other);polygon.append(index);cap.add(index)
                        elif d[a]==0:cap.add(a)
                    if len(polygon)>=3:new.append(polygon)
                if len(cap)>=3:
                    gradient=np.linalg.solve(p[tet[1:]]-p[tet[0]],local_scores[1:,owner]-local_scores[1:,other]-local_scores[0,owner]+local_scores[0,other])
                    new.append(polygon_order(cap,-gradient))
                polygons=new
            for face in polygons:add_polygon(owner,face)
    result={label:np.asarray(list(faces.values()),int).reshape(-1,3)for label,faces in enumerate(faces_by_label)}
    v=np.asarray(points,float)
    for f in result.values():
        q=v[f]
        if np.any(np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)==0):raise ValueError('Degenerate emitted region triangle')
    return dict(positions=v,cells=result,evidence=dict(calculation_tetrahedra=len(t),labels=s.shape[1],
        original_controls_literal=np.array_equal(v[:len(p)],p),original_positions=len(p),positions=len(v),
        raw_tetrahedral_face_counts=raw_counts.tolist(),region_triangle_counts={k:len(f)for k,f in result.items()},
        maximum_score_equation_residual=maximum_replay,maximum_source_simplex_weight_violation=maximum_weight_violation,
        identical_affine_score_priority=None if priority is None else order,
        explicit_domain_restriction_proofs=restriction_proofs,
        identical_affine_score_ownership=identical_ownership,
        no_distance_based_coordinate_welding=True,shared_position_bank=True,anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-started,deadline_seconds=seconds))
