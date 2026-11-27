"""Index complete, opposite-wound vascular source openings without welding them.

The index describes geometric correspondence only. Repeated/ambiguous mouths,
different vascular systems and defective carriers cannot silently become joins.
"""
import numpy as np
from scipy.spatial import cKDTree
from .tube_continuity import boundary_loops
from .tube_sections import TubeGeometryError


def extract_mouths(mesh, *, owner, system):
    if system not in {'arterial','venous'}:raise ValueError('Explicit vascular system required')
    positions,rep,inverse=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True)
    try:rings=boundary_loops(inverse[mesh.faces])
    except TubeGeometryError as exc:
        return dict(owner=owner,complete=False,error=str(exc),mouths=[],frame=mesh.frame)
    mouths=[]
    for number,ring in enumerate(rings):
        p=positions[ring];center=p.mean(0);normal=np.cross(p-center,np.roll(p,-1,axis=0)-center).sum(0);area=float(np.linalg.norm(normal))*.5
        if not np.isfinite(p).all() or not area>1e-20:
            return dict(owner=owner,complete=False,error='Non-finite or zero-area boundary',mouths=[],frame=mesh.frame)
        mouths.append(dict(owner=owner,system=system,frame=mesh.frame,number=number,ids=rep[ring].tolist(),points=p.tolist(),center=center.tolist(),normal=(normal/(2*area)).tolist(),area_m2=area))
    return dict(owner=owner,complete=True,error=None,mouths=mouths,frame=mesh.frame)


def match_mouths(mouths, *, tolerance_m=1e-6, search_m=.002, area_tolerance=.01, maximum_pairs=50000):
    if not 0<tolerance_m<=search_m or not 0<=area_tolerance<1:raise ValueError('Invalid explicit tolerances')
    if not mouths:return dict(complete=True,examined_pairs=0,passed=[],near=[],ambiguous=[],unmatched=[])
    if len({p['frame']for p in mouths})!=1:raise ValueError('Mouths must share one explicit frame')
    tree=cKDTree(np.asarray([p['center']for p in mouths]))
    counts=tree.query_ball_point(tree.data,search_m,return_length=True)
    if (int(counts.sum())-len(mouths))//2>maximum_pairs:raise ValueError('Candidate-pair budget exceeded; split the saved index spatially')
    pairs=sorted(tree.query_pairs(search_m))
    passed=[];near=[]
    for i,j in pairs:
        first,second=mouths[i],mouths[j]
        if first['owner']==second['owner'] or first['system']!=second['system']:continue
        area_error=abs(second['area_m2']/first['area_m2']-1)
        cosine=float(np.asarray(first['normal'])@np.asarray(second['normal']))
        if area_error>area_tolerance or cosine>=-.99:continue
        p,q=np.asarray(first['points']),np.asarray(second['points']);n=len(p)
        # A reverse cyclic correspondence preserves the complete boundary order;
        # unconstrained nearest-neighbour assignments could scramble a mouth.
        curve_proof=None
        if len(p)==len(q):
            orders=[np.roll(np.arange(n)[::-1],phase)for phase in range(n)]
            gap,order=min(((float(np.linalg.norm(p-q[order],axis=1).max()),order)for order in orders),key=lambda item:item[0])
        else:
            from .closed_curve_correspondence import closed_curve_correspondence
            curve_proof=closed_curve_correspondence(p,q,opposite=True)
            gap=curve_proof['maximumGap_m'];order=np.arange(len(q))
        if gap>search_m:continue
        row=dict(first=first['owner'],second=second['owner'],firstPort=first['number'],secondPort=second['number'],firstRing=first['ids'],secondRing=np.asarray(second['ids'])[order].tolist(),maximumGap_m=gap,relativeAreaError=area_error,normalCosine=cosine,system=first['system'],completeSourceBoundaries=True,geometricCorrespondenceOnly=True,anatomicalAttachmentEstablished=False,wholeLumenPatencyEstablished=False)
        if curve_proof is not None:row['curveCorrespondence']=curve_proof
        (passed if gap<=tolerance_m else near).append(row)
    uses={}
    for r in passed:
        for key in [(r['first'],r['firstPort']),(r['second'],r['secondPort'])]:uses[key]=uses.get(key,0)+1
    ambiguous=[r for r in passed if any(uses[k]>1 for k in [(r['first'],r['firstPort']),(r['second'],r['secondPort'])])]
    ambiguous_ids={id(r)for r in ambiguous};passed=[r for r in passed if id(r)not in ambiguous_ids]
    paired={(r[key],r[key+'Port'])for r in passed for key in ['first','second']}
    unmatched=[dict(owner=p['owner'],number=p['number'])for p in mouths if (p['owner'],p['number'])not in paired]
    return dict(complete=True,examined_pairs=len(pairs),passed=passed,near=near,ambiguous=ambiguous,unmatched=unmatched,tolerance_m=tolerance_m,search_m=search_m,area_tolerance=area_tolerance,geometryChanged=False)


def declared_seam_slab(first,second,first_ring,second_ring,first_face,second_face,*,maximum_width_m):
    """Bound a declared joining-face contact by opposed complete halfspaces.

    This does not erase the raw triangle predicate. It returns an explicit
    bound for source-rounding contact at an already declared complete mouth.
    A face entering the other carrier beyond the narrow rim slab is rejected.
    """
    from geometry_tools.source_port_dependencies import check_source_mouth
    if not 0<maximum_width_m<=1e-8:raise ValueError('A source-rounding bound no larger than 10 nanometres is required')
    proof=check_source_mouth(first,second,first_ring,second_ring)
    if not proof['passed']:return dict(qualified=False,reason='Complete source-mouth proof failed')
    if not set(first.faces[first_face]).intersection(first_ring)or not set(second.faces[second_face]).intersection(second_ring):return dict(qualified=False,reason='Faces are not incident to the declared rim')
    rim=np.r_[first.vertices[first_ring],second.vertices[second_ring]];center=rim.mean(0);_,_,vh=np.linalg.svd(rim-center,full_matrices=False);normal=vh[-1]
    padding=128*np.finfo(float).eps*max(1.,float(np.abs(rim).max()));half=float(abs((rim-center)@normal).max())+padding
    a=(first.triangles[first_face]-center)@normal;b=(second.triangles[second_face]-center)@normal
    opposed=(a.min()>=-half and b.max()<=half and a.max()>half and b.min()<-half)or(b.min()>=-half and a.max()<=half and b.max()>half and a.min()<-half)
    result=dict(qualified=bool(opposed and 2*half<=maximum_width_m),maximumIntersectionSlabWidth_m=2*half,declaredWidthBound_m=maximum_width_m,completeMouthPassed=True,rawTriangleContactRetained=True,anatomicalAcceptance=False,method='planar-rim-halfspaces')
    if result['qualified']:return result
    # A valid cut through a curved source branch need not be planar. Use the
    # explicitly corresponding rim coordinates as an interval reference, never
    # as edited model geometry. Reject ideal triangle crossings and certify an
    # actual separating-axis overlap no thicker than source rounding permits.
    from geometry_tools.distance import triangle_pairs
    rims=[first.vertices[first_ring],second.vertices[second_ring]]
    if len(rims[0])==len(rims[1]):
        joint=(rims[0]+rims[1])*.5;targets=[joint,joint]
    else:
        # A T subdivision is a real corner of the source boundary, even when
        # the other carrier's same edge has no vertex there. Project to the
        # complete coarser reference polyline only for this certificate. The
        # model arrays are never edited, and the complete curve discrepancy
        # must itself satisfy the nanometre source-rounding bound.
        from .closed_curve_correspondence import closed_curve_correspondence
        ordered=[]
        for mesh,ids in [(first,first_ring),(second,second_ring)]:
            unique,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True)
            selected=set(inverse[np.asarray(ids,int)])
            cycle=next(x for x in boundary_loops(inverse[mesh.faces])if set(x)==selected)
            ordered.append(unique[cycle])
        rims=ordered
        curve=closed_curve_correspondence(*rims,opposite=True)
        if curve['maximumGap_m']>maximum_width_m:return dict(result,reason='Complete unequal-rim discrepancy exceeds source rounding',curveProof=curve)
        coarse=min(rims,key=len);edge=np.roll(coarse,-1,axis=0)-coarse;length2=np.einsum('ij,ij->i',edge,edge);targets=[]
        for ring in rims:
            fraction=np.clip(np.einsum('kij,ij->ki',ring[:,None]-coarse,edge)/length2,0,1)
            points=coarse+fraction[:,:,None]*edge
            nearest=np.linalg.norm(ring[:,None]-points,axis=2).argmin(1)
            targets.append(points[np.arange(len(ring)),nearest])
    actual=[first.triangles[first_face],second.triangles[second_face]];ideal=[]
    for triangle,ring,target_ring in zip(actual,rims,targets):
        q=triangle.copy()
        for point,target in zip(ring,target_ring):q[np.all(triangle==point,axis=1)]=target
        ideal.append(q)
    correction=max(float(np.linalg.norm(x-y,axis=1).max())for x,y in zip(actual,ideal))
    if correction>maximum_width_m:return dict(result,reason='Declared rim discrepancy exceeds the source-rounding bound')
    relation=str(triangle_pairs(ideal[0][None],ideal[1][None])['relation'][0])
    if relation in {'crossing','coplanar-overlap'}:return dict(result,reason='Ideal adjoining source faces still intersect',idealRelation=relation)
    edges=[np.roll(q,-1,axis=0)-q for q in ideal];normals=[np.cross(e[0],e[1])for e in edges]
    axes=list(normals)+[np.cross(x,y)for x in edges[0]for y in edges[1]]+[np.cross(n,e)for n,es in zip(normals,edges)for e in es]
    origin=np.concatenate(ideal).mean(0);proofs=[]
    for axis in axes:
        length=np.linalg.norm(axis)
        if length<=1e-30:continue
        axis=axis/length;ia,ib=[(q-origin)@axis for q in ideal];aa,bb=[(q-origin)@axis for q in actual]
        if ia.max()<=ib.min()+padding:
            overlap=max(0.,float(aa.max()-bb.min()));gap=float(ib.min()-ia.max())
        elif ib.max()<=ia.min()+padding:
            overlap=max(0.,float(bb.max()-aa.min()));gap=float(ia.min()-ib.max())
        else:continue
        bound=overlap+2*padding
        if bound<=maximum_width_m:proofs.append((bound,axis,gap))
    if not proofs:return dict(result,reason='No source-rounding separating-axis certificate',idealRelation=relation)
    bound,axis,gap=min(proofs,key=lambda p:p[0])
    return dict(result,qualified=True,method='curved-rim-separating-axis',maximumIntersectionSlabWidth_m=bound,maximumJointCorrection_m=correction,idealRelation=relation,idealSeparationGap_m=gap,separatingAxis=axis.tolist(),referenceOrigin=origin.tolist(),modelGeometryModified=False)
