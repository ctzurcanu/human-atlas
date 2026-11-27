"""Remove one caller-declared internal disk at an existing shared profile.

This certifies disk topology, opposing source walls and preserved boundaries.
It never classifies an anatomical valve or chooses a disk from proximity.
"""
from collections import Counter,defaultdict
import numpy as np
from .tube_continuity import continuity_audit,boundary_loops


def remove_shared_disk(vertices,faces,root_ids,cap_faces,*,minimum_wall_side_m=1e-8):
    p=np.asarray(vertices,float);f=np.asarray(faces,int);ids=np.asarray(cap_faces,int)
    physical,rep,inv=np.unique(p,axis=0,return_index=True,return_inverse=True);ff=inv[f]
    root=inv[np.asarray(root_ids,int)];rootset=set(map(int,root))
    if len(root)!=6 or len(rootset)!=6 or len(ids)!=4 or len(set(map(int,ids)))!=4:
        raise ValueError('An explicit six-point shared profile and four source disk triangles required')
    if set(ff[ids].ravel())!=rootset:raise ValueError('Disk controls differ from the declared source profile')
    audit=continuity_audit(p,f[ids],expected_components=1,expected_mouths=1,expected_euler=1)
    if not audit['topological_continuity_pass']:raise ValueError('Declared source surface is not one oriented disk')
    loops=boundary_loops(ff[ids]);expected={tuple(sorted(map(int,[a,b])))for a,b in zip(root,np.roll(root,-1))}
    actual={tuple(sorted(map(int,[a,b])))for a,b in zip(loops[0],np.roll(loops[0],-1))}
    if actual!=expected:raise ValueError('Disk boundary is not the ordered shared source profile')
    owners=defaultdict(list)
    for fi,t in enumerate(ff):
        for a,b in zip(t,np.roll(t,-1)):owners[tuple(sorted(map(int,[a,b])))].append(fi)
    sides=[];caps=set(map(int,ids))
    for a,b in zip(root,np.roll(root,-1)):
        incident=owners[tuple(sorted(map(int,[a,b])))];disk=[i for i in incident if i in caps];walls=[i for i in incident if i not in caps]
        if len(disk)!=1 or len(walls)!=2:raise ValueError('Shared profile must have one disk and two existing wall faces on every edge')
        t=physical[ff[disk[0]]];normal=np.cross(t[1]-t[0],t[2]-t[0]);normal/=np.linalg.norm(normal)
        third=[next(int(v)for v in ff[i]if v not in [a,b])for i in walls];signed=(physical[third]-physical[a])@normal
        if signed[0]*signed[1]>=0 or np.abs(signed).min()<minimum_wall_side_m:
            raise ValueError('Existing source walls do not lie strictly on opposite sides of this disk')
        sides.append(dict(rootEdgeRawIds=rep[[a,b]].tolist(),diskFace=disk[0],wallFaces=walls,signedOppositeWallDistances_m=signed.tolist()))
    keep=np.setdiff1d(np.arange(len(f)),ids);after=ff[keep]
    for center in root:
        star=after[np.any(after==center,axis=1)];graph=defaultdict(set)
        for tri in star:
            a,b=[int(v)for v in tri if v!=center];graph[a].add(b);graph[b].add(a)
        todo=[min(graph)];seen=set()
        while todo:
            v=todo.pop()
            if v not in seen:seen.add(v);todo.extend(graph[v]-seen)
        if len(seen)!=len(graph)or any(len(x)!=2 for x in graph.values()):
            raise ValueError('Disk removal does not leave a single manifold link at every profile point')
    count=lambda faces:Counter(tuple(sorted(map(int,[a,b])))for t in faces for a,b in zip(t,np.roll(t,-1)))
    old,new=count(ff),count(after)
    if {e for e,n in old.items()if n==1}!={e for e,n in new.items()if n==1}:
        raise ValueError('An external boundary changed during shared disk removal')
    if sum(n>2 for n in old.values())-sum(n>2 for n in new.values())!=6:
        raise ValueError('Expected exactly six repaired nonmanifold profile edges')
    return keep,dict(complete=True,rootRawIds=list(map(int,root_ids)),removedSourceFaceIds=ids.tolist(),retainedSourceFaceIds=keep.tolist(),
        diskAudit=audit,oppositeSourceWallEvidence=sides,minimumAbsoluteOppositeWallDistance_m=min(abs(v)for row in sides for v in row['signedOppositeWallDistances_m']),
        allProfileVertexLinksManifold=True,externalBoundaryEdgesLiteral=True,sourceCoordinatesChanged=False,
        scope='Caller-declared source disk; no anatomical valve identity inferred',anatomicalAcceptance=False)
