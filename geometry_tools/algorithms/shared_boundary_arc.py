"""Reconcile an explicitly identified split boundary arc of two tube sheets.

Two annular source sheets already meet at the arc endpoints. Align only the
intermediate near-coincident controls. The union retains three real mouths;
there is no generated bridge, cap, replacement branch or coarse remeshing.
"""
import numpy as np
from .tube_continuity import boundary_loops
from ..model import Mesh


def topology(vertices,faces):
    used=np.unique(faces);v,inverse=np.unique(np.asarray(vertices)[used],axis=0,return_inverse=True)
    lookup=np.full(len(vertices),-1,int);lookup[used]=inverse;f=lookup[faces]
    edges=np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1)
    e,edge_ids,count=np.unique(edges,axis=0,return_inverse=True,return_counts=True)
    parent=np.arange(len(f))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=int(parent[i])
        return i
    first={}
    for i,k in enumerate(edge_ids):
        face=i//3
        if int(k)in first:parent[find(face)]=find(first[int(k)])
        else:first[int(k)]=face
    rings=boundary_loops(f)
    return dict(physicalVertices=len(v),edges=len(e),faces=len(f),eulerCharacteristic=len(v)-len(e)+len(f),
                boundaryLoops=len(rings),boundaryLoopLengths=sorted(map(len,rings)),nonmanifoldEdges=int((count>2).sum()),edgeConnectedComponents=len({find(i)for i in range(len(f))}))


def reconcile_shared_arc(mesh,first_faces,second_faces,first_ring,second_ring,*,maximum_motion_m):
    if not 0<maximum_motion_m<=5e-7:raise ValueError('An explicit source-rounding bound no larger than0.5micrometres is required')
    subsets=[np.asarray(x,int)for x in [first_faces,second_faces]]
    if not all(len(x)>0 for x in subsets)or not np.array_equal(np.sort(np.r_[subsets[0],subsets[1]]),np.arange(len(mesh.faces))):
        raise ValueError('Both complete disjoint source sheets must partition the actual carrier')
    v,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inverse[mesh.faces];rings=[];external_sizes=[]
    for ids,ring in zip(subsets,[first_ring,second_ring]):
        target=inverse[np.asarray(ring,int)];actual=boundary_loops(f[ids])
        if len(actual)!=2:raise ValueError('Each source sheet must retain two complete external mouths')
        candidates=[x for x in actual if set(x)==set(target)]
        if len(candidates)!=1 or len(target)!=6 or len(np.unique(target))!=6:raise ValueError('An explicit complete six-point source cycle is required')
        rings.append(np.asarray(candidates[0],int))
        external_sizes.append(len(next(x for x in actual if set(x)!=set(target))))
        source_topology=topology(mesh.vertices,mesh.faces[ids])
        if source_topology['eulerCharacteristic']!=0 or source_topology['edgeConnectedComponents']!=1:raise ValueError('Each complete source sheet must be one annulus without an extra handle or component')
    a,b=rings;distance=np.linalg.norm(v[a,None]-v[b][None],axis=2);matches={}
    for i in range(len(a)):
        close=np.flatnonzero(distance[i]<=maximum_motion_m)
        if len(close)>1:raise ValueError('Ambiguous source arc correspondence')
        if len(close):matches[i]=int(close[0])
    if len(matches)!=4 or len(set(matches.values()))!=4:raise ValueError('Exactly four uniquely corresponding source arc controls are required')
    exact=[i for i,j in matches.items()if np.array_equal(v[a[i]],v[b[j]])]
    if len(exact)!=2:raise ValueError('The two source arc endpoints must already be exactly shared')
    edges=[]
    for i in range(6):
        j=(i+1)%6
        if i in matches and j in matches:
            if (matches[j]+1)%6!=matches[i]:raise ValueError('Shared boundary arc winding must be opposite')
            edges.append((i,j))
    degree={i:sum(i in e for e in edges)for i in matches}
    if len(edges)!=3 or sorted(degree.values())!=[1,1,2,2]or {i for i,d in degree.items()if d==1}!=set(exact):
        raise ValueError('Correspondences must be one three-edge arc with the original shared endpoints')
    q=mesh.vertices.copy();moves=[]
    for i,j in matches.items():
        old,target=v[b[j]],v[a[i]];ids=np.flatnonzero(np.all(mesh.vertices==old,axis=1))
        if not np.array_equal(old,target):
            if np.any(np.isin(mesh.faces[subsets[0]],ids)):raise ValueError('A moving control belongs to the fixed source sheet')
            q[ids]=target;moves.append(dict(rawIds=ids.tolist(),from_m=old.tolist(),to_m=target.tolist(),distance_m=float(np.linalg.norm(old-target))))
    if len(moves)!=2:raise ValueError('Only two intermediate source controls may move')
    result=Mesh(q,mesh.faces,mesh.frame,mesh.attributes,mesh.source);proof=topology(q,mesh.faces)
    if proof['boundaryLoops']!=3 or proof['boundaryLoopLengths']!=sorted(external_sizes+[6])or proof['eulerCharacteristic']!=-1 or proof['nonmanifoldEdges']or proof['edgeConnectedComponents']!=1:
        raise ValueError('The result must retain the three original junction mouths, without a cap or extra handle')
    return result,dict(moves=moves,maximumMotion_m=max(x['distance_m']for x in moves),sharedArcEdges=3,
        fixedSourceSheetLiteral=True,facesLiteral=True,sourceControlCountUnchanged=True,
        topology=proof,newFaces=0,newCaps=0,anatomicalAcceptance=False)
