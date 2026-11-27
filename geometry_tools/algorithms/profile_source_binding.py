"""Rebind source section ownership across literal-control diagonal repairs.

This permits only two-triangle, four-control patches with identical oriented
boundaries. It never authorizes moving source controls, adding/removing faces,
welding points, changing an interface or accepting the resulting anatomy.
"""
from collections import Counter
import numpy as np


def _boundary(faces):
    edges=Counter((int(a),int(b))for f in faces for a,b in zip(f,np.roll(f,-1)))
    return Counter({e:n-edges[e[::-1]]for e,n in edges.items()if n>edges[e[::-1]]})


def verify_profile_topology_replay(before,after):
    a,b=np.asarray(before,int),np.asarray(after,int)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3:
        raise ValueError('Literal corresponding source face count required')
    changed=set(map(int,np.flatnonzero(np.any(a!=b,axis=1))));patches=[]
    while changed:
        i=min(changed);match=None
        for j in sorted(changed-{i}):
            ids=[i,j];old,new=a[ids],b[ids]
            vertices=np.unique(old)
            if len(vertices)!=4 or not np.array_equal(vertices,np.unique(new)):continue
            if len(_boundary(old))!=4 or _boundary(old)!=_boundary(new):continue
            match=j;patches.append(dict(faceIds=ids,oldFaces=old.tolist(),newFaces=new.tolist()));break
        if match is None:raise ValueError('Source profile cache differs beyond explicit oriented quadrilateral diagonal patches')
        changed.difference_update([i,match])
    return dict(complete=True,changedPatches=patches,changedFaces=2*len(patches),sourceControlIdentityMustBeVerifiedSeparately=True,sourcePatchBoundariesLiteral=True,anatomicalAcceptance=False)
