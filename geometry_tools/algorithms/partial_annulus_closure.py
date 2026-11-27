"""Complete one explicitly declared five-quad annulus for computation only.

The added source-gap faces must be independently forbidden in emitted anatomy.
This helper alone never identifies an anatomical opening or accepts a surface.
"""
from collections import Counter
import numpy as np
from .tube_continuity import boundary_loops,continuity_audit


def complete_five_quad_collar(vertices,faces,cell_faces,root_ids):
    p=np.asarray(vertices,float);f=np.asarray(faces,int);ids=np.asarray(cell_faces,int)
    _,rep,inv=np.unique(p,axis=0,return_index=True,return_inverse=True)
    root=set(map(int,inv[np.asarray(root_ids,int)]));ff=f[ids];physical=inv[ff]
    if len(root)!=6 or len(ids)!=10 or len(np.unique(physical))!=12:
        raise ValueError('One explicit six-control root and ten-face twelve-control partial collar required')
    loops=boundary_loops(physical)
    if len(loops)!=1 or len(loops[0])!=12 or not root<=set(physical.ravel()):
        raise ValueError('Expected one literal twelve-edge source gap boundary')
    directed=Counter((int(a),int(b))for tri in physical for a,b in zip(tri,np.roll(tri,-1)))
    rim=[e for e,n in directed.items()if n and directed[e[::-1]]==0]
    cross=[e for e in rim if (e[0]in root)!=(e[1]in root)]
    if len(cross)!=2:raise ValueError('Exactly two source gap side edges required')
    first=next((e for e in cross if e[0]in root),None)
    second=next((e for e in cross if e[1]in root),None)
    if first is None or second is None:raise ValueError('Source gap orientation is inconsistent')
    ra,oa=first;ob,rb=second
    root_cycle=inv[np.asarray(root_ids,int)]
    root_edges={tuple(sorted(map(int,[a,b])))for a,b in zip(root_cycle,np.roll(root_cycle,-1))}
    if tuple(sorted((ra,rb)))not in root_edges:
        raise ValueError('The gap endpoints do not form an edge of the declared root')
    raw={}
    for v in ff.ravel():raw.setdefault(int(inv[v]),int(v))
    added=np.array([[raw[oa],raw[ra],raw[rb]],[raw[oa],raw[rb],raw[ob]]],int)
    augmented=np.r_[ff,added]
    audit=continuity_audit(p,augmented,expected_components=1,expected_mouths=2,expected_euler=0)
    if not audit['topological_continuity_pass']:raise ValueError('Computational completion is not a valid annulus')
    completed_loops=boundary_loops(inv[augmented])
    if sorted(map(len,completed_loops))!=[6,6]or not any(set(l)==root for l in completed_loops):
        raise ValueError('Completed collar did not retain the literal six-control root and outer port')
    return added,dict(complete=True,sourceCellFaceIds=ids.tolist(),addedGapTriangles=added.tolist(),
        rootRawIds=list(map(int,root_ids)),sourcePointsMoved=0,computationOnly=True,
        emittedGapFacesPermitted=False,annulusAudit=audit,anatomicalAcceptance=False)
