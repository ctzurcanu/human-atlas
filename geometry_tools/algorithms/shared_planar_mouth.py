"""Join two construction envelopes through the same genuine planar mouth.

Only coincident planar closure faces are removed. Remaining boundaries must
agree geometrically; collinear subdivisions propagate to both sides. This
avoids zero-volume fragments from a Boolean between coplanar closure disks.
"""
from collections import defaultdict
import numpy as np
from .finite_tissue_boolean import CommonPositionBank
from .shared_triangle_edges import reconcile_boundaries,refine_surfaces,boundary

def stitch_planar_mouth(vertices_a,faces_a,vertices_b,faces_b,normal,offset,tolerance=2e-12):
    n=np.asarray(normal,float);n/=np.linalg.norm(n);bank=CommonPositionBank(np.asarray(vertices_a,float),tolerance)
    ids=np.array([bank.insert(p)for p in vertices_b],int);surfaces={};removed=[];source_faces={}
    for name,v,f in [('a',np.asarray(vertices_a),np.asarray(faces_a,int)),('b',np.asarray(vertices_b),np.asarray(faces_b,int))]:
        coplanar=np.all(abs(v[f]@n-offset)<=tolerance,axis=1)
        if not coplanar.any():raise ValueError('Both envelopes need a genuine planar closure at this mouth')
        side=f[~coplanar];source_faces[name]=np.flatnonzero(~coplanar);surfaces[name]=side if name=='a'else ids[side];removed.append(int(coplanar.sum()))
    plans=defaultdict(set);residual=reconcile_boundaries(bank.positions,surfaces['a'],surfaces['b'][:,[0,2,1]],plans,tolerance)
    surfaces,parents,count=refine_surfaces(bank,surfaces,plans);joined=np.r_[surfaces['a'],surfaces['b']]
    if boundary(joined):raise ValueError('The joined envelope is not closed')
    return bank.positions,joined,dict(removedConstructionDiskTriangles=removed,refinedSideFaces=count,maximumBoundaryResidualMeters=residual,sourceAParents=source_faces['a'][parents['a']],sourceBParents=source_faces['b'][parents['b']])
