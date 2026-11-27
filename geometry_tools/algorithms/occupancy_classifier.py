"""Analysis-only closing of literal open vascular rims for extent-volume queries.

Classifier caps are never anatomical surfaces or delivered lumen plugs. The
result describes an auxiliary occupied outer extent, not finite wall tissue.
"""
from collections import defaultdict
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology


def close_literal_rims(mesh):
    vertices,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True)
    faces=inverse[mesh.faces];owners=defaultdict(list)
    for face in faces:
        for a,b in zip(face,np.roll(face,-1)):owners[tuple(sorted((int(a),int(b))))].append((int(a),int(b)))
    if any(len(x)>2 for x in owners.values()):raise ValueError('Nonmanifold source cannot define an occupancy classifier')
    boundary=[x[0]for x in owners.values()if len(x)==1];outgoing=defaultdict(list);incoming=defaultdict(list)
    for a,b in boundary:outgoing[a].append(b);incoming[b].append(a)
    if any(len(outgoing[v])!=1 or len(incoming[v])!=1 for v in set(outgoing)|set(incoming)):
        raise ValueError('Source rim is branched or inconsistently oriented')
    remaining=set(outgoing);v=vertices.tolist();caps=[];loops=[]
    while remaining:
        first=min(remaining);ring=[];current=first
        while current not in ring:
            ring.append(current);remaining.discard(current);current=outgoing[current][0]
        if current!=first or len(ring)<3:raise ValueError('Source rim does not form a simple oriented cycle')
        points=vertices[ring];center=points.mean(0);ci=len(v);v.append(center.tolist())
        # A centre fan is permitted only if it lies consistently inside the
        # projected convex rim. Nonconvex/degenerate source mouths are explicit.
        normal=np.cross(points-center,np.roll(points,-1,axis=0)-center);total=normal.sum(0)
        if np.linalg.norm(total)==0 or np.any(normal@total<=0):raise ValueError('Nonconvex or degenerate rim needs an explicit classifier construction')
        for a,b in zip(ring,ring[1:]+ring[:1]):caps.append([b,a,ci])
        loops.append(dict(sourceRawPhysicalIds=ring,center_m=center.tolist(),points=len(ring),bounds=[points.min(0).tolist(),points.max(0).tolist()]))
    result=Mesh(np.asarray(v),np.r_[faces,np.asarray(caps,int).reshape(-1,3)],mesh.frame)
    audit=topology(result)
    if not audit['closed_oriented_manifold']:raise ValueError('Auxiliary classifier is not a closed oriented manifold')
    if audit['oriented_triangle_integral_m3']<=0:raise ValueError('Auxiliary classifier does not have positive source orientation')
    return result,dict(analysisOnly=True,sourceGeometrySha256=mesh.geometry_sha,sourceFacesRetained=True,
        sourceFaces=len(faces),classifierFaces=len(caps),literalRims=loops,closed_oriented_manifold=True,
        finiteHistologicalWallEstablished=False,deliveryPermitted=False)
