"""Recover an actual open tube's ordered rings from its mesh incidence.

Two real boundary cycles seed graph-distance layers. Every layer must itself
be a cycle of the same size, and every original face must connect adjacent
layers. No assumed quad pairing, thinning, curve invention or mesh replacement.
The resulting course is a guide, not proof of anatomical placement.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from geometry_tools.model import Mesh
from geometry_tools.oriented_rays import boundary_loops


def annular_mesh_course(mesh,root_direction):
    direction=np.asarray(root_direction,float)
    if direction.shape!=(3,)or not np.isfinite(direction).all()or np.linalg.norm(direction)<1e-12:
        raise ValueError('An explicit nonzero source-root direction required')
    used=np.unique(mesh.faces);v,representative,inv=np.unique(mesh.vertices[used],axis=0,return_index=True,return_inverse=True)
    raw_to_local=np.full(len(mesh.vertices),-1,int);raw_to_local[used]=inv;f=raw_to_local[mesh.faces]
    _,loops,_,topology=boundary_loops(Mesh(v,f,mesh.frame))
    if len(loops)!=2 or topology['nonmanifold_edges']or topology['inconsistent_shared_edge_winding']:
        raise ValueError('Exactly two actual manifold source mouths required')
    if len(loops[0])!=len(loops[1])or len(loops[0])<3:raise ValueError('Corresponding source boundary cycles required')
    root=max(loops,key=lambda s:float(v[s].mean(0)@direction));size=len(root)
    edge=np.unique(np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
    graph=coo_matrix((np.ones(2*len(edge)),(np.r_[edge[:,0],edge[:,1]],np.r_[edge[:,1],edge[:,0]])),shape=(len(v),len(v))).tocsr()
    distance=dijkstra(graph,directed=False,indices=root,unweighted=True,min_only=True)
    if not np.isfinite(distance).all()or not np.array_equal(distance,np.rint(distance)):
        raise ValueError('A single integral source-incidence tube required')
    layer=distance.astype(int);rings=[np.flatnonzero(layer==k)for k in range(int(layer.max())+1)]
    if any(len(s)!=size for s in rings)or len(f)!=2*size*(len(rings)-1):
        raise ValueError('Source does not have complete equal-sized annular layers')
    if np.any(np.ptp(layer[f],axis=1)!=1):raise ValueError('A source face skips or lies wholly in a course layer')
    for s in rings:
        local=graph[s][:,s]
        if not np.all(np.diff(local.indptr)==2):raise ValueError('A source layer is not a simple boundary cycle')
    centres=np.array([v[s].mean(0)for s in rings]);radius=np.array([np.linalg.norm(v[s]-c,axis=1).mean()for s,c in zip(rings,centres)])
    return dict(centres=centres,radius_m=radius,ring_representatives=np.array([used[representative[s]]for s in rings])),dict(
        method='Actual source mouths and complete graph-distance annuli',rings=len(rings),controlsPerRing=size,
        everySourceFaceAccountedFor=True,sourceGeometrySha256=mesh.geometry_sha,rootDirection=direction.tolist(),
        sourceFacesReplaced=False,anatomicalAcceptance=False)
