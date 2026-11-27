"""One linear source-edge field for all explicitly paired mouths of an owner."""
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.csgraph import dijkstra
from scipy.sparse.linalg import spsolve
from .tube_continuity import boundary_loops

def mouth_displacement_basis(mesh,declarations,*,support_m):
    if not 0<support_m<=.02 or not declarations:raise ValueError('Explicit finite support and source-mouth declarations required')
    p,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces];loops=boundary_loops(f);root=set();values={}
    for d in declarations:
        raw=np.asarray(d['sourceRawIds'],int);target=np.asarray(d['targetPoints_m'],float)
        if raw.ndim!=1 or len(raw)<3 or raw.min()<0 or raw.max()>=len(mesh.vertices)or target.shape!=(len(raw),3)or not np.isfinite(target).all():raise ValueError('Complete corresponding source-mouth controls required')
        ids=inv[raw]
        if len(np.unique(ids))!=len(ids)or len([loop for loop in loops if set(loop)==set(ids)])!=1:raise ValueError('A declared root is not a complete actual source mouth')
        if root.intersection(ids):raise ValueError('Declared root curves share physical ownership')
        root.update(map(int,ids))
        for i,q in zip(ids,target):values[int(i)]=q-p[i]
    roots=np.array(sorted(root));edges=np.unique(np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0);length=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)
    if np.any(length<=0):raise ValueError('Degenerate source edge')
    a,b=edges.T;graph=coo_matrix((np.r_[length,length],(np.r_[a,b],np.r_[b,a])),shape=(len(p),len(p))).tocsr();distance=dijkstra(graph,directed=False,indices=roots,min_only=True)
    other=np.array(sorted(set(v for loop in loops for v in loop)-root),dtype=int);fixed=np.unique(np.r_[roots,other,np.flatnonzero(distance>=support_m)]).astype(int);free=np.setdiff1d(np.arange(len(p)),fixed)
    field=np.zeros_like(p);field[roots]=np.array([values[int(i)]for i in roots]);w=1/length;adj=coo_matrix((np.r_[w,w],(np.r_[a,b],np.r_[b,a])),shape=(len(p),len(p))).tocsr();L=diags(np.asarray(adj.sum(1)).ravel())-adj
    if len(free):field[free]=spsolve(L[free][:,free].tocsc(),-L[free][:,fixed]@field[fixed])
    if not np.isfinite(field).all():raise ValueError('Invalid source Dirichlet field')
    np.testing.assert_array_equal(field[other],np.zeros((len(other),3)))
    return field[inv],dict(complete=True,support_m=support_m,rootPhysicalIds=roots.tolist(),fixedOtherBoundaryPhysicalIds=other.tolist(),freePhysicalControls=len(free),sourceAdjacencyLiteral=True,sourceGraphWeightsFixed=True,fieldIsLinearInDeclaredTargetDisplacement=True,anatomicalAcceptance=False)
