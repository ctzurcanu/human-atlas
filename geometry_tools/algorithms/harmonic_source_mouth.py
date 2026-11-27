"""Align one explicit existing mouth, with finite graph support and fixed ends.

No nearest-point correspondence is inferred. Source edges define a positive
Dirichlet displacement field; all other complete source openings stay fixed.
This is a proposal only. Facet strain, source profiles, crossings, recipients
and packed coordinates require separate caller checks.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.csgraph import dijkstra
from scipy.sparse.linalg import spsolve
from geometry_tools.model import Mesh
from .tube_continuity import boundary_loops
from .surface_normals import transport_vertex_words

def align_source_mouth(mesh,raw_ids,target,*,support_m,maximum_motion_m,
                       maximum_mouth_area_change,maximum_facet_rotation_degrees):
    raw=np.asarray(raw_ids,int);target=np.asarray(target,float)
    if raw.ndim!=1 or len(raw)<3 or target.shape!=(len(raw),3)or not np.isfinite(target).all():raise ValueError('Complete corresponding source and target mouth required')
    if not 0<support_m<=.02 or not 0<maximum_motion_m<=.001:raise ValueError('Finite local support and motion budgets required')
    if raw.min()<0 or raw.max()>=len(mesh.vertices):raise ValueError('Mouth source index outside geometry')
    p,rep,inv=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True);f=inv[mesh.faces];root=inv[raw]
    loops=boundary_loops(f);matches=[l for l in loops if set(l)==set(root)]
    if len(np.unique(root))!=len(root)or len(matches)!=1:raise ValueError('Declared source controls are not one complete actual opening')
    area=lambda q:float(np.linalg.norm(np.cross(q-q.mean(0),np.roll(q,-1,axis=0)-q.mean(0)).sum(0))*.5)
    change=abs(area(target)/area(mesh.vertices[raw])-1)
    if change>maximum_mouth_area_change:raise ValueError('Target mouth exceeds source area budget')
    edges=np.unique(np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0);length=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)
    if np.any(length<=0):raise ValueError('Degenerate physical edge')
    a,b=edges.T;distance_graph=coo_matrix((np.r_[length,length],(np.r_[a,b],np.r_[b,a])),shape=(len(p),len(p))).tocsr()
    distances=dijkstra(distance_graph,directed=False,indices=root,min_only=True)
    other=np.unique(np.concatenate([l for l in loops if set(l)!=set(root)]))if len(loops)>1 else np.empty(0,int)
    fixed=np.unique(np.r_[root,other,np.flatnonzero(distances>=support_m)]);free=np.setdiff1d(np.arange(len(p)),fixed)
    displacement=np.zeros_like(p);displacement[root]=target-p[root]
    weights=1/length;adj=coo_matrix((np.r_[weights,weights],(np.r_[a,b],np.r_[b,a])),shape=(len(p),len(p))).tocsr();L=diags(np.asarray(adj.sum(1)).ravel())-adj
    if len(free):displacement[free]=spsolve(L[free][:,free].tocsc(),-L[free][:,fixed]@displacement[fixed])
    if not np.isfinite(displacement).all():raise ValueError('Unanchored or invalid displacement field')
    q=mesh.vertices+displacement[inv];q[np.isin(inv,root)]=target[np.array([np.flatnonzero(root==i)[0]for i in inv[np.isin(inv,root)]])]
    motion=float(np.linalg.norm(q-mesh.vertices,axis=1).max())
    if motion>maximum_motion_m:raise ValueError('Source-mouth motion budget exceeded')
    np.testing.assert_array_equal(q[np.isin(inv,other)],mesh.vertices[np.isin(inv,other)])
    attrs={k:v.copy()for k,v in mesh.attributes.items()};attrs['normals'],normals=transport_vertex_words(mesh.vertices,q,mesh.faces,attrs['normals'],maximum_facet_rotation_degrees=maximum_facet_rotation_degrees)
    return Mesh(q,mesh.faces,mesh.frame,attrs,mesh.source),dict(complete=True,method='positive source-edge Dirichlet displacement with explicit finite support',sourceRawIds=raw.tolist(),targetPoints_m=target.tolist(),support_m=support_m,maximumMotion_m=motion,relativeMouthAreaChange=change,fixedOtherMouths=len(loops)-1,sourceFacesRetained=True,exactCoordinateAliasesCoupled=True,sourceNormalTransport=normals,changedRawIds=np.flatnonzero(np.any(q!=mesh.vertices,axis=1)).tolist(),anatomicalAcceptance=False)
