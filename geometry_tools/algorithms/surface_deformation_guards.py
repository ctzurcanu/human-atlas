"""Exact source-tangent displacement-gradient bounds for cell motion.

For each triangle, ||D d||_F < 1 keeps its tangent map noncollapsed and its
projection orientation positive. This is a local guarantee only; global
self-contact, tissue packing and anatomical appropriateness remain separate.
"""
import numpy as np


def surface_deformation_guards(vertices,faces,cells,movable_cells,*,maximum_norm=.9):
    v=np.asarray(vertices,float);f=np.asarray(faces,int);cell=np.asarray(cells,int)
    if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or cell.shape!=(len(v),)or not np.isfinite(v).all()or np.any(f<0)or np.any(f>=len(v))or not 0<maximum_norm<1:
        raise ValueError('Finite source triangles, source cells and a strict local norm bound below one required')
    ids=cell[f];moving=np.isin(ids,np.array(sorted(set(map(int,movable_cells))),int));keep=np.any(moving,axis=1)&np.any(ids!=ids[:,[0]],axis=1);source_faces=np.flatnonzero(keep);ids=ids[keep];p=v[f[keep]]
    e1=p[:,1]-p[:,0];e2=p[:,2]-p[:,0];length=np.linalg.norm(e1,axis=1);height=np.linalg.norm(np.cross(e1,e2),axis=1)/length
    if not np.isfinite(height).all()or np.any(length<=0)or np.any(height<=0):raise ValueError('Degenerate source faces must be repaired before gradient construction')
    projection=np.einsum('ij,ij->i',e2,e1/length[:,None]);weights=np.zeros((len(ids),3,2));weights[:,1,0]=1/length;weights[:,0,0]=-1/length;weights[:,1,1]=-projection/(length*height);weights[:,2,1]=1/height;weights[:,0,1]=-weights[:,1:,1].sum(1)
    return dict(cell_ids=ids,gradient_weights=weights,source_face_ids=source_faces,maximum_norm=np.full(len(ids),maximum_norm),minimum_source_height_m=np.array(float(height.min())if len(height)else 0.))


def evaluate_surface_guards(delta,guards):
    d=np.asarray(delta,float);ids=np.asarray(guards['cell_ids'],int);w=np.asarray(guards['gradient_weights'],float);gradient=np.einsum('nka,nkb->nab',d[ids][:,:,[0,2]],w);norm=np.linalg.norm(gradient,axis=(1,2));limits=np.asarray(guards['maximum_norm'],float)
    return dict(maximumNorm=float(norm.max(initial=0.)),maximumViolation=float((norm-limits).max(initial=0.)),violatingFaces=int(np.count_nonzero(norm>limits+1e-8)),norms=norm)
