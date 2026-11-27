"""Optional native winding-number queries for explicitly oriented triangle soups.

Uses libigl's implementation of Barill et al., Fast Winding Numbers for Soups
and Clouds (2018). An inside classification is a numerical reconstruction
choice, not proof of watertight topology or anatomical tissue identity.
"""
import numpy as np

def winding_numbers(vertices,faces,points):
    import igl
    v=np.asarray(vertices,dtype='<f8');f=np.asarray(faces,dtype='<i8');q=np.asarray(points,dtype='<f8')
    if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or q.ndim!=2 or q.shape[1]!=3:raise ValueError('Explicit 3D vertices, faces and query points required')
    if not np.isfinite(v).all()or not np.isfinite(q).all()or not len(f)or f.min()<0 or f.max()>=len(v):raise ValueError('Invalid winding-number source or queries')
    origin=v.mean(0)
    return igl.fast_winding_number(v-origin,f,q-origin)
