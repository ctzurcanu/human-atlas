"""Normalize a shared boundary complex on its actual common Float32 bank.

Only literal coordinate aliases, exactly zero-area facets and reciprocal
calculation faces cancel. There is no distance weld or independent per-cell
geometry repair. Full emitted-cell and occupied-volume audits are mandatory.
"""
import numpy as np

def normalize_shared_boundary(positions,cells,origin=(0.,0.,0.)):
 p=np.asarray(positions,float);o=np.asarray(origin,float)
 if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all()or o.shape!=(3,):raise ValueError('Finite shared positions and origin required')
 relative=(p-o).astype('<f4');v,remap=np.unique(relative,axis=0,return_inverse=True);v=v.astype(float)+o;result={};counts={}
 for name,source in cells.items():
  f=remap[np.asarray(source,int)];q=v[f];area=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1);zero=area==0;f=f[~zero]
  if len(f):
   # Accumulate the complete oriented boundary chain before rejecting a
   # duplicate. +,+,- has one surviving face, independently of input order.
   first=np.argmin(f,axis=1);cyclic=np.take_along_axis(f,(first[:,None]+np.arange(3))%3,axis=1);sign=np.where(cyclic[:,1]<cyclic[:,2],1,-1)
   canonical,inv=np.unique(np.sort(f,axis=1),axis=0,return_inverse=True);coefficient=np.bincount(inv,weights=sign,minlength=len(canonical)).astype('int64')
   if np.any(abs(coefficient)>1):raise ValueError('Identically oriented duplicate tissue facet after reciprocal cancellation')
   retained=coefficient!=0;emitted=canonical[retained].copy();reverse=coefficient[retained]<0;emitted[reverse]=emitted[reverse][:,[0,2,1]]
  else:emitted=np.empty((0,3),dtype='int64')
  cancelled=len(f)-len(emitted);result[name]=emitted;counts[name]=dict(source_triangles=len(source),exact_zero_area_removed=int(zero.sum()),reciprocal_calculation_triangles_cancelled=cancelled,emitted_triangles=len(result[name]))
 return dict(positions=v,cells=result,source_position_map=remap,evidence=dict(common_float32_origin=o.tolist(),source_positions=len(p),stored_positions=len(v),literal_aliases_welded=len(p)-len(v),maximum_position_error_m=float(np.linalg.norm(v[remap]-p,axis=1).max()),cells=counts,no_distance_welding=True,anatomical_acceptance=False))
