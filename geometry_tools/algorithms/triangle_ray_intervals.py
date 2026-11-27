"""Bounded-memory first ray/surface hit on the actual Float64 triangles.

No distant body silhouette, hole filling, vertex welding or broad-phase proxy
is substituted for a local surface. A miss remains explicit. This is floating
point evaluation, not exact arithmetic or a closed-volume proof.
"""
import numpy as np

def first_positive_triangle_hits(origins,directions,triangles,*,maximum_distance_m,batch_size=8):
 o=np.asarray(origins,float);d=np.asarray(directions,float);t=np.asarray(triangles,float)
 if o.ndim!=2 or o.shape[1]!=3 or d.shape!=o.shape or t.ndim!=3 or t.shape[1:]!=(3,3)or not np.isfinite(np.r_[o.ravel(),d.ravel(),t.ravel()]).all():raise ValueError('Finite 3D rays and actual triangles required')
 length=np.linalg.norm(d,axis=1)
 if np.any(length==0)or not np.isfinite(maximum_distance_m)or maximum_distance_m<=0 or batch_size<1:raise ValueError('Nonzero rays, positive finite travel cap and batch size required')
 d=d/length[:,None];a=t[:,0];e=t[:,1]-a;f=t[:,2]-a;distance=np.full(len(o),np.nan);face=np.full(len(o),-1,int);points=np.full(o.shape,np.nan);weights=np.full(o.shape,np.nan)
 for start in range(0,len(o),batch_size):
  q=o[start:start+batch_size];direction=d[start:start+batch_size];p=np.cross(direction[:,None],f);det=np.einsum('fi,bfi->bf',e,p);valid=det!=0;inv=np.divide(1.,det,out=np.zeros_like(det),where=valid);relative=q[:,None]-a;u=np.einsum('bfi,bfi->bf',relative,p)*inv;cross=np.cross(relative,e);v=np.einsum('bi,bfi->bf',direction,cross)*inv;travel=np.einsum('fi,bfi->bf',f,cross)*inv;valid&=(u>=0)&(v>=0)&(u+v<=1)&(travel>0)&(travel<=maximum_distance_m);travel[~valid]=np.inf;j=travel.argmin(1);row=np.arange(len(q));found=np.isfinite(travel[row,j]);target=start+row[found];jj=j[found];tt=travel[row[found],jj];distance[target]=tt;face[target]=jj;weights[target]=np.c_[1-u[row[found],jj]-v[row[found],jj],u[row[found],jj],v[row[found],jj]];points[target]=np.einsum('ni,nij->nj',weights[target],t[jj])
 hit=face>=0;residual=np.full(len(o),np.nan);residual[hit]=np.linalg.norm(points[hit]-o[hit]-distance[hit,None]*d[hit],axis=1)
 return dict(hit=hit,face_ids=face,distance_m=distance,points=points,barycentric=weights,ray_residual_m=residual,all_actual_triangles_tested=True,maximum_distance_m=float(maximum_distance_m),closed_volume_proven=False)
