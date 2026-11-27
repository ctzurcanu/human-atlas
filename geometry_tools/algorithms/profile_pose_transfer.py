"""Transfer section motion through an explicit source correspondence.

An orthogonal direction map may reverse handedness, but conjugating each
proper section rotation preserves recipient handedness and section shape.
The caller must establish control correspondence from source provenance;
this function never matches or welds nearby points.
"""
import numpy as np


def transfer_profile_poses(before,after,recipient,profiles,direction_map,*,maximum_source_residual_m):
 a,b,c=[np.asarray(x,float)for x in [before,after,recipient]];M=np.asarray(direction_map,float)
 if a.ndim!=2 or a.shape[1]!=3 or a.shape!=b.shape or a.shape!=c.shape or not np.isfinite(np.r_[a.ravel(),b.ravel(),c.ravel()]).all():raise ValueError('Finite explicitly corresponding source and recipient controls required')
 if M.shape!=(3,3)or not np.isfinite(M).all()or not np.allclose(M.T@M,np.eye(3),atol=1e-12,rtol=0):raise ValueError('Explicit orthogonal direction map required')
 if not np.isfinite(maximum_source_residual_m)or not 0<=maximum_source_residual_m<=1e-5:raise ValueError('Explicit finite source pose residual tolerance required')
 centres=[];rotations=[];errors=[];owned=np.zeros(len(a),bool)
 for raw in profiles:
  p=np.asarray(raw,int)
  if p.ndim!=1 or len(p)<3 or np.any(p<0)or np.any(p>=len(a))or len(np.unique(p))!=len(p)or owned[p].any():raise ValueError('Disjoint corresponding source section controls required')
  owned[p]=True;ac,bc,cc=[x[p].mean(0)for x in [a,b,c]]
  A=a[p]-ac;B=b[p]-bc
  if np.linalg.matrix_rank(A,tol=1e-12)<2:raise ValueError('Noncollinear section controls required')
  u,_,vt=np.linalg.svd(B.T@A);D=np.eye(3);D[-1,-1]=np.linalg.det(u@vt);R=u@D@vt
  error=float(np.linalg.norm(A@R.T-B,axis=1).max())
  if error>maximum_source_residual_m:raise ValueError('Source section motion is not rigid within the declared tolerance')
  centres.append(cc+M@(bc-ac));rotations.append(M@R@M.T);errors.append(error)
 return dict(target_centres=np.asarray(centres).reshape(-1,3),rotations=np.asarray(rotations).reshape(-1,3,3),source_pose_residual_m=np.asarray(errors),direction_map_determinant=float(np.linalg.det(M)),source_controls_matched_by_proximity=False,recipient_section_shape_preserved=True,anatomical_acceptance=False)
