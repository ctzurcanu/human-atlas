"""Restore a numerically near-feasible conic displacement without relaxing bounds."""
import numpy as np
def retract_origin_feasible(displacement,hard,rhs,*,maximum_motion_m,distortion=None,distortion_bound=.05,maximum_relative_retraction=1e-6):
 d=np.asarray(displacement,float);b=np.asarray(rhs,float)
 if d.ndim!=2 or d.shape[1]!=3 or not np.isfinite(d).all() or not np.isfinite(b).all() or maximum_motion_m<=0 or not 0<maximum_relative_retraction<=1e-4:raise ValueError('Finite bounded near-feasible displacement required')
 motion=float(np.linalg.norm(d,axis=1).max(initial=0));tangent=float(np.linalg.norm((distortion@d.ravel()).reshape(-1,6),axis=1).max(initial=0)) if distortion is not None else 0.
 factor=1.
 if motion>maximum_motion_m:factor=min(factor,maximum_motion_m/(motion+1e-14))
 if distortion is not None and tangent>distortion_bound:factor=min(factor,distortion_bound/(tangent+1e-11))
 if factor<1:
  if np.any(b>0):raise ValueError('The original zero displacement is not feasible for all hard inequalities')
  if 1-factor>maximum_relative_retraction:raise ValueError('Numerical violation is larger than the declared tiny retraction budget')
 result=d*factor;residual=float(np.maximum(b-hard@result.ravel(),0).max(initial=0));after_motion=float(np.linalg.norm(result,axis=1).max(initial=0));after_tangent=float(np.linalg.norm((distortion@result.ravel()).reshape(-1,6),axis=1).max(initial=0)) if distortion is not None else 0.
 if residual>1e-10 or after_motion>maximum_motion_m+1e-12 or after_tangent>distortion_bound+1e-12:raise ValueError('Retracted point still fails the original independent bounds')
 return result,dict(applied=factor<1,factor=factor,maximumRelativeRetraction=maximum_relative_retraction,allHardRhsNonpositive=bool(np.all(b<=0)),maximumOriginalHardResidual_m=residual,maximumMotionBefore_m=motion,maximumMotionAfter_m=after_motion,maximumTangentBefore=tangent,maximumTangentAfter=after_tangent,originalBoundsRelaxed=False,maximumCorrection_m=float(np.linalg.norm(result-d,axis=1).max(initial=0)))
