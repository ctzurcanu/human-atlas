"""Convex tangent bounds about rotated, unscaled original facet frames.

If ||F-U||_F <= e and U has orthonormal columns, every singular value of F
lies in [1-e,1+e]. U is computed from a reviewed pose, but the whitening and
translation variables always refer to the original source. Changing the bend
frame therefore cannot reset or accumulate a physical stretch allowance.
"""
import numpy as np
from .source_facet_distortion import facet_tangent_rows


def corotated_tangent_rows(reference, reviewed, stations, station_count):
    np.testing.assert_array_equal(reference.faces, reviewed.faces)
    assert reference.frame == reviewed.frame
    t,q=reference.triangles,reviewed.triangles
    E=(t[:,1:]-t[:,:1]).transpose(0,2,1)
    Q=(q[:,1:]-q[:,:1]).transpose(0,2,1)
    _,s,v=np.linalg.svd(E,full_matrices=False)
    assert np.isfinite(s).all() and np.all(s[:,1]>0)
    R=np.swapaxes(v,1,2)/s[:,None,:]
    F=Q@R
    u,values,w=np.linalg.svd(F,full_matrices=False)
    assert np.all(values[:,1]>0), 'Reviewed degenerate facet cannot define a bend frame'
    target=u@w
    np.testing.assert_allclose(target.transpose(0,2,1)@target,np.broadcast_to(np.eye(2),(len(t),2,2)),atol=2e-14)
    original=E@R
    original3=np.concatenate([original,np.cross(original[:,:,0],original[:,:,1])[:,:,None]],axis=2)
    target3=np.concatenate([target,np.cross(target[:,:,0],target[:,:,1])[:,:,None]],axis=2)
    rotations=target3@original3.transpose(0,2,1)
    assert np.all(np.linalg.det(rotations)>1-2e-13), 'Reflection is forbidden'
    offset=(original-target).transpose(0,2,1).reshape(-1,6)
    C,proof=facet_tangent_rows(reference,stations,station_count)
    return C,offset,dict(**proof,reviewedFrameOnlyNotReferenceReset=True,
        properUnscaledBendFrames=True,maximumReviewedOriginalSingularStrain=float(np.abs(values-1).max()),
        maximumReviewedCorotatedResidual=float(np.linalg.norm(F-target,axis=(1,2)).max()),
        physicalGuarantee='Frobenius distance to a proper rotated original unit tangent bounds original singular strain')
