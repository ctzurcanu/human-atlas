"""Orient one explicit complete source section to its ordered course bisector.

The caller supplies the two neighboring centres and the actual ordered ring.
No calibre, centre, source correspondence or anatomical branch is inferred.
"""
import numpy as np


def orient_section(vertices, ring, previous_center, next_center, *,
                   maximum_rotation_degrees, maximum_motion_m):
    v=np.asarray(vertices,float); ring=np.asarray(ring,int)
    if ring.ndim!=1 or len(ring)<3 or len(np.unique(ring))!=len(ring):
        raise ValueError('One complete ordered source ring required')
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all():
        raise ValueError('Finite corresponding source controls required')
    q=v[ring]; centre=q.mean(0)
    incoming=centre-np.asarray(previous_center,float)
    outgoing=np.asarray(next_center,float)-centre
    lengths=np.linalg.norm([incoming,outgoing],axis=1)
    if np.any(lengths<1e-9): raise ValueError('Distinct ordered source centres required')
    incoming/=lengths[0]; outgoing/=lengths[1]
    bisector=incoming+outgoing
    if np.linalg.norm(bisector)<1e-5:
        raise ValueError('A reversing course needs a new rounded path, not one rotated section')
    bisector/=np.linalg.norm(bisector)
    normal=np.cross(q-centre,np.roll(q,-1,axis=0)-centre).sum(0)
    if np.linalg.norm(normal)<1e-20: raise ValueError('Source profile area must be positive')
    normal/=np.linalg.norm(normal)
    # Preserve the directed incoming source winding. Choosing a sign only to
    # minimize an angle could reverse correspondence through a folded course.
    if abs(normal@incoming)<1e-4:
        raise ValueError('Incoming source winding is ambiguous')
    target=bisector*np.sign(normal@incoming)
    cosine=float(np.clip(normal@target,-1,1))
    angle=float(np.degrees(np.arccos(cosine)))
    if angle>maximum_rotation_degrees or cosine<=-1+1e-10:
        raise ValueError('Explicit source section rotation bound exceeded')
    axis=np.cross(normal,target)
    K=np.array([[0.,-axis[2],axis[1]],[axis[2],0.,-axis[0]],[-axis[1],axis[0],0.]])
    R=np.eye(3)+K+K@K/(1+cosine)
    assert np.linalg.det(R)>1-1e-12
    # Every exact original alias receives the same proper transformation.
    ids=np.flatnonzero(np.any(np.all(v[:,None]==q[None],axis=2),axis=1))
    new=v.copy();new[ids]=centre+(v[ids]-centre)@R.T
    movement=float(np.linalg.norm(new-v,axis=1).max())
    if movement>maximum_motion_m: raise ValueError('Explicit section movement bound exceeded')
    before=np.linalg.norm(q[:,None]-q[None],axis=2)
    after=np.linalg.norm(new[ring,None]-new[None,ring],axis=2)
    error=float(abs(after-before).max())
    if error>1e-12: raise ValueError('Proper section shape was not preserved')
    np.testing.assert_allclose(new[ring].mean(0),centre,atol=1e-14,rtol=0)
    return new,dict(complete=True,ring=ring.tolist(),rawIds=ids.tolist(),
                    previousCenter_m=np.asarray(previous_center).tolist(),nextCenter_m=np.asarray(next_center).tolist(),
                    center_m=centre.tolist(),rotation=R.tolist(),rotationDegrees=angle,
                    maximumMotion_m=movement,maximumSectionDistanceError_m=error,
                    originalSectionShapeAreaAndCentrePreserved=True,
                    noCalibreReduction=True,noAnatomicalAttachmentInferred=True)
