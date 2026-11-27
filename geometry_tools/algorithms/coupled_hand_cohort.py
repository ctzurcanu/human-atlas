"""Five corresponding digit poses blended as one complete hand control field.

Finite source profiles and shared mouths remain rigid and coincident. This
initializer does not certify tissue clearance, patency or anatomical ownership.
"""
import numpy as np
from scipy.spatial.transform import Rotation
from .local_bone_cohort import local_bone_cohort,polyline_distance


def coupled_hand_cohort(vertices,sections,source_bones,target_bones,source_normal,target_normal,
                        baseline_vertices,*,radial_fields=None):
    if set(source_bones)!=set(range(1,6))or set(target_bones)!=set(range(1,6)):
        raise ValueError('All five corresponding named digits required')
    if radial_fields is not None and set(radial_fields)!=set(range(1,6)):
        raise ValueError('Every digit requires its measured radial field')
    v=np.asarray(vertices,float);base=np.asarray(baseline_vertices,float)
    if base.shape!=v.shape:raise ValueError('One corresponding installed baseline control bank required')
    def weights(points):
        d=np.column_stack([polyline_distance(points,source_bones[k])for k in range(1,6)])
        # The cubic compact support vanishes with zero first derivative.
        # A digit's radial map is never extrapolated into unrelated digits.
        support=np.maximum(0.,1-(d-d.min(1)[:,None])/.012)**3
        w=support*(d+.002)**-4
        return w/w.sum(1)[:,None]
    poses=[];rotations=[];evidence=[]
    for digit in range(1,6):
        q,R,proof=local_bone_cohort(v,sections,source_bones[digit],target_bones[digit],source_normal,target_normal,
            [source_bones[k]for k in range(1,6)if k!=digit],baseline_vertices=base,
            radial_scale=None if radial_fields is None else radial_fields[digit],neighbour_partition=False,
            pose_support=lambda p:weights(p)[:,digit-1]>0)
        poses.append(q);rotations.append(R);evidence.append(dict(digit=digit,**proof))
    poses=np.stack(poses);rotations=np.stack(rotations)
    q=np.einsum('ni,inj->nj',weights(v),poses)
    unchanged=np.all(poses==base[None],axis=(0,2));q[unchanged]=base[unchanged]
    for section in sections:
        if np.all(unchanged[section]):continue
        c=v[section].mean(0);w=weights(c[None])[0]
        centres=np.array([p[section].mean(0)for p in poses]);target=w@centres
        quat=Rotation.from_matrix(rotations[:,section[0]]).as_quat()
        # A quaternion outer-product mean is invariant to q versus -q.
        M=np.einsum('n,ni,nj->ij',w,quat,quat);_,vectors=np.linalg.eigh(M)
        R=Rotation.from_quat(vectors[:,-1]).as_matrix()
        q[section]=(v[section]-c)@R.T+target
    return q,dict(method='One normalized five-digit corresponding bone/skin field; source profiles use proper quaternion rotation means',
        namedDigits=[1,2,3,4,5],cohortMembers=evidence,sourceFacesRemoved=0,skeletonChanged=False,
        installedBaselineOutsideSupportPreserved=True,anatomicalAcceptance=False,
        actualSourceAndRecipientTriangleAuditRequired=True)
