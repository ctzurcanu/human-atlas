"""Replace only local paired-vessel plane surrogates, retaining a face ledger.

Omitting a surrogate proves no clearance. The caller must recheck every
originally separated pair on exact emitted surfaces, including omitted pairs.
All nonlocal and partly fixed triangle pairs retain the original halfspaces.
"""
import time
import numpy as np
from scipy import sparse
from .triangle_cross_candidates import TriangleCrossIndex
from .section_pair_separation import separated_triangle_axes,triangle_axis_sparse_rows


def regional_pair_preservation(first,second,first_station,second_station,station_count,*,
        first_region,second_region,maximum_motion_m,deadline=None):
    a,b=np.asarray(first_station),np.asarray(second_station)
    ra,rb=np.asarray(first_region,bool),np.asarray(second_region,bool)
    if ra.shape!=(len(first.vertices),) or rb.shape!=(len(second.vertices),):
        raise ValueError('Explicit original source-control correction-region masks required')
    moving=np.flatnonzero(np.any(a[first.faces]>=0,axis=1))
    blocks=[];bounds=[];omitted=[];old_gap=[];kept=0;prior=0;tested=0
    padding=maximum_motion_m*(2 if np.any(b>=0) else 1)
    for local in TriangleCrossIndex(second.triangles).batches(first.triangles[moving],padding_m=padding,
            batch_size=128,maximum_proposals=2000000):
        if deadline is not None and time.monotonic()>=deadline:
            raise TimeoutError('Declared regional pair-preservation deadline reached')
        pair=np.c_[moving[local[:,0]],local[:,1]]
        if not len(pair):continue
        tested+=len(pair)
        axis,gap=separated_triangle_axes(first.triangles[pair[:,0]],second.triangles[pair[:,1]])
        positive=np.isfinite(gap)&(gap>1e-12);prior+=int((~positive).sum())
        fa,fb=first.faces[pair[:,0]],second.faces[pair[:,1]]
        local_pair=(np.all(ra[fa],axis=1)&np.all(rb[fb],axis=1)
                    &np.all(a[fa]>=0,axis=1)&np.all(b[fb]>=0,axis=1))
        replace=positive&local_pair
        omitted.append(pair[replace]);old_gap.append(gap[replace])
        retain=positive&~local_pair;kept+=int(retain.sum())
        if not retain.any():continue
        A,rhs,_=triangle_axis_sparse_rows(first,second,a,b,pair[retain],axis[retain],
            station_count=station_count,gap_m=np.minimum(2e-8,gap[retain]*.5))
        live=np.asarray(A.getnnz(1)).ravel()>0
        if np.any(rhs[~live]>1e-10):raise ValueError('Unresolved fixed positive pair interval')
        blocks.append(A[live]);bounds.append(rhs[live])
    A=sparse.vstack(blocks,format='csr') if blocks else sparse.csr_matrix((0,3*station_count))
    rhs=np.concatenate(bounds) if bounds else np.array([])
    pairs=np.concatenate(omitted) if omitted else np.empty((0,2),int)
    gap=np.concatenate(old_gap) if old_gap else np.array([])
    return A,rhs,dict(pairs=pairs,original_separating_gap_m=gap),dict(complete=True,
        testedPairs=tested,retainedOriginalPositivePairs=kept,replacedLocalPositivePairSurrogates=len(pairs),
        preexistingOverlappingOrTouchingPairs=prior,rows=len(rhs),
        nonlocalAndPartlyFixedPairIntervalsRetained=True,
        actualEveryOriginalPairVerificationRequired=True,clearanceEstablished=False,
        sourceGeometryModified=False)
