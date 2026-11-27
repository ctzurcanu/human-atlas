"""Preserve actual triangle intervals for an explicitly supported motion field."""
import time
import numpy as np
from scipy import sparse
from .triangle_cross_candidates import TriangleCrossIndex
from .section_pair_separation import separated_triangle_axes,triangle_axis_sparse_rows


def preservation_rows(first,second,first_station,second_station,station_count,*,maximum_motion_m,self_pair=False,deadline=None):
    a,b=np.asarray(first_station),np.asarray(second_station)
    moving=np.flatnonzero(np.any(a[first.faces]>=0,axis=1))
    blocks=[];bounds=[];tested=0;retained=0;prior_overlap=0
    padding=maximum_motion_m*(2 if np.any(b>=0)else 1)
    for local in TriangleCrossIndex(second.triangles).batches(first.triangles[moving],padding_m=padding,batch_size=128,maximum_proposals=2000000):
        if deadline is not None:assert time.monotonic()<deadline,'Finite preservation-row stage'
        pair=np.c_[moving[local[:,0]],local[:,1]]
        if self_pair:
            # Include both changed-versus-fixed orders; deduplicate only when
            # the reciprocal first face is also in the moving subset.
            reciprocal=np.isin(pair[:,1],moving)
            pair=pair[~reciprocal|(pair[:,0]<pair[:,1])]
            fa,fb=first.faces[pair[:,0]],second.faces[pair[:,1]]
            shared=np.any(np.all(first.vertices[fa][:,:,None]==second.vertices[fb][:,None],axis=3),axis=(1,2))
            pair=pair[~shared]
        if not len(pair):continue
        tested+=len(pair)
        axis,gap=separated_triangle_axes(first.triangles[pair[:,0]],second.triangles[pair[:,1]])
        keep=np.isfinite(gap)&(gap>1e-12)
        prior_overlap+=int((~keep).sum());retained+=int(keep.sum())
        if not keep.any():continue
        A,rhs,_=triangle_axis_sparse_rows(first,second,a,b,pair[keep],axis[keep],station_count=station_count,gap_m=np.minimum(2e-8,gap[keep]*.5))
        live=np.asarray(A.getnnz(1)).ravel()>0
        assert np.all(rhs[~live]<=1e-10)
        blocks.append(A[live]);bounds.append(rhs[live])
    A=sparse.vstack(blocks,format='csr')if blocks else sparse.csr_matrix((0,3*station_count))
    rhs=np.concatenate(bounds)if bounds else np.empty(0)
    return A,rhs,dict(complete=True,sourceMovingFaceIds=moving.tolist(),testedPairs=tested,preservedIntervals=retained,
                      preexistingOverlappingOrNumericallyTouchingPairs=prior_overlap,rows=len(rhs),
                      strictPositiveGapsPreserved=True,wholeAnatomicalDisjointnessEstablished=False)
