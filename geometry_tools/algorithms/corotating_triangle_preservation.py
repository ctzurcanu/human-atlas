"""Finite, source-bound preparation of rotating pair-plane derivatives."""
import time
import numpy as np
from scipy import sparse
from .triangle_cross_candidates import TriangleCrossIndex
from .section_pair_separation import separated_triangle_axes
from .corotating_pair_intervals import corotating_vertex_pair_rows

def rotating_pair_preservation(first,second,first_groups,second_groups,centers,levers,*,maximum_motion_m,deadline):
    ga,gb=np.asarray(first_groups),np.asarray(second_groups)
    moving=np.flatnonzero(np.any(ga[first.faces]>=0,axis=1))
    matrices=[];bounds=[];ledgers=[];gaps=[];tested=0
    for local in TriangleCrossIndex(second.triangles).batches(first.triangles[moving],padding_m=2*maximum_motion_m,batch_size=128,maximum_proposals=2000000):
        if time.monotonic()>=deadline:raise TimeoutError('Declared rotating-plane preparation deadline')
        pairs=np.c_[moving[local[:,0]],local[:,1]];tested+=len(pairs)
        axes,gap=separated_triangle_axes(first.triangles[pairs[:,0]],second.triangles[pairs[:,1]])
        keep=np.isfinite(gap)&(gap>1e-12);pairs=pairs[keep];axes=axes[keep];gap=gap[keep]
        if not len(pairs):continue
        fa,fb=first.faces[pairs[:,0]],second.faces[pairs[:,1]]
        va=np.broadcast_to(fa[:,:,None],(len(pairs),3,3)).ravel()
        vb=np.broadcast_to(fb[:,None,:],(len(pairs),3,3)).ravel()
        # Six source controls, including fixed stars, each own 1/6 of the
        # plane rotation. Fixed controls contribute zero angular velocity.
        owner=np.concatenate([ga[fa],gb[fb]],axis=1)
        plane_groups=np.repeat(owner,9,axis=0)
        weights=np.full(plane_groups.shape,1/6.)
        A,current=corotating_vertex_pair_rows(first.vertices[va],second.vertices[vb],ga[va],gb[vb],
            centers,levers,np.repeat(axes,9,axis=0),plane_groups,weights)
        rhs=np.repeat(np.minimum(2e-8,gap*.5),9)-current
        live=np.asarray(A.getnnz(1)).ravel()>0
        if np.any(rhs[~live]>1e-10):raise ValueError('Invalid fixed source pair')
        matrices.append(A[live]);bounds.append(rhs[live]);ledgers.append(pairs);gaps.append(gap)
    A=sparse.vstack(matrices,format='csr') if matrices else sparse.csr_matrix((0,6*len(centers)))
    rhs=np.concatenate(bounds) if bounds else np.empty(0)
    return A,rhs,dict(pairs=np.concatenate(ledgers) if ledgers else np.empty((0,2),int),
        original_separating_gap_m=np.concatenate(gaps) if gaps else np.empty(0)),dict(complete=True,
        testedPairs=tested,originalPositivePairs=sum(map(len,ledgers)),rows=len(rhs),
        finiteGeometryAccepted=False,actualEveryOriginalPairVerificationRequired=True,
        sourceGeometryModified=False,planeRotation='Mean angular velocity of six literal source controls; fixed controls stay fixed')
