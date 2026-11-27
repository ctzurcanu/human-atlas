"""Conservative cross-mesh triangle-box proposals using bounded sphere buckets.

The native KD tree is only a broad phase. Each triangle is enclosed by a
centroid sphere, bucket radii are rounded up, and actual Float64 boxes filter
the proposals. The caller must apply the existing triangle predicates.
"""
import numpy as np
from scipy.spatial import cKDTree


class TriangleCrossIndex:
    def __init__(self, triangles):
        self.triangles=np.asarray(triangles,float)
        if self.triangles.ndim!=3 or self.triangles.shape[1:]!=(3,3) or not np.isfinite(self.triangles).all():
            raise ValueError('Finite literal triangles required')
        t=self.triangles;self.low=t.min(1);self.high=t.max(1);centers=t.mean(1)
        radii=np.linalg.norm(t-centers[:,None],axis=2).max(1)
        # Buckets prevent one large source triangle inflating every query.
        powers=np.floor(np.log2(np.maximum(radii,np.finfo(float).tiny)))
        self.buckets=[]
        for power in np.unique(powers):
            ids=np.flatnonzero(powers==power)
            self.buckets.append((ids,cKDTree(centers[ids]),float(radii[ids].max())))

    def batches(self, triangles, *, batch_size=1024, maximum_proposals=4000000, padding_m=0.):
        t=np.asarray(triangles,float)
        if t.ndim!=3 or t.shape[1:]!=(3,3) or not np.isfinite(t).all():
            raise ValueError('Finite literal triangles required')
        if not np.isfinite(padding_m)or padding_m<0:raise ValueError('Finite nonnegative triangle-box padding required')
        centers=t.mean(1);radii=np.linalg.norm(t-centers[:,None],axis=2).max(1)
        low,high=t.min(1),t.max(1);total=0
        for start in range(0,len(t),batch_size):
            stop=min(start+batch_size,len(t))
            for ids,tree,radius in self.buckets:
                # A small scale-aware padding covers floating arithmetic in
                # centroid/radius evaluation; it cannot remove true proposals.
                pad=64*np.finfo(float).eps*(1+np.abs(centers[start:stop]).max(1)+radius)
                neighbours=tree.query_ball_point(centers[start:stop],radii[start:stop]+radius+pad+padding_m,workers=2)
                lengths=np.fromiter(map(len,neighbours),dtype=int,count=stop-start)
                if not lengths.sum():continue
                first=np.repeat(np.arange(start,stop),lengths)
                second=ids[np.concatenate(neighbours).astype(int)]
                keep=np.all(high[first]+padding_m>=self.low[second],axis=1)&np.all(self.high[second]+padding_m>=low[first],axis=1)
                pairs=np.c_[first[keep],second[keep]];total+=len(pairs)
                if total>maximum_proposals:raise RuntimeError('Finite cross-triangle proposal budget exceeded')
                if len(pairs):yield pairs
