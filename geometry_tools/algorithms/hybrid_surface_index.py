"""Native regular-face queries plus Float64 fallback for every collapsed face.

This partitions query acceleration only. The literal source face bank is
unchanged; tiny or degenerate source faces are never silently discarded.
"""
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import library
from geometry_tools.spatial import BVH,box_distance_squared
from geometry_tools.distance import point_triangle,triangle_pairs


class HybridSurfaceIndex:
    def __init__(self,mesh):
        self.mesh=mesh;used,inv=np.unique(mesh.faces,return_inverse=True);v=mesh.vertices[used]
        self.origin=(v.min(0)+v.max(0))*.5;self.scale=float(np.ptp(v,axis=0).max())
        if self.scale<=0:self.scale=1.
        local=((v-self.origin)/self.scale).astype(np.float32);f=inv.reshape(-1,3).astype(np.uint32)
        t=local[f].astype(float);bad=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)==0
        self.native_ids=np.flatnonzero(~bad);self.fallback_ids=np.flatnonzero(bad);self.scene=None;self.fallback=None
        if len(self.native_ids):
            self.scene=library().t.geometry.RaycastingScene(nthreads=2)
            self.scene.add_triangles(library().core.Tensor(local),library().core.Tensor(f[~bad]))
        if len(self.fallback_ids):
            self.fallback=BVH(Mesh(mesh.vertices,mesh.faces[bad],mesh.frame))
            self.fallback_low=self.fallback.triangles.min(1);self.fallback_high=self.fallback.triangles.max(1)

    def closest_points(self,points):
        q=np.asarray(points,float);n=len(q)
        if q.ndim!=2 or q.shape[1]!=3 or not np.isfinite(q).all():raise ValueError('Finite metre queries required')
        if self.scene is not None and n:
            raw=self.scene.compute_closest_points(library().core.Tensor(((q-self.origin)/self.scale).astype(np.float32)),nthreads=2)
            ids=self.native_ids[raw['primitive_ids'].numpy().astype(int)]
            d,p,w=point_triangle(q,self.mesh.triangles[ids])
        else:
            ids=np.full(n,-1,int);d=np.full(n,np.inf);p=np.zeros((n,3));w=np.zeros((n,3))
        tested=0
        if self.fallback is not None:
            # A literal face-box lower bound limits the slow exact path. It
            # cannot discard a fallback face closer than the native witness.
            for first in range(0,n,2048):
                batch=q[first:first+2048]
                lower=box_distance_squared(batch[:,None],batch[:,None],self.fallback_low[None],self.fallback_high[None]).min(1)
                for k in np.flatnonzero(lower<=d[first:first+len(batch)]):
                    i=first+int(k);hit=self.fallback.closest_point(q[i]);tested+=hit['tested_triangles']
                    if hit['distance_m']**2<d[i]:
                        ids[i]=self.fallback_ids[hit['face_id']];d[i]=hit['distance_m']**2
                        p[i]=hit['point'];w[i]=hit['barycentric']
        return dict(distance_m=np.sqrt(d),points=p,face_ids=ids,barycentric=w),dict(
            nativeFaces=len(self.native_ids),float64FallbackFaces=len(self.fallback_ids),
            float64FallbackTests=tested,allSourceFacesRetained=True,nearestTriangleVerified=False)

    def segment_clear(self,a,b):
        a=np.asarray(a,float);b=np.asarray(b,float);clear=np.ones(len(a),bool)
        if self.scene is not None:
            rays=np.c_[(a-self.origin)/self.scale,(b-a)/self.scale]
            hit=self.scene.cast_rays(library().core.Tensor(rays.astype(np.float32)),nthreads=2)['t_hit'].numpy()
            clear&=~((hit>=0)&(hit<=1))
        if self.fallback is not None:
            for first in range(0,len(a),1024):
                aa=a[first:first+1024];bb=b[first:first+1024]
                low=np.minimum(aa,bb);high=np.maximum(aa,bb)
                pairs=np.argwhere(np.all(high[:,None]>=self.fallback_low[None],axis=2)&np.all(low[:,None]<=self.fallback_high[None],axis=2))
                for start in range(0,len(pairs),4096):
                    chosen=pairs[start:start+4096];segment=np.stack([aa[chosen[:,0]],bb[chosen[:,0]],bb[chosen[:,0]]],axis=1)
                    hit=triangle_pairs(segment,self.fallback.triangles[chosen[:,1]])
                    clear[first+chosen[hit['distance_m']<=1e-12,0]]=False
        return clear
