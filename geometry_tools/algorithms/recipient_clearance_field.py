"""Reusable fixed-surface queries with a union of independently closed tissues.

Combining closed meshes into one parity query can incorrectly mark their
overlap as outside. Each closed source keeps its own occupancy query. Bounds
only cull points that cannot be inside; original triangles determine queries.
Open surfaces have no invented solid volume. Final Float64 contact replay is
mandatory; these native distances are route proposals, not certificates.
"""
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex,library
from .hybrid_surface_index import HybridSurfaceIndex


class RecipientClearanceField:
    def __init__(self,recipients,*,margin_m=.00002):
        if not recipients or not 0<margin_m<=.001:raise ValueError('Explicit recipients and finite positive margin required')
        self.margin_m=float(margin_m);vertices=[];faces=[];offset=0;self.closed=[];self.labels=[];face_ends=[]
        frame=recipients[0][0].frame
        for mesh,closed,label in recipients:
            if mesh.frame!=frame:raise ValueError('Recipient query frames differ')
            vertices.append(mesh.vertices);faces.append(mesh.faces+offset);offset+=len(mesh.vertices)
            self.labels.append(label);face_ends.append((face_ends[-1]if face_ends else 0)+len(mesh.faces))
            if closed:
                self.closed.append((SurfaceIndex(mesh),mesh.vertices.min(0),mesh.vertices.max(0),label))
        self.surface=HybridSurfaceIndex(Mesh(np.concatenate(vertices),np.concatenate(faces),frame))
        self.face_ends=np.asarray(face_ends,int)
        self.query_points=0;self.occupancy_points=0

    def clearance(self,points):
        q=np.asarray(points,float)
        if q.ndim!=2 or q.shape[1]!=3 or not np.isfinite(q).all():raise ValueError('Finite query points required')
        self.query_points+=len(q);hit,_=self.surface.closest_points(q);gap=hit['distance_m']-self.margin_m
        for index,low,high,label in self.closed:
            chosen=np.flatnonzero(np.all(q>=low,axis=1)&np.all(q<=high,axis=1))
            if not len(chosen):continue
            self.occupancy_points+=len(chosen)
            tensor=library().core.Tensor(((q[chosen]-index.origin)/index.scale).astype(np.float32))
            inside=index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
            # Negative means blocked, not a claimed penetration-depth estimate.
            ids=chosen[inside];gap[ids]=np.minimum(gap[ids],-self.margin_m)
        return gap

    def segment_clear(self,a,b):
        return self.surface.segment_clear(a,b)

    def witnesses(self,points):
        q=np.asarray(points,float);hit,_=self.surface.closest_points(q);owner=np.searchsorted(self.face_ends,hit['face_ids'],side='right')
        rows=[dict(point_m=point.tolist(),nearestPart=self.labels[int(k)],distance_m=float(d),insideClosedSources=[])
              for point,k,d in zip(q,owner,hit['distance_m'])]
        for index,low,high,label in self.closed:
            chosen=np.flatnonzero(np.all(q>=low,axis=1)&np.all(q<=high,axis=1))
            if not len(chosen):continue
            tensor=library().core.Tensor(((q[chosen]-index.origin)/index.scale).astype(np.float32))
            inside=index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
            for i in chosen[inside]:rows[int(i)]['insideClosedSources'].append(label)
        return rows
