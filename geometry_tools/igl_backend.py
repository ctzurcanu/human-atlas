"""Optional native Float64 nearest-surface queries on unchanged source triangles."""
import time
import numpy as np
from .distance import point_triangle


class SurfaceIndex:
    """Reusable native double-precision AABB over immutable copied coordinates."""

    def __init__(self, mesh):
        import igl
        if not len(mesh.faces):
            raise ValueError('A nonempty source surface is required')
        self.vertices=np.array(mesh.vertices,dtype=np.float64,copy=True)
        self.faces=np.array(mesh.faces,dtype=np.int64,copy=True)
        self.triangles=self.vertices[self.faces]
        self.geometry_sha=mesh.geometry_sha
        self.tree=igl.AABB();self.tree.init(self.vertices,self.faces)

    def closest_points(self, points):
        return _closest_points(self,points)


def closest_points(mesh, points):
    return SurfaceIndex(mesh).closest_points(points)


def _closest_points(index, points):
    """Return libigl nearest faces with original Float64 triangle witnesses.

    This is numerical double precision, not an exact-arithmetic certificate.
    Unlike Float32 ray-scene backends, it retains very thin source triangles.
    """
    started = time.monotonic()
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('Finite N by 3 surface queries required')
    squared, ids, native = index.tree.squared_distance(index.vertices,index.faces,points)
    ids = np.asarray(ids, dtype=np.int64)
    if np.any(ids<0)or np.any(ids>=len(index.faces)):
        raise ValueError('Native nearest-surface query returned an invalid face')
    d2, witnesses, barycentric = point_triangle(points, index.triangles[ids])
    if not np.isfinite(d2).all():
        raise ValueError('Native nearest-surface query returned an invalid witness')
    return dict(distance_m=np.sqrt(d2), points=witnesses, face_ids=ids,
                barycentric=barycentric), dict(
        backend='libigl', precision='Float64', source_geometry_sha256=index.geometry_sha,
        query_count=len(points), original_source_triangles=True,
        maximum_native_witness_difference_m=float(np.linalg.norm(witnesses-native, axis=1).max(initial=0)),
        elapsed_seconds=time.monotonic()-started, exact_arithmetic=False)
