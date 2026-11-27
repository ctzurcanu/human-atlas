"""Render-only native scene, preserving every source triangle and its index.

Float32 projection is display evidence, not a geometric clearance predicate.
Even subpixel/degenerate projected faces remain in the source-index sequence.
Numerical SurfaceIndex queries keep their stricter collapse rejection.
"""
import numpy as np
from geometry_tools.open3d_backend import library

class SourceTriangleRender:
 def __init__(self,mesh,*,nthreads=2):
  if not len(mesh.faces):raise ValueError('Rendering requires source triangles')
  used,inv=np.unique(mesh.faces,return_inverse=True);v=mesh.vertices[used];self.mesh=mesh;self.origin=(v.min(0)+v.max(0))*.5;self.scale=float(np.ptp(v,axis=0).max())
  if not self.scale>0:raise ValueError('Point-only source cannot be rendered')
  local=((v-self.origin)/self.scale).astype(np.float32);faces=inv.reshape(-1,3).astype(np.uint32);native=local[faces].astype(float);original=v[inv.reshape(-1,3)];area=lambda t:np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1);zero=area(original)==0;collapsed=area(native)==0;o3d=library();self.scene=o3d.t.geometry.RaycastingScene(nthreads=nthreads);self.scene.add_triangles(o3d.core.Tensor(local),o3d.core.Tensor(faces));self.render_evidence=dict(sourceTriangleCount=len(mesh.faces),sourceTriangleOrderPreserved=True,originalZeroAreaFaces=int(zero.sum()),float32NewlyCollapsedFaces=int((collapsed&~zero).sum()),maximumVertexRounding_m=float(np.linalg.norm(local.astype(float)*self.scale+self.origin-v,axis=1).max()),geometryModified=False,renderOnly=True,clearancePredicate=False)
