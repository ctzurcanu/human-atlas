"""True surface-bound covering prototype, never a nodal distance extrapolation."""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';source=base/'whole-wall-v25/whole-wall.npz';out=base/'surface-covering-v43';out.mkdir(exist_ok=False);m=Mesh.load(source);o=np.array([0.,1.375,.03125]);a=md.Manifold(md.Mesh64(m.vertices-o,m.faces.astype('uint64')));assert str(a.status())=='Error.NoError';thickness=.00004;kernel=md.Manifold.hull_points(np.r_[np.eye(3),-np.eye(3)]*thickness);print('Erode fixed wall with eight-face convex kernel',flush=True);core=a.minkowski_difference(kernel);assert str(core.status())=='Error.NoError';volume=core.volume();print('Erosion complete',round(time.monotonic()-start,2),'s',core.num_tri(),'triangles',flush=True);coating=a-core
for name,solid in [('core',core),('covering',coating)]:
 z=solid.to_mesh64();Mesh(z.vert_properties[:,:3]+o,z.tri_verts.astype('int64'),m.frame).save(out/(name+'.npz'))
report=dict(constructed=True,source_sha256=digest_file(source),kernel='Convex octahedron',kernel_axis_radius_m=thickness,source_volume_mL=a.volume()*1e6,core_volume_mL=volume*1e6,covering_volume_mL=coating.volume()*1e6,core_triangles=core.num_tri(),covering_triangles=coating.num_tri(),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
