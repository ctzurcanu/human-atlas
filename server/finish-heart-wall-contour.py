"""Review-stage removal of closed numerical voids, with explicit bounds."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
import manifold3d as md
base=ROOT/'.local-models/geometry-workbench/heart';source=base/'whole-wall-v8/dense-wall.npz';out=base/'whole-wall-v9';out.mkdir(exist_ok=False);start=time.monotonic();m=Mesh.load(source);origin=np.array([0.,1.38,.03]);solid=md.Manifold(md.Mesh64(m.vertices-origin,m.faces.astype('uint64')));assert str(solid.status())=='Error.NoError'
pieces=solid.decompose();vol=np.array([q.volume()for q in pieces]);main=int(vol.argmax());voids=[]
for i,q in enumerate(pieces):
 if i==main:continue
 b=np.array(q.bounding_box());assert vol[i]<0 and abs(vol[i])*1e6<.1 and np.linalg.norm(b[3:]-b[:3])<.015,'Non-numerical component needs separate review'
 voids.append(dict(volume_mL=vol[i]*1e6,bounds_m=(b+np.r_[origin,origin]).tolist()))
assert sum(abs(v['volume_mL'])for v in voids)<.1
simple=pieces[main].simplify(.00015);encoded=simple.to_mesh64();v=encoded.vert_properties[:,:3]+origin;f=encoded.tri_verts.astype('int64');mesh=Mesh(v,f,m.frame);mesh.save(out/'whole-wall.npz')
result=dict(constructed=True,source_sha256=digest_file(source),explicit_void_limit_mL=.1,filled_closed_numerical_voids=voids,total_filled_mL=sum(abs(v['volume_mL'])for v in voids),simplification_tolerance_m=.00015,volume_mL=simple.volume()*1e6,genus=int(simple.genus()),components=len(simple.decompose()),vertices=len(v),triangles=len(f),elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False)
atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
