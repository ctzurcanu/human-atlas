"""Avoid critical zero-valued grid corners with a measured 10 µm isolevel.

This source-bound contour offset is below the .55 mm acquisition spacing;
it is recorded, not described as an anatomical myocardial thickness.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.regular_grid_isosurface import extract
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';p=b/'whole-wall-v7/field.npz';out=b/'whole-wall-v11';out.mkdir(exist_ok=False);z=np.load(p);o=z['origin'];h=float(z['spacing']);v,f,stats=extract(z['distance']+.00001,o,h);print('Stable contour',len(f),'faces',round(time.monotonic()-start,1),'s',flush=True);solid=md.Manifold(md.Mesh64(v-o,f.astype('uint64')));assert str(solid.status())=='Error.NoError';pieces=solid.decompose();vol=np.array([q.volume()for q in pieces]);main=int(vol.argmax());assert all(x<0 and abs(x)*1e6<.1 for i,x in enumerate(vol)if i!=main);assert abs(np.delete(vol,main)).sum()*1e6<.1
simple=pieces[main].set_tolerance(.000001).simplify(.00015);packed=simple.to_mesh64();relative=packed.vert_properties[:,:3];faces=packed.tri_verts.astype('int64');mesh=Mesh(relative.astype('float32').astype(float)+o,faces,'model:local-ta98:assembled-metres');mesh.save(out/'whole-wall.npz');print('Encoded wall',len(faces),'faces',round(time.monotonic()-start,1),'s',flush=True);proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof)
result=dict(passed=proof['passed'],constructed=True,source_field_sha256=digest_file(p),isolevel_offset_m=.00001,closing_radius_m=1.5*h,numerical_tolerance_m=.000001,simplification_tolerance_m=.00015,filled_numerical_voids_mL=(-np.delete(vol,main)*1e6).tolist(),volume_mL=simple.volume()*1e6,genus=int(simple.genus()),triangles=len(faces),vertices=len(relative),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
