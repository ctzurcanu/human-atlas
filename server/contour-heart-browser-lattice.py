"""Contour on a well-conditioned lattice already exactly representable by GPU.

Two signed distance levels on each side bound contour interpolation error by
one quarter of the source voxel. Quarter-voxel edges cannot create critical
zero-corner slivers. The final simplifier retains a subset of the input
vertices, all on one binary Float32 lattice; no later rounding is hidden.
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
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'whole-wall-v7/field.npz';out=b/'whole-wall-v14';out.mkdir(exist_ok=True);assert not any(out.iterdir());z=np.load(source);d=z['distance'];h=float(z['spacing']);o=np.array([0.,1.375,.03125]);q=np.where(d<0,-1.,1.)*(.5+np.minimum(np.floor(abs(d)/h),1))*h;assert np.all(q!=0)
v,f,stats=extract(q,z['origin'],h);quantum=2.**-24;relative=np.rint((v-o)/quantum)*quantum;assert np.array_equal(relative,relative.astype('float32').astype(float));solid=md.Manifold(md.Mesh64(relative,f.astype('uint64')));assert str(solid.status())=='Error.NoError';print('Conditioned contour',len(f),'faces',round(time.monotonic()-start,1),'s',flush=True);pieces=solid.decompose();vol=np.array([p.volume()for p in pieces]);main=int(vol.argmax());extra=np.delete(vol,main);assert abs(extra).sum()*1e6<.12,'Unreviewed non-numerical tissue component'
simple=pieces[main].set_tolerance(.000001).simplify(.00012);packed=simple.to_mesh64();relative=packed.vert_properties[:,:3];faces=packed.tri_verts.astype('int64');assert np.array_equal(relative,relative.astype('float32').astype(float)),'Simplifier introduced non-source controls';mesh=Mesh(relative+o,faces,'model:local-ta98:assembled-metres');mesh.save(out/'whole-wall.npz');print('Browser wall',len(faces),'faces',round(time.monotonic()-start,1),'s',flush=True);proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof)
report=dict(passed=proof['passed'],source_field_sha256=digest_file(source),field_interpolation_error_bound_m=h/4,position_lattice_m=quantum,position_quantization_error_bound_m=3**.5*quantum/2,browser_origin=o.tolist(),post_simplification_rounding_error_m=0,closing_radius_m=1.5*h,simplification_tolerance_m=.00012,removed_small_component_volumes_mL=(extra*1e6).tolist(),volume_mL=simple.volume()*1e6,genus=int(simple.genus()),triangles=len(faces),vertices=len(relative),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
