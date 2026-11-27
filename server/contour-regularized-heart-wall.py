"""Conforming contour of the saved regularized field; no QEM topology loss."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.regular_grid_isosurface import extract
import manifold3d as md
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);start=time.monotonic();source=b/'whole-wall-v7/field.npz';z=np.load(source);d=z['distance'];origin=z['origin'];h=float(z['spacing']);v,f,stats=extract(d,origin,h)
Mesh(v,f,'model:local-ta98:assembled-metres').save(out/'dense-wall.npz');print('Conforming contour saved',len(f),'faces',round(time.monotonic()-start,1),'s',flush=True)
solid=md.Manifold(md.Mesh64((v-origin).astype(float),f.astype('uint64')));status=str(solid.status());atomic_json(out/'contour.json',dict(status=status,extraction=stats,elapsed_seconds=time.monotonic()-start));assert status=='Error.NoError',status
pieces=solid.decompose();vol=np.array([p.volume()for p in pieces]);main=int(vol.argmax());assert all(abs(x)<h**3*8 for i,x in enumerate(vol)if i!=main),'Distinct source tissue component needs anatomical review'
solid=pieces[main];print('Native manifold valid; genus',solid.genus(),flush=True);simple=solid.simplify(.00015);packed=simple.to_mesh64();v=packed.vert_properties[:,:3]+origin;f=packed.tri_verts.astype('int64');mesh=Mesh(v,f,'model:local-ta98:assembled-metres');mesh.save(out/'whole-wall.npz')
result=dict(constructed=True,source_field_sha256=digest_file(source),closing_radius_m=1.5*h,simplification_tolerance_m=.00015,discarded_numerical_component_volumes_mL=np.delete(vol,main).tolist(),tissue_volume_mL=simple.volume()*1e6,genus=int(simple.genus()),components=len(simple.decompose()),vertices=len(v),triangles=len(f),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
