"""Repair a pre-existing local native bile-wall crossing, preserving mouths."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json
from geometry_tools.audit import audit
b=ROOT/'.local-models/geometry-workbench/liver-completion';src=b/'female-inlets-v247-biliary';out=b/'female-inlets-v249-biliary';out.mkdir(exist_ok=False);start=time.monotonic();q=json.loads((src/'junction.json').read_text());check=json.loads((src/'parent-check.json').read_text());m=Mesh.load(src/'parent.npz');v,ix=np.unique(m.vertices,axis=0,return_inverse=True);f=ix[m.faces];tri=v[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal=np.zeros_like(v)
for j in range(3):np.add.at(normal,f[:,j],fn)
normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-30)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(path.read_text());row=next(p for p in a['parts']if p['id']=='LOCAL:female:bile_duct:1031');native=AtlasReader().read_part(a,row,frame='model:local-female-ta98:assembled-metres',catalogue_path=path);ids=np.r_[tuple(np.arange(g['start']//3,(g['start']+g['count'])//3)for g in row['tissueDrawGroups']if g.get('neighbor'))];protected=(native.triangles[ids].reshape(-1,3)-q['translation_m'])/q['scale'];origin=np.array([0,1.25,0]);protected=(protected-origin).astype('<f4').astype(float)+origin;keys={tuple(x)for x in protected};fixed=np.array([tuple(x)in keys for x in v]);fixed[np.unique(f[q['parent_interface_faces']])]=True
examples=np.asarray(check['self_intersections']['examples'],int);delta=np.zeros_like(v);protected_mouth=np.unique(f[q['parent_interface_faces']])
for outer in np.unique(examples[:,0]):
 inside=np.unique(f[examples[examples[:,0]==outer,1]]);inside=np.setdiff1d(inside,np.r_[f[outer],protected_mouth]);n=fn[outer]/np.linalg.norm(fn[outer]);signed=(v[inside]-v[f[outer,0]])@n;amount=np.maximum(0,signed+.00003);delta[inside]-=amount[:,None]*n
new=Mesh(((v+delta-origin).astype('<f4').astype(float)+origin),f,m.frame);new.save(out/'parent.npz');Mesh.load(src/'child.npz').save(out/'child.npz');result=audit(new,require_closed=True,max_pairs=1000000,max_tests=12000000);atomic_json(out/'parent-check.json',result);changed=np.flatnonzero(np.any(new.vertices!=v,axis=1));atomic_json(out/'shared-control-map.json',dict(source=v[changed].tolist(),destination=new.vertices[changed].tolist(),frame=m.frame,scale=q['scale'],translation_m=q['translation_m']));q.update(passed=result['passed'],native_wall_repair=dict(preexisting_crossing_pairs=check['self_intersections']['crossing_face_pairs'],maximum_surface_offset_m=float(np.linalg.norm(delta,axis=1).max()),shared_control_changes=int(np.count_nonzero(fixed[changed])),shared_control_map=str(out/'shared-control-map.json')),elapsed_seconds=time.monotonic()-start);atomic_json(out/'junction.json',q);print(result['passed'],result['failures'],flush=True)
if not result['passed']:raise ValueError('Native wall repair failed; not installed')
