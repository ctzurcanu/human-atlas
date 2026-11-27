"""Keep finite source attachment collars when opening ventricular outflows.

Guards are computational volumes only. They retain existing host tissue and
cannot create a new tissue solid, plug or duplicate valve.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,sys,time
import numpy as np
import manifold3d as md
from scipy.spatial import ConvexHull
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,classifier_caps_outside
ap=argparse.ArgumentParser();ap.add_argument('--bank',required=True);ap.add_argument('--fluids',required=True);ap.add_argument('--stage',required=True)
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
bank=h/args.bank;e=json.loads((bank/'emission.json').read_text());assert digest_file(bank/'cells.bin')==e['bank_sha256']
src=h/e['source_stage'];meta=json.loads((src/'ownership.json').read_text());z=np.load(src/'arrangement.npz');p=z['positions'].astype('<f4').astype(float)
surfaces={r['id']:z[r['faces']] for r in meta['surfaces']}
fields=json.loads((h/args.fluids/'fields.json').read_text());assert fields['passed'] and fields['source_bank_sha256']==e['bank_sha256'];hosts=set(fields['wall_part_ids'])
groups=defaultdict(list);contact_ids=set()
for pid in hosts:
 for use in meta['cells'][pid]:
  peer=use.get('neighbor')
  if peer in meta['cells'] and peer not in hosts:
   sid=use['surface'];contact_ids.add(sid)
   for triangle in p[surfaces[sid]]:
    key=(sid,*np.floor(triangle.mean(0)/.001).astype(int))
    groups[key].append(triangle)
directions=np.array([v for v in np.ndindex(3,3,3) if v!=(1,1,1)],float)-1
directions/=np.linalg.norm(directions,axis=1)[:,None]
guards=[]
for i,triangles in enumerate(groups.values()):
 if time.monotonic()-start>35:raise TimeoutError('Attachment guard construction deadline')
 vertices=np.unique(np.concatenate(triangles),axis=0)
 points=(vertices[:,None]+.0002*directions).reshape(-1,3)
 hull=ConvexHull(points);fs=hull.simplices.copy();t=points[fs]
 flip=np.einsum('ti,ti->t',np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),hull.equations[:,:3])<0
 fs[flip]=fs[flip][:,[0,2,1]]
 guards.append(make_solid(points,fs)[0])
print(len(guards),'local attachment guards; 0.2 mm construction envelope',flush=True)
guard=md.Manifold.batch_boolean(guards,md.OpType.Add)
assert guard.status()==md.Error.NoError
host=Mesh.load(h/args.fluids/'host-boundary.npz');host.save(out/'host-boundary.npz');records=[]
for i,old in enumerate(fields['records']):
 m=Mesh.load(h/args.fluids/f'fluid-{i}.npz');fluid,oid=make_solid(m.vertices,m.faces)
 final=fluid-guard;mm=final.to_mesh64();q=mm.vert_properties[:,:3];ff=np.asarray(mm.tri_verts)
 original,parent=provenance(mm);oldcaps=np.load(h/args.fluids/f'fluid-{i}-caps.npz')['cap_face_ids']
 caps=np.flatnonzero((original==oid)&np.isin(parent,oldcaps))
 proof=classifier_caps_outside(q,ff,caps,host.vertices,host.faces)
 mesh=Mesh(q,ff,host.frame);mesh.save(out/f'fluid-{i}.npz');atomic_npz(out/f'fluid-{i}-caps.npz',cap_face_ids=caps)
 records.append(dict(old,geometry_sha256=mesh.geometry_sha,cap_outside_tissue_proof=proof,
   attachment_guard_volume_removed_m3=float(fluid.volume()-final.volume()),surviving_computational_caps=len(caps)))
fields.update(records=records,source_fluid_stage=args.fluids,protected_source_charts=sorted(contact_ids),
 attachment_guard_envelope_m=.0002,guard_groups=len(guards),elapsed_seconds=time.monotonic()-start)
atomic_json(out/'fields.json',fields)
print(json.dumps(dict(passed=True,guard_groups=len(guards),elapsed_seconds=time.monotonic()-start)),flush=True)
