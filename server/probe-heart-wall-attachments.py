"""Current-frame, full-triangle contact evidence at declared cardiac attachments."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args()
start=time.monotonic();out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());pp={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';meshes={};rows=[]
def mesh(pid):
 if pid not in meshes:
  meshes[pid]=reader.read_part(a,pp[pid],frame=frame,catalogue_path=path)
 return meshes[pid]
for pid in a['finiteCardiacIntegration']['papillaryPartIds']:
 ventricle='left'if 'left ventricle'in pp[pid]['name'].lower()else'right'
 for host in [f'TA98HEART:male:myocardium:{ventricle}-ventricle','TA98HEART:male:endocardium']:
  if time.monotonic()-start>45:raise TimeoutError('Attachment probe stage deadline')
  first,second=mesh(pid),mesh(host)
  if np.any(np.maximum(first.bounds[0],second.bounds[0])>np.minimum(first.bounds[1],second.bounds[1])+.0003):
   row=dict(part_id=pid,host=host,bounds_disjoint_at_300um=True,attachment_inferred=False)
  else:
   witnesses,contact=pair_contacts(first,second,.0003,max_pairs=600000,max_tests=6000000)
   assert contact['complete_within_pair_budget'],(pid,host)
   atomic_npz(out/f'contact-{len(rows)}.npz',**witnesses)
   row=dict(part_id=pid,host=host,contact=contact,first_geometry_sha256=first.geometry_sha,second_geometry_sha256=second.geometry_sha,attachment_inferred=False)
  rows.append(row);atomic_json(out/'progress.json',dict(rows=rows,elapsed_seconds=time.monotonic()-start),overwrite=True)
  print(json.dumps(dict(part=pid,host=host,crossings=row.get('contact',{}).get('crossing_face_pairs'),coplanar=row.get('contact',{}).get('coplanar_overlap_face_pairs'),bounds_disjoint=row.get('bounds_disjoint_at_300um',False))),flush=True)
for i,(pid,m)in enumerate(meshes.items()):m.save(out/f'source-{i}.npz')
atomic_json(out/'contacts.json',dict(catalogue_sha256=digest_file(path),rows=rows,sources=[dict(part_id=pid,path=f'source-{i}.npz',geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha)for i,(pid,m)in enumerate(meshes.items())],elapsed_seconds=time.monotonic()-start,geometry_changed=False,anatomical_acceptance=False))
