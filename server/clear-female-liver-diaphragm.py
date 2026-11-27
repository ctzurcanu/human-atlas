"""Bounded caudal placement trials; never move the fixed female diaphragm."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_json
base=ROOT/'.local-models/geometry-workbench/liver-completion/female-transfer-v202'
a=json.loads((base/'atlas.json').read_text());reader=AtlasReader();frame='model:local-female-ta98:assembled-metres'
p=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Liver:capsule');full=reader.read_part(a,p,frame=frame)
outer=Mesh(full.vertices,full.faces[:p['tissueDrawGroups'][0]['count']//3],frame)
fixed_ids=set(a['finiteDiaphragmIntegration']['partIds']);targets=[]
for p in a['parts']:
 if p['id']in fixed_ids and not p.get('suppressed'):
  mesh=reader.read_part(a,p,frame=frame);targets.append((p['id'],mesh,BVH(mesh)))
trials=[];chosen=None
for dy in [-.004,-.008,-.012,-.016,-.020]:
 q=Mesh(outer.vertices+[0,dy,0],outer.faces,frame);b=BVH(q);rows=[]
 for pid,m,tb in targets:
  if np.any(m.bounds[1]<q.bounds[0])or np.any(m.bounds[0]>q.bounds[1]):continue
  _,r=pair_contacts(q,m,0,first_bvh=b,second_bvh=tb,max_pairs=100000,max_tests=2000000)
  rows.append(dict(part=pid,**r))
  if r['truncated']:raise RuntimeError('Finite contact budget exhausted')
 count=sum(r['crossing_face_pairs']+r['coplanar_overlap_face_pairs']for r in rows)
 trial=dict(translationDeltaMetres=[0,dy,0],crossings=count,checks=rows);trials.append(trial)
 print(json.dumps(dict(caudalMetres=-dy,crossings=count)),flush=True)
 if not count:chosen=trial;break
record=dict(trials=trials,chosen=chosen,passed=chosen is not None,fixedSkeleton=True,fixedDiaphragm=True,orientationAndScaleUnchanged=True)
atomic_json(base/'diaphragm-clearance-search.json',record)
if chosen is None:raise ValueError('No tested bounded caudal placement clears the fixed diaphragm')
