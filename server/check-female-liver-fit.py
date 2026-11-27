"""Measure the actual transplanted exterior against fixed female diaphragm."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
base=ROOT/'.local-models/geometry-workbench/liver-completion/female-transfer-v206'
a=json.loads((base/'atlas.json').read_text());reader=AtlasReader();frame='model:local-female-ta98:assembled-metres'
host=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Liver:capsule');full=reader.read_part(a,host,frame=frame)
outer=Mesh(full.vertices,full.faces[:host['tissueDrawGroups'][0]['count']//3],frame);bvh=BVH(outer)
rows=[];fixed_ids=set(a['finiteDiaphragmIntegration']['partIds'])
fixed_ids.update(p['id']for p in a['parts']if p['id'].startswith('LOCAL:female:')and p['system']=='skeletal'and not p.get('suppressed'))
for p in a['parts']:
 if p.get('suppressed')or p['id']not in fixed_ids:continue
 mesh=reader.read_part(a,p,frame=frame)
 if np.any(mesh.bounds[1]<outer.bounds[0])or np.any(mesh.bounds[0]>outer.bounds[1]):continue
 _,r=pair_contacts(outer,mesh,0,first_bvh=bvh,max_pairs=100000,max_tests=2000000)
 rows.append(dict(part=p['id'],**r));print(p['id'],r['crossing_face_pairs'],r['coplanar_overlap_face_pairs'],r['truncated'],flush=True)
result=dict(fixedDiaphragm=True,fixedSkeleton=True,checks=rows,passed=bool(rows)and all(not r['crossing_face_pairs']and not r['coplanar_overlap_face_pairs']and not r['truncated']for r in rows),anatomicalAcceptance=False)
atomic_json(base/'fixed-context-fit.json',result);print(json.dumps(dict(passed=result['passed'],checks=len(rows))))
