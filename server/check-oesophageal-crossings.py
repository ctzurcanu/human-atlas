"""Measure the oesophageal course against fixed adjoining organ tissues."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);args=ap.parse_args();sex=args.sex;b=ROOT/'.local-models/geometry-workbench/liver-completion/final-ready-v265';a=json.loads((b/(sex+'.json')).read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres';start=time.monotonic();out=b/(sex+'-oesophagus-crossings');out.mkdir(exist_ok=False);rows=[]
targets=[parts['TA98PART:ZA:Liver:capsule']]+[parts[x]for x in a['finiteDiaphragmIntegration']['partIds']if ':Muscle:'in x or x.endswith(':central')]
for pid in a['oesophagealWallGeometry']['partIds']:
 if not pid.endswith(':external-connective-coat'):continue
 row=parts[pid];whole=r.read_part(a,row,frame=frame);count=row['tissueDrawGroups'][0]['count']//3;outer=Mesh(whole.vertices,whole.faces[:count],frame);bv=BVH(outer)
 for p in targets:
  bb=np.array(p['bounds'])
  if np.any(bb[1]<outer.bounds[0])or np.any(bb[0]>outer.bounds[1]):continue
  peer=r.read_part(a,p,frame=frame)
  if p['id']=='TA98PART:ZA:Liver:capsule':peer=Mesh(peer.vertices,peer.faces[:p['liverReferenceTriangleCount']],frame)
  pairs,q=pair_contacts(outer,peer,0,first_bvh=bv,max_pairs=200000,max_tests=5000000);name=pid.replace(':','_')+'--'+p['id'].replace(':','_');atomic_json(out/(name+'.json'),q);rows.append(dict(source=pid,target=p['id'],crossings=q['crossing_face_pairs'],coplanar=q['coplanar_overlap_face_pairs'],truncated=q['truncated'],examples=q.get('examples',[])));print(pid,p['id'],q['crossing_face_pairs'],flush=True)
atomic_json(out/'checks.json',dict(passed=all(not x['crossings']and not x['truncated']for x in rows),checks=rows,diaphragm_fixed=True,elapsed_seconds=time.monotonic()-start))
