"""Bounded end-fixed course changes with the connected gastric wall."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_json
from geometry_tools.oesophageal_clearance_field import transform
ap=argparse.ArgumentParser();ap.add_argument('--sex',required=True,choices=['male','female']);args=ap.parse_args();sex=args.sex;b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/('oesophageal-course-v274-'+sex);out.mkdir(exist_ok=False);a=json.loads((b/'final-ready-v265'/(sex+'.json')).read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres';start=time.monotonic();moving=[];fixed=[]
for pid in a['oesophagealWallGeometry']['partIds']+[pid for pid in a['gastricCompleteWall']['partIds']if ':019:'in pid]:
 if not(pid.startswith('GASTRIC:')or pid.endswith(':external-connective-coat')):continue
 p=parts[pid];m=r.read_part(a,p,frame=frame);moving.append((pid,Mesh(m.vertices,m.faces[:p['tissueDrawGroups'][0]['count']//3],frame)))
for pid in ['TA98PART:ZA:Liver:capsule']+a['finiteDiaphragmIntegration']['partIds']:
 p=parts[pid];m=r.read_part(a,p,frame=frame)
 if pid.endswith(':Liver:capsule'):m=Mesh(m.vertices,m.faces[:p['liverReferenceTriangleCount']],frame)
 fixed.append((pid,m,BVH(m)))
trials=[];chosen=None
for displacement in [[0,0,-.003],[.002,0,-.006],[.004,0,-.009],[.006,0,-.012]]:
 rows=[]
 for pid,m in moving:
  vertices,J=transform(m.vertices,displacement,sex=sex);q=Mesh(vertices,m.faces,frame);bv=BVH(q)
  for fid,peer,pb in fixed:
   if np.any(peer.bounds[1]<q.bounds[0])or np.any(peer.bounds[0]>q.bounds[1]):continue
   _,proof=pair_contacts(q,peer,0,first_bvh=bv,second_bvh=pb,max_pairs=200000,max_tests=5000000);rows.append(dict(moving=pid,fixed=fid,crossings=proof['crossing_face_pairs'],coplanar=proof['coplanar_overlap_face_pairs'],truncated=proof['truncated']))
 row=dict(displacement_m=displacement,checks=rows,constant_y_section_area_preserved=True,determinant=1.);trials.append(row);print(sex,displacement,[(q['moving'],q['fixed'],q['crossings'])for q in rows if q['crossings']],flush=True)
 if all(not q['crossings']and not q['coplanar']and not q['truncated']for q in rows):chosen=row;break
atomic_json(out/'fit.json',dict(passed=chosen is not None,chosen=chosen,trials=trials,diaphragm_fixed=True,hepatic_geometry_fixed=True,elapsed_seconds=time.monotonic()-start))
if chosen is None:raise ValueError('No current course change clears the fixed adjoining tissues')
