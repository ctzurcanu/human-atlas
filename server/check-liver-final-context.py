"""Replay final reference witnesses and fixed stomach/diaphragm clearance."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.liver_completion import validate_liver_completion
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--candidate');args=ap.parse_args();sex=args.sex;b=ROOT/'.local-models/geometry-workbench/liver-completion/final-ready-v265';candidate=Path(args.candidate)if args.candidate else b/(sex+'.json');a=json.loads(candidate.read_text());before=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json').read_text());bp={p['id']:p for p in before['parts']};r=AtlasReader();frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres';start=time.monotonic();witnesses=validate_liver_completion(a,reader=r,frame=frame);host=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Liver:capsule');m=r.read_part(a,host,frame=frame);outer=Mesh(m.vertices,m.faces[:host['liverReferenceTriangleCount']],frame);bv=BVH(outer);rows=[]
for p in a['parts']:
 if p.get('suppressed')or p.get('defaultHidden')or p['id']in a.get('defaultHidden',[]):continue
 target=p['id']in a['finiteDiaphragmIntegration']['partIds']or p['id'].startswith('GASTRIC:'+sex+':019:')or p.get('system')=='skeletal'or(p['id'].startswith('OESOPHWALL:'+sex+':')and p['id'].endswith(':external-connective-coat'))
 if not target:continue
 if not(p['id'].startswith('GASTRIC:male:')and a['gastricCompleteWall'].get('sizeCorrection')):assert p==bp[p['id']],('A fixed oesophagus/skeleton/diaphragm record moved',p['id'])
 bb=np.asarray(p['bounds'])
 if np.any(bb[1]<outer.bounds[0])or np.any(bb[0]>outer.bounds[1]):continue
 peer=r.read_part(a,p,frame=frame)
 if p['id'].startswith('OESOPHWALL:'):peer=Mesh(peer.vertices,peer.faces[:p['tissueDrawGroups'][0]['count']//3],frame)
 _,q=pair_contacts(outer,peer,0,first_bvh=bv,max_pairs=200000,max_tests=4000000);rows.append(dict(part=p['id'],crossings=q['crossing_face_pairs'],coplanar=q['coplanar_overlap_face_pairs'],truncated=q['truncated']));print(p['id'],q['crossing_face_pairs'],flush=True)
passed=bool(rows)and all(not x['crossings']and not x['coplanar']and not x['truncated']for x in rows);atomic_json(candidate.parent/(sex+'-context.json'),dict(passed=passed,checks=rows,reference_witnesses=len(witnesses),fixed_skeleton_oesophagus_diaphragm=True,male_stomach_size_adjustment=bool(a['gastricCompleteWall'].get('sizeCorrection')),elapsed_seconds=time.monotonic()-start));print(sex,passed,'context checks',len(rows),flush=True)
if not passed:raise ValueError('Final context check failed; do not install')
