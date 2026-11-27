"""Fit the liver impression around the fixed oesophagus, preserving native mouths."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.liver_impression_field import transform,posterior_impression,slab_impression
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--posterior',action='store_true');ap.add_argument('--down-probe',action='store_true');ap.add_argument('--slab',action='store_true');ap.add_argument('--stage');args=ap.parse_args();sex=args.sex;b=ROOT/'.local-models/geometry-workbench/liver-completion/final-ready-v265';out=b/(args.stage or sex+('-down-probe'if args.down_probe else'-posterior-probe'if args.posterior else'-impression-probe'));out.mkdir(exist_ok=False);a=json.loads((b/(sex+'.json')).read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres';S=a['liverTransferIntegration']['scale']if sex=='female'else 1.;T=np.array(a['liverTransferIntegration']['translation_m'])if sex=='female'else np.zeros(3);start=time.monotonic();p=parts['TA98PART:ZA:Liver:capsule'];m=r.read_part(a,p,frame=frame);source=(m.vertices-T)/S;faces=m.faces[:p['liverReferenceTriangleCount']];targets=[]
ids=list(a['finiteDiaphragmIntegration']['partIds'])+[p['id']for p in a['parts']if p['id'].startswith('GASTRIC:'+sex+':019:')]+[p['id']for p in a['parts']if p['id'].startswith('OESOPHWALL:'+sex+':')and p['id'].endswith(':external-connective-coat')]
for pid in ids:
 p=parts[pid];peer=r.read_part(a,p,frame=frame)
 if pid.startswith('OESOPHWALL:'):peer=Mesh(peer.vertices,peer.faces[:p['tissueDrawGroups'][0]['count']//3],frame)
 targets.append((pid,peer,BVH(peer)))
trials=[];chosen=None
for amplitude,down in ([(.012 if sex=='male'else .024,x)for x in [.003,.006,.009,.012,.015]]if args.slab else[(.012,x)for x in [.003,.006,.009,.012]]if args.down_probe else[(x,0.)for x in ([.006,.012,.018,.024,.030]if args.posterior else[.008,.016,.024,.032,.040])]):
 v,J=(slab_impression(source,amplitude,down=down)if args.slab else posterior_impression(source,amplitude,down=down)if args.posterior or args.down_probe else transform(source,amplitude));outer=Mesh(v*S+T,faces,frame);bv=BVH(outer);rows=[]
 for pid,peer,pb in targets:
  if np.any(peer.bounds[1]<outer.bounds[0])or np.any(peer.bounds[0]>outer.bounds[1]):continue
  _,q=pair_contacts(outer,peer,0,first_bvh=bv,second_bvh=pb,max_pairs=200000,max_tests=5000000);rows.append(dict(part=pid,crossings=q['crossing_face_pairs'],coplanar=q['coplanar_overlap_face_pairs'],truncated=q['truncated']))
 row=dict(amplitude_m=amplitude,down_m=down,checks=rows,minimum_jacobian=float(np.linalg.det(J).min()),maximum_displacement_m=float(np.linalg.norm(v-source,axis=1).max()));trials.append(row);print(sex,amplitude,down,[(q['part'],q['crossings'])for q in rows if q['crossings']],flush=True)
 if all(not q['crossings']and not q['coplanar']and not q['truncated']for q in rows):chosen=row;break
atomic_json(out/'fit.json',dict(passed=chosen is not None,chosen=chosen,trials=trials,fixed=['oesophagus','diaphragm','stomach','skeleton'],source_field='geometry_tools/liver_impression_field.py',field_kind='posterior'if args.posterior or args.down_probe else'lateral',elapsed_seconds=time.monotonic()-start))
if chosen is None:raise ValueError('Bounded impression field did not clear all contacts; not installed')
