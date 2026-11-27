"""Audit displayed chest tissues against literal installed/candidate heart parts.

The catalogue-level hidden donor list is authoritative. Results retain both
part identities; a coat-only pass is never labelled a whole-heart pass.
This is a surface-crossing check, not a containment or anatomical certificate.
"""
from pathlib import Path
import argparse, json, sys, time
from collections import Counter
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_json, atomic_npz, digest_file

ap=argparse.ArgumentParser()
ap.add_argument('--catalogue',default='.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json')
ap.add_argument('--stage',required=True)
ap.add_argument('--part-set',choices=['exposed','walls','roots','all','other-cardiac','pulmonary','packing-vascular'],default='exposed')
ap.add_argument('--context',choices=['chest-wall','lungs','bronchi','native-lungs','oesophagus','diaphragm'],required=True)
args=ap.parse_args();start=time.monotonic()
path=ROOT/args.catalogue;atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']}
marker=atlas['finiteCardiacIntegration'];sex=marker['sex'];hidden=set(atlas.get('defaultHidden',[]))
all_ids=set(sum([marker.get(k,[])for k in ['wallPartIds','valvePartIds','papillaryPartIds','chordalPartIds','rootPartIds']],[]))
if args.part_set=='exposed':ids=[f'TA98HEART:{sex}:epicardium',*marker['rootPartIds']]
elif args.part_set=='walls':ids=marker['wallPartIds']
elif args.part_set=='roots':ids=marker['rootPartIds']
elif args.part_set=='other-cardiac':ids=[p['id']for p in atlas['parts']if p['system']=='cardiac'and p['id']not in all_ids and p['id']not in hidden and not p.get('suppressed')and not p.get('defaultHidden')]
elif args.part_set=='pulmonary':ids=atlas['pulmonarySegmentIntegration']['partIds']+atlas['pulmonaryBranchIntegration']['partIds']
elif args.part_set=='packing-vascular':ids=sorted(set(marker['rootPartIds'])|{p['id']for p in atlas['parts']if p.get('localRootPacking',{}).get('field')=='aortic-arch-arterial-cohort' or p.get('cardiacRootRegistration')})
else:ids=sorted(all_ids)
assert all(not parts[i].get('suppressed')for i in ids)
reader=AtlasReader();frame=f'model:local-{sex}-ta98:assembled-metres'
def read(pid):return reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=path)
meshes=[read(i)for i in ids]
if not meshes:
 out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
 report=dict(catalogue_sha256=digest_file(path),part_set=args.part_set,context=args.context,moving_ids=[],empty_moving_set=True,complete_within_pair_budget=True,crossing_face_pairs=0,coplanar_overlap_face_pairs=0,part_pairs=[],elapsed_seconds=time.monotonic()-start)
 atomic_json(out/'report.json',report);print('No additional active cardiac source meshes; explicit empty-set audit.',flush=True);raise SystemExit(0)
lo=np.min([m.bounds[0]for m in meshes],axis=0);hi=np.max([m.bounds[1]for m in meshes],axis=0)
def active(p):return p['id']not in hidden and not p.get('defaultHidden')and not p.get('suppressed')
if args.context=='chest-wall':
 fixed_ids=[p['id']for p in atlas['parts']if active(p)and p['system']in ['skeletal','muscular']and p['id']not in all_ids and np.all(np.array(p['bounds'])[1]>=lo)and np.all(np.array(p['bounds'])[0]<=hi)]
elif args.context=='lungs':fixed_ids=atlas['pulmonarySegmentIntegration']['partIds']
elif args.context=='bronchi':fixed_ids=atlas['pulmonaryBranchIntegration']['partIds']
elif args.context=='native-lungs':fixed_ids=['LOCAL:female:l_lung:1034','LOCAL:female:r_lung:1044']if sex=='female'else ['LOCAL:male:l_lung:802','LOCAL:male:r_lung:808']
elif args.context=='diaphragm':fixed_ids=atlas['finiteDiaphragmIntegration']['partIds']
else:fixed_ids=atlas['oesophagealWallGeometry']['partIds']
assert fixed_ids
out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
def combine(ms):
 offsets=np.r_[0,np.cumsum([len(m.vertices)for m in ms])]
 return Mesh(np.concatenate([m.vertices for m in ms]),np.concatenate([m.faces+offsets[i]for i,m in enumerate(ms)]),frame),np.repeat(np.arange(len(ms)),[len(m.faces)for m in ms])
moving,mowners=combine(meshes);fixed,fowners=combine([read(i)for i in fixed_ids])
atomic_json(out/'inputs.json',dict(catalogue_sha256=digest_file(path),part_set=args.part_set,context=args.context,moving_ids=ids,fixed_ids=fixed_ids,hidden_donors_excluded=True,moving_faces=len(moving.faces),fixed_faces=len(fixed.faces)))
result,proof=pair_contacts(moving,fixed,0.,max_pairs=500000,max_tests=4000000)
atomic_npz(out/'contacts.npz',**result)
counts=Counter((ids[mowners[i]],fixed_ids[fowners[j]])for(i,j),kind in zip(result['pairs'],result['relation'])if kind in ['crossing','coplanar-overlap'])
report=dict(**proof,part_pairs=[dict(heart=a,context=b,name=parts[b]['name'],crossings=n)for(a,b),n in counts.most_common()],catalogue_sha256=digest_file(path),elapsed_seconds=time.monotonic()-start)
atomic_json(out/'report.json',report)
summary={k:v for k,v in report.items()if k in ['crossing_face_pairs','coplanar_overlap_face_pairs','complete_within_pair_budget','elapsed_seconds']};summary['part_pair_count']=len(report['part_pairs']);summary['first_part_pairs']=report['part_pairs'][:12]
print(json.dumps(summary,indent=2),flush=True)
