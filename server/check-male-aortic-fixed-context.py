"""Bounded grouped contacts against an explicit immutable context catalogue."""
from pathlib import Path
import sys,json,time,argparse,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts

args=argparse.ArgumentParser();args.add_argument('--catalogue',type=Path,required=True);args.add_argument('--group',choices=['heart','lungs','skeleton-muscles','viscera','vessels','nerves'],required=True);args=args.parse_args();T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-coronary-ostia-v2';O=D/'fixed-context-v1';O.mkdir(exist_ok=True);out=O/(args.group+'.json');assert not out.exists();a=json.loads(args.catalogue.read_text());source_sha=digest_file(args.catalogue);reader=AtlasReader();manifest=json.loads((D/'owned-bank-v1/tissue.json').read_text());tri=[]
for name,row in manifest['parts'].items():
 m=Mesh.load(D/'owned-bank-v1'/(row['id'].replace(':','_')+'.npz'));keep=np.concatenate([np.arange(g['start']//3,(g['start']+g['count'])//3)for g in row['tissueDrawGroups']if not g.get('neighbor')]);tri.append(m.triangles[keep])
tri=np.concatenate(tri);v,iv=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);candidate=Mesh(v,iv.reshape(-1,3),'model:local-ta98:assembled-metres')
if args.group=='heart':candidate.save(O/'actual-decoded-assembled-wall.npz')
physical={r['part']for r in json.loads((ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/preview-proof.json').read_text())['physical_tissues']};own={'ZA:Ascending aorta','ZA:Left coronary artery','ZA:Right coronary artery'}
def wanted(p):
 pid=p['id'];system=p.get('system');primary=pid.startswith('ZA:');
 if args.group=='heart':return pid in physical
 if args.group=='lungs':return pid.startswith('TA98PULMONARY:male:')and system=='respiratory'
 if args.group=='skeleton-muscles':return primary and system in ['skeletal','muscular']
 if args.group=='viscera':return (primary and system in ['respiratory','digestive'])or pid=='OESOPHWALL:male:native:external-connective-coat'
 if args.group=='vessels':return primary and system in ['arterial','venous']and pid not in own
 if args.group=='nerves':return primary and system=='nervous'
contexts=[]
for p in a['parts']:
 if p.get('suppressed')or not wanted(p):continue
 bounds=np.asarray(p['bounds'])
 if np.any(bounds[1]<candidate.bounds[0])or np.any(bounds[0]>candidate.bounds[1]):continue
 contexts.append(p)
records=[]
for p in contexts:
 m=reader.read_part(a,p,frame=candidate.frame)
 if args.group=='heart':
  groups=[np.arange(g['start']//3,(g['start']+g['count'])//3)for g in p.get('tissueDrawGroups',[])if g.get('neighbor')not in physical]
  if not groups:continue
  m=Mesh(m.vertices,m.faces[np.concatenate(groups)],m.frame)
 pairs,summary=pair_contacts(candidate,m,0.,max_pairs=100000,max_tests=1000000);wrong=np.isin(pairs['relation'],['crossing','coplanar-overlap']);example=pairs['pairs'][wrong][:12];records.append(dict(part=p['id'],context_geometry_sha256=m.geometry_sha,source_part_sha256=hashlib.sha256(json.dumps(p,sort_keys=True,separators=(',',':')).encode()).hexdigest(),summary=summary,example_face_pairs=example.tolist()));assert summary['complete_within_pair_budget']
report=dict(schema='male-aortic-actual-decoded-fixed-context/v1',source_catalogue=str(args.catalogue.resolve()),source_catalogue_sha256=source_sha,owned_bank_sha256=manifest['bank_sha256'],candidate_decoded_geometry_sha256=candidate.geometry_sha,group=args.group,context_records=records,context_meshes_checked=len(records),crossing_parts=[r['part']for r in records if r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs']],anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(out,report);print(json.dumps(dict(group=args.group,elapsed_seconds=report['elapsed_seconds'],context_meshes_checked=len(records),crossings=[dict(part=r['part'],crossings=r['summary']['crossing_face_pairs'],overlaps=r['summary']['coplanar_overlap_face_pairs'])for r in records if r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs']])))
