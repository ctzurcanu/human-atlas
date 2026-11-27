"""Exact source-face contact delta for the one bounded aortic outer repair."""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import ROOT,atomic_json,digest_file,atomic_npz
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

args=argparse.ArgumentParser();args.add_argument('--catalogue',type=Path,required=True);args=args.parse_args();T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-coronary-ostia-v2';O=D/'contact-delta-v1';O.mkdir(exist_ok=False);a=json.loads(args.catalogue.read_text());rows={p['id']:p for p in a['parts']};reader=AtlasReader();before=Mesh.load(ROOT/'.local-models/geometry-workbench/heart/male-native-aortic-source-defects-v2/native-ascending-aorta.npz');after=Mesh.load(ROOT/'.local-models/geometry-workbench/heart/male-aortic-section-transport-v1/candidate-outer.npz');assert np.array_equal(before.faces,after.faces);physical={r['part']for r in json.loads((ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/preview-proof.json').read_text())['physical_tissues']};contexts=[]
for group in ['heart','vessels']:
 source=json.loads((D/'fixed-context-v1'/(group+'.json')).read_text());contexts.extend(rows[r['part']]for r in source['context_records'])
proof=[];arrays={};totals=dict(before_pairs=0,after_pairs=0,introduced_pairs=0,removed_pairs=0)
for i,p in enumerate(contexts):
 m=reader.read_part(a,p,frame=before.frame);keep=np.arange(len(m.faces))
 if p['id']in physical:
  groups=[np.arange(g['start']//3,(g['start']+g['count'])//3)for g in p.get('tissueDrawGroups',[])if g.get('neighbor')not in physical]
  if not groups:continue
  keep=np.concatenate(groups);m=Mesh(m.vertices,m.faces[keep],m.frame)
 tree=BVH(m);sets=[]
 for label,candidate in [('before',before),('after',after)]:
  actual,summary=pair_contacts(candidate,m,0.,second_bvh=tree,max_pairs=100000,max_tests=1000000);assert summary['complete_within_pair_budget'];wrong=np.isin(actual['relation'],['crossing','coplanar-overlap']);pairs=actual['pairs'][wrong].copy();pairs[:,1]=keep[pairs[:,1]];sets.append(set(map(tuple,pairs)));arrays[f'context_{i}_{label}_pairs']=pairs
 old,new=sets;introduced=new-old;removed=old-new;totals['before_pairs']+=len(old);totals['after_pairs']+=len(new);totals['introduced_pairs']+=len(introduced);totals['removed_pairs']+=len(removed);proof.append(dict(part=p['id'],before_source_face_pairs=len(old),after_source_face_pairs=len(new),introduced_source_face_pairs=len(introduced),removed_source_face_pairs=len(removed),introduced_context_faces=sorted({int(k[1])for k in introduced}),introduced_aortic_faces=sorted({int(k[0])for k in introduced}),examples_introduced=[list(map(int,k))for k in sorted(introduced)[:8]]))
report=dict(schema='male-aortic-exact-source-face-contact-delta/v1',source_context_sha256=digest_file(args.catalogue),original_aortic_geometry_sha256=before.geometry_sha,repaired_aortic_geometry_sha256=after.geometry_sha,unchanged_source_topology=True,totals=totals,parts=proof,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(O/'evidence.json',report);atomic_npz(O/'exact-source-face-pairs.npz',**arrays);print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'],totals=totals,nonzero=[{k:r[k]for k in ['part','before_source_face_pairs','after_source_face_pairs','introduced_source_face_pairs','removed_source_face_pairs']}for r in proof if r['before_source_face_pairs']or r['after_source_face_pairs']])))
