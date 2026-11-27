"""Apply one declared proper similarity to the complete owned cardiac closure."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
from scipy.spatial.transform import Rotation
from importlib import util
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.geometry_instance import instance_pose
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--proposal',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--catalogue');args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
path=ROOT/(args.catalogue or f'.local-models/ta98-review/bowel-wall-local-atlas/{args.sex}/atlas.json');original=json.loads(path.read_text());a=copy.deepcopy(original);proposal=json.loads((h/args.proposal/'proposal.json').read_text());assert proposal['source_catalogue_sha256']==digest_file(path)
mat=np.array(proposal['matrix']);scale=float(np.cbrt(np.linalg.det(mat[:3,:3])));r=mat[:3,:3]/scale;assert scale>0 and np.max(abs(r.T@r-np.eye(3)))<1e-12 and np.linalg.det(r)>0
spec=util.spec_from_file_location('heart_forest',ROOT/'scripts/build-spatial-dependency-forest.py');module=util.module_from_spec(spec);spec.loader.exec_module(module);forest=module.build_forest(args.sex,original);groups=[g for root in forest['forest']for g in root['children']if g['id']=='assembly:heart'];assert len(groups)==1;moving=groups[0]['activePartIds'];forest['sourceSha256']=digest_file(path);atomic_json(out/'source-forest.json',forest)
parts={p['id']:p for p in a['parts']};old={p['id']:p for p in original['parts']};reader=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres'
def point(v):return (np.array(v)@mat[:3,:3].T+mat[:3,3]).tolist()
for pid in moving:
 p=parts[pid];prior=instance_pose(p);old_r,old_scale=prior if prior is not None else (np.eye(3),1.)
 p['geometryOrigin']=point(p.get('geometryOrigin',[0,0,0]));p['geometryInstance']=dict(rotation=Rotation.from_matrix(r@old_r).as_quat().tolist(),scale=old_scale*scale)
 p['cardiacThoracicFit']=dict(schema='cardiac-thoracic-packing/v2',sourceCatalogueSha256=digest_file(path),matrix=mat.tolist())
 mesh=reader.read_part(a,p,frame=frame,catalogue_path=path);source=reader.read_part(original,old[pid],frame=frame,catalogue_path=path)
 assert np.max(abs(mesh.vertices-(source.vertices@mat[:3,:3].T+mat[:3,3])))<2e-15
 p['bounds']=mesh.bounds.tolist()
for g in a.get('cardiacFeatureIntegration',{}).get('termGroups',[]):g['anchor']['point']=point(g['anchor']['point'])
marker=a['finiteCardiacIntegration'];marker['geometryInstance']=copy.deepcopy(parts[marker['wallPartIds'][0]]['geometryInstance']);marker['thoracicFit']=dict(schema='cardiac-thoracic-packing/v2',sourceCatalogueSha256=digest_file(path),matrix=mat.tolist(),movingPartIds=moving,independentRecipientJoinsVerified=False,clinicalSizeEstablished=False)
for pid,p in parts.items():
 if pid not in moving:assert p==old[pid]
assert original['chunks']==a['chunks']
atomic_json(out/'atlas.json',a);atomic_json(out/'packing.json',dict(source_catalogue_sha256=digest_file(path),moving_parts=moving,relative_scale=scale,matrix=mat.tolist(),all_nonmembers_literal=True,banks_modified=False,internal_coordinates_and_interfaces_preserved=True,installed=False,elapsed_seconds=time.monotonic()-start));print('Prepared',args.sex,len(moving),'owned cardiac parts; fixed skeleton and recipient tissues unchanged',flush=True)
