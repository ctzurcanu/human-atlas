"""Prepare a source-guided common pulmonary similarity, without installation."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
from scipy.spatial.transform import Rotation
from importlib import util
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.geometry_instance import instance_pose
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);choice=ap.add_mutually_exclusive_group(required=True);choice.add_argument('--factor',type=float);choice.add_argument('--proposal');args=ap.parse_args();start=time.monotonic()
if args.factor is not None:assert 1<=args.factor<=1.45
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';original=json.loads(path.read_text());a=copy.deepcopy(original);parts={p['id']:p for p in a['parts']};source_parts={p['id']:p for p in original['parts']};reader=AtlasReader();frame='model:local-male-ta98:assembled-metres'
def read(pid):return reader.read_part(original,source_parts[pid],frame=frame,catalogue_path=path)
# Match the centres and superior extent of the retained native male lungs.
# Native mesh data is a placement reference, never replacement geometry.
source=[];target=[]
for side,native in [('left','LOCAL:male:l_lung:802'),('right','LOCAL:male:r_lung:808')]:
 ids=[i for i in a['pulmonarySegmentIntegration']['partIds']if i.startswith('TA98PULMONARY:male:'+side+'-')];v=np.unique(np.concatenate([read(i).vertices for i in ids]),axis=0);n=read(native).vertices
 source.append((v.min(0)+v.max(0))*.5);target.append((n.min(0)+n.max(0))*.5)
source=np.array(source);target=np.array(target);pivot=source.mean(0)
# Retain the source orientation: only common positive scale and translation.
if args.proposal:
 proposal=json.loads((ROOT/args.proposal).read_text());assert proposal['source_catalogue_sha256']==digest_file(path);mat=np.array(proposal['lung_matrix']);scale=float(np.cbrt(np.linalg.det(mat[:3,:3])));assert 1<=scale<=1.45 and np.max(abs(mat[:3,:3]-np.eye(3)*scale))<1e-12;delta=mat[:3,3]-pivot*(1-scale)
else:
 scale=args.factor;delta=target.mean(0)-pivot;mat=np.eye(4);mat[:3,:3]*=scale;mat[:3,3]=pivot*(1-scale)+delta
spec=util.spec_from_file_location('pulmonary_forest',ROOT/'scripts/build-spatial-dependency-forest.py');module=util.module_from_spec(spec);spec.loader.exec_module(module);forest=module.build_forest('male',original);g=[g for root in forest['forest']for g in root['children']if g['id']=='assembly:pulmonary'];assert len(g)==1;moving=g[0]['activePartIds'];forest['sourceSha256']=digest_file(path);atomic_json(out/'source-forest.json',forest)
for pid in moving:
 p=parts[pid];prior=instance_pose(p);r,s=prior if prior is not None else(np.eye(3),1.);p['geometryOrigin']=(np.array(p.get('geometryOrigin',[0,0,0]))*scale+mat[:3,3]).tolist();p['geometryInstance']=dict(rotation=Rotation.from_matrix(r).as_quat().tolist(),scale=s*scale)
 actual=reader.read_part(a,p,frame=frame,catalogue_path=path);assert np.max(abs(actual.vertices-(read(pid).vertices*scale+mat[:3,3])))<2e-15;p['bounds']=actual.bounds.tolist()
 p['pulmonaryThoracicFit']=dict(schema='pulmonary-thoracic-packing/v1',sourceCatalogueSha256=digest_file(path),factor=scale,matrix=mat.tolist())
for marker in ['pulmonarySegmentIntegration','pulmonaryBranchIntegration']:a[marker]['thoracicFit']=dict(schema='pulmonary-thoracic-packing/v1',factor=scale,matrix=mat.tolist(),sourceCatalogueSha256=digest_file(path),nativeWholeLungsReferenceOnly=True)
for p,q in zip(original['parts'],a['parts']):
 if p['id']not in moving:assert p==q
atomic_json(out/'atlas.json',a);atomic_json(out/'packing.json',dict(source_catalogue_sha256=digest_file(path),moving_parts=moving,matrix=mat.tolist(),relative_scale=scale,translation_m=delta.tolist(),native_reference_centres=target.tolist(),source_centres=source.tolist(),banks_modified=False,all_nonmembers_literal=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start));print('Prepared',len(moving),'pulmonary parts; scale',scale,'translation mm',delta*1000,flush=True)
