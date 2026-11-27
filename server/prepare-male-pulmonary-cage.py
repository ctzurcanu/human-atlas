"""Explicit lower-thorax cage; upper lungs and airway intake stay literal.

One positive-Jacobian field carries all 19 territories and 27 airway pieces.
Only a diagnostic candidate is emitted. Actual tissue/contact gates follow.
"""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
from scipy.interpolate import PchipInterpolator
from importlib import util
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.algorithms.longitudinal_organ_cage import LongitudinalOrganCage
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--width',type=float,required=True);ap.add_argument('--depth',type=float,required=True);ap.add_argument('--dx-mm',type=float,default=0);ap.add_argument('--dz-mm',type=float,default=0);ap.add_argument('--middle-lift-mm',type=float,default=0);args=ap.parse_args();start=time.monotonic()
assert 1<=args.width<=1.2 and .9<=args.depth<=1.1 and abs(args.dx_mm)<=5 and abs(args.dz_mm)<=5
assert 0<=args.middle_lift_mm<=12
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';original=json.loads(path.read_text());a=copy.deepcopy(original);parts={p['id']:p for p in a['parts']};old={p['id']:p for p in original['parts']};reader=AtlasReader();frame='model:local-male-ta98:assembled-metres'
spec=util.spec_from_file_location('pulmonary_cage_forest',ROOT/'scripts/build-spatial-dependency-forest.py');module=util.module_from_spec(spec);spec.loader.exec_module(module);forest=module.build_forest('male',original);g=[g for root in forest['forest']for g in root['children']if g['id']=='assembly:pulmonary'];assert len(g)==1;moving=g[0]['activePartIds'];forest['sourceSha256']=digest_file(path);atomic_json(out/'source-forest.json',forest)
fractions=np.array([0.,.05,.95,1.]);source_cov=np.broadcast_to(np.diag([2.,1.]),(4,2,2)).copy();target_cov=source_cov.copy();target_cov[:2]=np.diag([2*args.width**2,args.depth**2]);centres=np.zeros((4,2));targets=centres.copy();targets[:2]=[args.dx_mm/1000,args.dz_mm/1000]
cage=LongitudinalOrganCage(np.array([1.42,1.475]),np.array([1.42,1.475]),fractions,centres,targets,source_cov,target_cov,minimum_stretch=.9)
ys=np.array([1.28,1.30,1.32,1.36,1.40,1.44,1.46,1.475,1.50]);yt=ys+args.middle_lift_mm/1000*np.array([0,0,0,.8,1.,.4,0,0,0]);ywarp=PchipInterpolator(ys,yt);yd=ywarp.derivative();assert yd(np.linspace(1.30,1.475,4096)).min()>.5
buffers=[];offset=0;layout={};records=[];chunk=len(a['chunks']);origin=np.array([0.,1.3521,0.]);bank=out/'pulmonary-cage.bin'
def put(arr):
 global offset
 pad=(-offset)%4
 if pad:buffers.append(bytes(pad));offset+=pad
 value=offset;data=arr.tobytes();buffers.append(data);offset+=len(data);return value
for pid in moving:
 p=parts[pid];assert p.get('geometryOrigin')==origin.tolist()and not p.get('geometryInstance')and not p.get('sourceOffset')and p.get('geometryScale',1)==1
 mesh=reader.read_part(original,old[pid],frame=frame,catalogue_path=path);q=mesh.vertices[:,[1,0,2]];mapped_yxz=cage.map(q);active=(q[:,0]>1.30)&(q[:,0]<1.475);mapped_yxz[active,0]=ywarp(q[active,0]);j_yxz=cage.jacobians(q);j_yxz[active,0,0]=yd(q[active,0]);mapped=mapped_yxz[:,[1,0,2]];j=j_yxz[:,[1,0,2]][:,:,[1,0,2]];det=np.linalg.det(j);assert det.min()>0
 positions=(mapped-origin).astype('<f4');actual=positions.astype(float)+origin;upper=mesh.vertices[:,1]>=1.475;assert np.array_equal(actual[upper],mesh.vertices[upper]),'Upper airway/source controls moved'
 row=dict(positions=put(positions));attrs=mesh.attributes
 for k,value in attrs.items():
  if k=='normals':
   normals=np.linalg.solve(j.transpose(0,2,1),value.astype(float)[...,None])[...,0];length=np.linalg.norm(normals,axis=1);assert np.all(length>0);normals/=length[:,None];words=np.rint(normals*32767).astype('<i2');words[upper]=value[upper].astype('<i2');value=words
  elif k=='uvs':value=value.astype('<f4')
  elif k=='normalSigns':value=value.astype('i1')
  elif k=='colors':value=value.astype('u1')
  else:raise ValueError(k)
  row[k]=put(value)
 row['indices']=put(mesh.faces.astype('<u4'));p.update(row,chunk=chunk,vertexCount=len(mesh.vertices),indexCount=int(mesh.faces.size),triangles=len(mesh.faces),bounds=np.array([actual.min(0),actual.max(0)]).tolist());p['pulmonaryThoracicCage']=dict(schema='longitudinal-pulmonary-packing/v1',sourceCatalogueSha256=digest_file(path),upperFixedYMetres=1.475,fullGrowthBelowYMetres=1.42,lowerWidthScale=args.width,lowerDepthScale=args.depth,lowerTranslationXZMetres=[args.dx_mm/1000,args.dz_mm/1000],middleLiftMetres=args.middle_lift_mm/1000,minimumContinuousJacobian=float(det.min()),upperSourceControlsLiteral=True)
 records.append(dict(id=pid,minimum_jacobian=float(det.min()),maximum_storage_error_m=float(np.abs(actual-mapped).max()),upper_unchanged_controls=int(upper.sum())))
atomic_bytes(bank,b''.join(buffers));a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size))
for p in a['parts']:
 if p['id']not in moving:assert p==old[p['id']]
for k in ['pulmonarySegmentIntegration','pulmonaryBranchIntegration']:a[k]['thoracicCage']=dict(schema='longitudinal-pulmonary-packing/v1',bankSha256=digest_file(bank),all46PartsInOneField=True,upperIntakeUnchanged=True,anatomicalAcceptance=False)
atomic_json(out/'atlas.json',a);atomic_json(out/'packing.json',dict(source_catalogue_sha256=digest_file(path),moving_parts=moving,all_nonmembers_literal=True,source_banks_literal=True,replacement_bank=str(bank),replacement_bank_sha256=digest_file(bank),records=records,upper_fixed_y_m=1.475,lower_full_y_m=1.42,width_scale=args.width,depth_scale=args.depth,translation_xz_m=[args.dx_mm/1000,args.dz_mm/1000],middle_lift_m=args.middle_lift_mm/1000,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start));print('Prepared',len(moving),'coordinated pulmonary parts; upper airway literal; positive cage Jacobian',flush=True)
