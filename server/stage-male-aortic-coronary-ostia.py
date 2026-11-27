"""Bounded communicating aortic/LCA/proximal-RCA wall from actual source ports.

No heart or far coronary control moves. Temporary root classifiers must be
buried and removed; four external blood ports receive annular wall rims only.
"""
from pathlib import Path
import sys,time,json,hashlib
import numpy as np
import manifold3d as md
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,CommonPositionBank
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim,continuity_audit
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.audit import audit

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-coronary-ostia-v2';D.mkdir(exist_ok=False);cp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';source_bytes=cp.read_bytes();source_sha=hashlib.sha256(source_bytes).hexdigest();a=json.loads(source_bytes);rows={p['id']:p for p in a['parts']};frozen_path=ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/atlas.json';assert digest_file(frozen_path)=='a0d7b46e033d145617a8937e1aff489db90b7a2808aa76ec90192f155ad95518';frozen=json.loads(frozen_path.read_text());old={p['id']:p for p in frozen['parts']};owned={p['part']for p in json.loads(frozen_path.with_name('preview-proof.json').read_text())['physical_tissues']}|{'ZA:Ascending aorta','ZA:Left coronary artery','ZA:Right coronary artery'}
for pid in owned:
 assert {k:v for k,v in rows[pid].items()if k!='chunk'}=={k:v for k,v in old[pid].items()if k!='chunk'},pid
 assert a['chunks'][rows[pid]['chunk']]==frozen['chunks'][old[pid]['chunk']],pid
atomic_json(D/'input-rebase.json',dict(consumed_catalogue_sha256=source_sha,frozen_fitted_heart_catalogue_sha256=digest_file(frozen_path),all_owned_part_records_literal=True,all_owned_chunk_descriptors_literal=True,owned_part_ids=sorted(owned)))
reader=AtlasReader();frame='model:local-ta98:assembled-metres';parent=Mesh.load(ROOT/'.local-models/geometry-workbench/heart/male-aortic-section-transport-v1/candidate-outer.npz');sources=[parent,reader.read_part(a,rows['ZA:Left coronary artery'],frame=frame)];rca=reader.read_part(a,rows['ZA:Right coronary artery'],frame=frame);rb=max(recover_station_rings(rca.vertices,rca.faces),key=lambda b:max(rca.vertices[b['rings'][0],1].mean(),rca.vertices[b['rings'][-1],1].mean()));rcarings=list(rb['rings']);assert rca.vertices[rcarings[0],1].mean()>rca.vertices[rcarings[-1],1].mean();selected=np.concatenate(rcarings[:13]);keep=np.all(np.isin(rca.faces,selected),axis=1);sources.append(Mesh(rca.vertices,rca.faces[keep],frame));V=[];F=[];R=[];C=[];offset=0
for source in sources:
 used,inv=np.unique(source.faces,return_inverse=True);v=source.vertices[used];f=inv.reshape(-1,3);branch=recover_station_rings(v,f);assert len(branch)==1;rings=list(branch[0]['rings']);
 if len(V)==1 and v[rings[0],1].mean()<v[rings[-1],1].mean():rings=rings[::-1]
 if len(V)==2 and v[rings[0],1].mean()<v[rings[-1],1].mean():rings=rings[::-1]
 V.append(v);F.append(f);R.append(rings);C.append(np.array([v[r].mean(0)for r in rings]));offset+=len(v)
assert len(R[0])==25 and len(R[1])==25 and len(R[2])==13
thickness=[.0009,.0003,.0003];ORIGIN=np.array([0.,1.4,.03125]);external=[(0,0),(0,24),(1,24),(2,12)];targets={}
for i in [1,2]:
 station=int(np.argmin(np.linalg.norm(C[0]-C[i][0],axis=1)));center=C[0][station];targets[i]=center+.60*(C[i][0]-center)
plan=dict(schema='male-aortic-coronary-root-construction/v1',source_catalogue_sha256=source_sha,heart_source_unchanged=True,original_native_LCA_and_RCA_positions_unchanged=True,source_aortic_outer_repair='male-aortic-section-transport-v1',root_sleeves=[dict(child=i,native_source_root=C[i][0].tolist(),buried_classifier_center=targets[i].tolist(),length_m=float(np.linalg.norm(targets[i]-C[i][0])))for i in [1,2]],external_ports=[dict(component=i,station=s,controls=12)for i,s in external],illustrative_thickness_m=thickness,heart_and_arch_connections_pending=True,anatomical_acceptance=False,installed=False);atomic_json(D/'construction-plan.json',plan)
inputs={};sheets=[];ancestry_rows=[];closed_blood=None
for inner in [False,True]:
 solids=[];native_fields=[]
 for i in range(3):
  native=V[i].copy()
  if inner:
   for r,c in zip(R[i],C[i]):
    delta=V[i][r]-c;native[r]=V[i][r]-delta/np.linalg.norm(delta,axis=1)[:,None]*thickness[i]
  native_fields.append(native);v=native.copy();f=F[i].copy();original_faces=len(f);internal_loop=None
  if i in [1,2]:
   root=next(r for r in boundary_loops(f)if set(r)==set(R[i][0]));internal_loop=np.arange(len(v),len(v)+len(root));sleeve=targets[i]+.8*(v[root]-v[root].mean(0));v=np.r_[v,sleeve];f=np.r_[f,annular_rim(root,internal_loop)]
  external_caps=[];internal_caps=[]
  for loop in boundary_loops(f):
   center=len(v);v=np.r_[v,v[loop].mean(0)[None]];ids=np.arange(len(f),len(f)+len(loop));f=np.r_[f,np.array([[center,int(loop[(j+1)%len(loop)]),int(loop[j])]for j in range(len(loop))])];(internal_caps if internal_loop is not None and set(loop)==set(internal_loop)else external_caps).extend(ids.tolist())
  solid,oid=make_solid((v-ORIGIN)*1000,f);assert solid.status()==md.Error.NoError;inputs[oid]=dict(inner=inner,component=i,source_count=original_faces,external=external_caps,internal=internal_caps);solids.append(solid)
 joined=md.Manifold.batch_boolean(solids,md.OpType.Add);assert joined.status()==md.Error.NoError;raw=joined.to_mesh64();vp=np.asarray(raw.vert_properties)[:,:3]/1000+ORIGIN;f=np.asarray(raw.tri_verts,int);oid,originface=provenance(raw);internal=np.array([int(fi)in inputs[int(o)]['internal']for o,fi in zip(oid,originface)]);external_caps=np.array([int(fi)in inputs[int(o)]['external']for o,fi in zip(oid,originface)]);atomic_json(D/(('inner'if inner else'outer')+'-classifier-proof.json'),dict(internal_classifier_faces_retained=int(internal.sum()),external_classifier_faces=int(external_caps.sum())));assert not internal.any(),('unburied_internal_classifier',inner,int(internal.sum()))
 bank=CommonPositionBank(np.concatenate(native_fields),tolerance=2e-12);mapping=np.array([bank.insert(q)for q in vp]);sheet=Mesh(bank.positions,mapping[f[~external_caps]],frame);label='inner'if inner else'outer';sheet.save(D/(label+'-sheet.npz'));check=continuity_audit(sheet.vertices,sheet.faces,expected_mouths=4,expected_components=1);atomic_json(D/(label+'-continuity.json'),check);assert check['topological_continuity_pass'];sheets.append(sheet);ancestry_rows.append(dict(component=np.array([inputs[int(o)]['component']for o in oid[~external_caps]]),source_face=originface[~external_caps]));
 if inner:closed_blood=Mesh(bank.positions,mapping[f],frame);closed_blood.save(D/'closed-blood-classifier.npz')
outer,inner=sheets;offset=len(outer.vertices);v=np.r_[outer.vertices,inner.vertices];faces=[outer.faces,inner.faces[:,::-1]+offset];outer_loops=boundary_loops(outer.faces);inner_loops=boundary_loops(inner.faces);ports=[]
for component,station in external:
 expected=V[component][R[component][station]];center=expected.mean(0);ol=min(outer_loops,key=lambda r:np.linalg.norm(outer.vertices[r].mean(0)-center));il=min(inner_loops,key=lambda r:np.linalg.norm(inner.vertices[r].mean(0)-center));target=outer.vertices[ol]-(outer.vertices[ol]-center)/np.linalg.norm(outer.vertices[ol]-center,axis=1)[:,None]*thickness[component];order=np.argmin(np.linalg.norm(inner.vertices[il][None]-target[:,None],axis=2),axis=1);assert len(ol)==len(il)==len(set(order))==12;assert np.linalg.norm(inner.vertices[il[order]]-target,axis=1).max()<1e-11;assert all(np.any(np.all(outer.vertices[ol]==point,axis=1))for point in expected);rim=annular_rim(ol,il[order]+offset);ports.append(dict(component=component,station=station,outer_indices=ol.tolist(),inner_indices=(il[order]+offset).tolist(),annular_face_start=int(sum(map(len,faces))),annular_faces=len(rim),literal_native_outer_mouth=True));faces.append(rim)
wall=Mesh(v,np.concatenate(faces),frame);wall.save(D/'finite-wall.npz');check=audit(wall,require_closed=True,expected_components=1,expected_genus=[3]);atomic_json(D/'numerical-audit.json',check);assert check['passed'],check['failures'];atomic_npz(D/'surfaces-and-ancestry.npz',vertices=v,outer_faces=outer.faces,inner_faces=inner.faces[:,::-1]+offset,outer_source_component=ancestry_rows[0]['component'],inner_source_component=ancestry_rows[1]['component'],outer_source_face=ancestry_rows[0]['source_face'],inner_source_face=ancestry_rows[1]['source_face'],frame=np.asarray(frame),**{f'annulus_{i}_faces':r for i,r in enumerate(faces[2:])})
report=dict(schema='male-aortic-coronary-ostia-finite-network/v1',source_catalogue_sha256=source_sha,geometry_sha256=wall.geometry_sha,wall_faces=len(wall.faces),external_blood_ports=4,communicating_coronary_origins=2,internal_classifier_caps=0,filled_lumen_disks=0,finite_wall_passed=check['passed'],all_native_coronary_source_controls_literal=True,heart_unchanged=True,heart_and_arch_interfaces_pending=True,illustrative_thickness_m=thickness,ports=ports,elapsed_seconds=time.monotonic()-T,anatomical_acceptance=False,installable=False);atomic_json(D/'construction-proof.json',report);print(json.dumps({k:v for k,v in report.items()if k!='ports'}))
