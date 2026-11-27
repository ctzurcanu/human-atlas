"""Saved bounded stages: all preferred current geometry, two actual root AABBs."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';C=P/'three-course-finite-educational-v1';OUT=C/'current-context-v1';FRAME='model:local-ta98:assembled-metres';MASTER=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';RESOLVED=B/'proposal-v2/resolved-installed.json'

def prepare():
 start=time.monotonic();OUT.mkdir(exist_ok=False);raw=MASTER.read_bytes();atlas=json.loads(raw);availability=json.loads(RESOLVED.read_bytes());assert digest_file(MASTER)==availability['sourceSha256'];(OUT/'source-atlas.json').write_bytes(raw);preferred=[r['id']for r in availability['effectivePartVisibility']if r['available']and not r['suppressed']];excluded=['ZA:Middle cerebral artery (M3 segment).l','ZA:Postcentral arterial branch.l'];preferred=sorted(set(preferred)-set(excluded));parts={p['id']:p for p in atlas['parts']};preferred=[pid for pid in preferred if parts[pid].get('indexCount',0)>0];z=np.load(F/'three-course-source-charts.npz');g=json.loads((C/'generation-complete.json').read_text());q=np.load(C/'assembly-v1/wall-with-source-and-derived-corner-normals.npz');rows={r['original_id']:r for r in g['classifiers']};origin=np.array([0.,1.75,0.]);emitted=(q['vertices']-origin).astype('<f4').astype(float)+origin;domains=[]
 for index,allowed in enumerate([{'M3':set(z['M3_root_patch_face_ids'].tolist()),'main':set(z['main_root_to_M3_patch_face_ids'].tolist())},{'main':set(z['main_secondary_patch_face_ids'].tolist()),'twig':set(z['twig_root_patch_face_ids'].tolist())}]):
  root=OUT/('root'+str(index));root.mkdir();oldv=np.r_[z['M3_vertices'],z['branch_vertices']];oldf=[];mask=[];outer=[]
  for course,ids in allowed.items():oldf.extend(z['M3_faces'][sorted(ids)]if course=='M3'else z['branch_faces'][sorted(ids)]+len(z['M3_vertices']))
  for i,(oid,fid)in enumerate(zip(q['original_id'],q['source_face_id'])):
   row=rows[int(oid)]
   if row['course']not in allowed or int(fid)>=row['source_faces']:continue
   if row['source_face_ids'][int(fid)]not in allowed[row['course']]:continue
   mask.append(i)
   if row['role'].endswith('outer'):outer.append(i)
  old=Mesh(oldv,np.asarray(oldf),FRAME);newouter=Mesh(emitted,q['faces'][outer],FRAME);newwall=Mesh(emitted,q['faces'][mask],FRAME);old.save(root/'old-local-outer.npz');newouter.save(root/'emitted-local-outer.npz');newwall.save(root/'emitted-local-wall.npz');points=np.r_[old.triangles.reshape(-1,3),newwall.triangles.reshape(-1,3)];bounds=[(points.min(0)-1e-6).tolist(),(points.max(0)+1e-6).tolist()];domains.append(dict(index=index,bounds_m=bounds,old_outer_faces=len(old.faces),new_outer_faces=len(newouter.faces),new_wall_faces=len(newwall.faces),numerical_broadphase_margin_m=1e-6))
 batches=[preferred[i:i+256]for i in range(0,len(preferred),256)];plan=dict(schema='current-preferred-two-root-context-scan/v1',source_sha256=availability['sourceSha256'],resolved_source=str(RESOLVED),resolved_sha256=digest_file(RESOLVED),preferred_meshes=len(preferred),excluded_replaced_source_ids=excluded,visibility_rule='available && !suppressed in current frontend export; hidden/retired unavailable alternatives excluded. The 46 available default-hidden hepatic layers are conservatively included too.',batches=batches,domains=domains,frame=FRAME,elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'scan-plan.json',plan);print(json.dumps(dict(preferred_meshes=len(preferred),batches=len(batches),domains=domains,elapsed_seconds=plan['elapsed_seconds'])))

def scan(index):
 start=time.monotonic();plan=json.loads((OUT/'scan-plan.json').read_text());a=json.loads((OUT/'source-atlas.json').read_bytes());assert digest_file(OUT/'source-atlas.json')==plan['source_sha256'];parts={p['id']:p for p in a['parts']};reader=AtlasReader();hits=[];count=0
 for pid in plan['batches'][index]:
  mesh=reader.read_part(a,parts[pid],frame=FRAME);tri=mesh.triangles;lo=tri.min(1);hi=tri.max(1);count+=1
  for domain in plan['domains']:
   bounds=np.array(domain['bounds_m']);selected=np.flatnonzero(np.all(hi>=bounds[0],axis=1)&np.all(lo<=bounds[1],axis=1))
   if not len(selected):continue
   used,inv=np.unique(mesh.faces[selected],return_inverse=True);local=Mesh(mesh.vertices[used],inv.reshape(-1,3),FRAME,{k:v[used]for k,v in mesh.attributes.items()},dict(mesh.source,source_face_ids=selected.tolist()));path=OUT/('scan'+str(index)+'-root'+str(domain['index'])+'-peer'+str(len(hits))+'.npz');local.save(path);hits.append(dict(root=domain['index'],id=pid,path=str(path),actual_triangle_aabb_candidates=len(selected),source_face_ids=selected.tolist(),full_geometry_sha256=mesh.geometry_sha,full_attributes_sha256=mesh.attributes_sha,bank_sha256=mesh.source['bank_sha256']))
 report=dict(batch=index,parts_decoded=count,candidates=hits,method='Every current preferred active mesh decoded with actual transforms; actual triangle AABB intersection is a conservative candidate filter, no metadata-only bounds exclusion.',elapsed_seconds=time.monotonic()-start);atomic_json(OUT/('scan'+str(index)+'.json'),report);print(json.dumps(dict(batch=index,decoded=count,root_peer_candidates=len(hits),elapsed_seconds=report['elapsed_seconds'])))

def contact(index):
 start=time.monotonic();plan=json.loads((OUT/'scan-plan.json').read_text());report=json.loads((OUT/('scan'+str(index)+'.json')).read_text());meshes={};results=[]
 for candidate in report['candidates']:
  peer=Mesh.load(candidate['path']);root=candidate['root'];directory=OUT/('root'+str(root));row=dict(**candidate)
  for role,file in [('old_outer','old-local-outer.npz'),('new_outer','emitted-local-outer.npz'),('new_wall','emitted-local-wall.npz')]:
   key=(root,role)
   if key not in meshes:meshes[key]=Mesh.load(directory/file)
   arrays,summary=pair_contacts(meshes[key],peer,0.);assert summary['complete_within_pair_budget'];cross=arrays['pairs'][np.isin(arrays['relation'],['crossing','coplanar-overlap'])];row[role]=dict(summary,source_peer_face_ids=sorted(set(int(candidate['source_face_ids'][j])for j in cross[:,1])))
  row['new_outer_peer_faces_absent_from_source_baseline']=sorted(set(row['new_outer']['source_peer_face_ids'])-set(row['old_outer']['source_peer_face_ids']));row['new_wall_peer_faces_absent_from_source_baseline']=sorted(set(row['new_wall']['source_peer_face_ids'])-set(row['old_outer']['source_peer_face_ids']));results.append(row)
  atomic_json(OUT/('contacts'+str(index)+'-peer'+str(len(results)-1)+'.json'),row)
 result=dict(batch=index,root_peer_candidates=len(results),results=results,elapsed_seconds=time.monotonic()-start);atomic_json(OUT/('contacts'+str(index)+'.json'),result);print(json.dumps(dict(batch=index,candidates=len(results),crossing_peers=[dict(root=r['root'],id=r['id'],old=r['old_outer']['crossing_face_pairs'],new=r['new_outer']['crossing_face_pairs'],wall=r['new_wall']['crossing_face_pairs'],new_outer_peer_faces=len(r['new_outer_peer_faces_absent_from_source_baseline']),new_wall_peer_faces=len(r['new_wall_peer_faces_absent_from_source_baseline']))for r in results if r['old_outer']['crossing_face_pairs']or r['new_outer']['crossing_face_pairs']or r['new_wall']['crossing_face_pairs']],elapsed_seconds=result['elapsed_seconds'])))

def finish():
 plan=json.loads((OUT/'scan-plan.json').read_text());scans=[json.loads((OUT/('scan'+str(i)+'.json')).read_text())for i in range(len(plan['batches']))];assert sum(r['parts_decoded']for r in scans)==plan['preferred_meshes'];rows=[]
 for index,scan in enumerate(scans):
  if scan['candidates']:
   report=json.loads((OUT/('contacts'+str(index)+'.json')).read_text());assert len(report['results'])==len(scan['candidates']);rows.extend(report['results'])
 baseline=[dict(root=r['root'],id=r['id'],old_outer_pairs=r['old_outer']['crossing_face_pairs'],new_outer_pairs=r['new_outer']['crossing_face_pairs'],new_wall_pairs=r['new_wall']['crossing_face_pairs'])for r in rows if r['old_outer']['crossing_face_pairs']]
 new=[dict(root=r['root'],id=r['id'],new_outer_peer_faces=r['new_outer_peer_faces_absent_from_source_baseline'],new_wall_peer_faces=r['new_wall_peer_faces_absent_from_source_baseline'])for r in rows if r['new_outer_peer_faces_absent_from_source_baseline']or r['new_wall_peer_faces_absent_from_source_baseline']]
 report=dict(schema='complete-current-preferred-local-root-context/v1',installed=False,source_sha256=plan['source_sha256'],resolved_sha256=plan['resolved_sha256'],preferred_meshes_actually_decoded=plan['preferred_meshes'],bounded_scan_batches=len(scans),all_eligible_current_meshes_examined=True,source_geometry_and_attributes_unchanged=True,conservative_actual_triangle_AABB_root_peer_candidates=len(rows),all_candidate_contacts_complete=True,source_baseline_crossing_peers=baseline,new_peer_face_contacts=new,no_new_outer_peer_faces=not any(r['new_outer_peer_faces_absent_from_source_baseline']for r in rows),local_context_gate_passed=not new,installation_ready=False,blocking_reason='The newly finite inner wall crosses six faces of a separate native temporal MCA branch at the first root. That interface is not declared or open, and native cerebral source crossings remain. No installation is permitted until the actual arterial port and fixed-recipient packing are resolved.',scope='Two declared changed-junction envelopes only. Source outside faces are literal. Full-course postcentral/gyrus source baseline contacts are retained separately in native-mca-controls-v3/context-and-collection-proof.json; upstream M2 and remaining native cut-mouth continuity are unresolved.',scan_elapsed_seconds_sum=sum(r['elapsed_seconds']for r in scans));atomic_json(OUT/'context-gate.json',report);print(json.dumps(report))

if __name__=='__main__':
 if len(sys.argv)<2:raise SystemExit('Choose prepare, scan BATCH, contact BATCH or finish')
 if sys.argv[1]=='prepare':prepare()
 elif sys.argv[1]=='scan':scan(int(sys.argv[2]))
 elif sys.argv[1]=='contact':contact(int(sys.argv[2]))
 elif sys.argv[1]=='finish':finish()
 else:raise SystemExit('Unknown stage')
