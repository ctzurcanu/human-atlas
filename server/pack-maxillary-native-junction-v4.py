"""Pack an already-passed six-source maxillary exterior for local review only.

No geometry fitting, live atlas mutation, database job or installation occurs.
"""
import argparse, copy, hashlib, importlib.util, json, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops
W=ROOT/'.local-models/geometry-workbench/circulatory-completion/maxillary-native-v1'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,default=lambda a:a.tolist()if hasattr(a,'tolist')else str(a))+'\n')
def run(candidate,out,storage_failure_report=None,prepared_storage_packet=None):
 t=time.monotonic();proofpath=candidate/'proof.json';p=json.loads(proofpath.read_text());assert p['schema']=='maxillary-local-root-completion-candidate/v4'and p['complete']and p['numerical_gate_pass'],'Only completed passing v4 may be packed'
 assert H(CAT)==p['source_master_sha256']and H(candidate/'cohort.npz')==p['candidate_sha256']
 assert not p['self_contact']['defect_pairs']and p['self_contact']['complete_within_pair_budget']and not p['unexpected_geometric_boundary_positions_m']and not p['nonforward_normal_faces']
 assert p['all_classifier_caps_removed']and p['all_output_faces_have_one_owner']and p['original_control_displacement_m']==0
 assert len(p['external_ports'])==22 and all(all(port[k]for k in ['positions_literal','normal_words_literal','uv_literal','all_on_boundary'])for port in p['external_ports'])
 planpath=candidate.parent/'plan.json';plan=json.loads(planpath.read_text());assert plan['source_master_sha256']==H(CAT)
 invpath=W/'complete-authority-inventory.json';inv=json.loads(invpath.read_text());assert inv['term_count']==len(inv['terms'])==38 and inv['source_catalogue_sha256']==H(CAT)
 rows={r['source']:r for r in plan['sources']};selected=plan['selected_whole_source_indices'];assert selected==[0,1,8,9,10,11];ids=[rows[i]['id']for i in selected]
 source=json.loads(CAT.read_text());records={x['id']:x for x in source['parts']};reader=AtlasReader();mesh=Mesh.load(candidate/'cohort.npz');a=np.load(candidate/'ancestry.npz',allow_pickle=False);owner=a['source_indices'];storage_packet_proof=None
 if prepared_storage_packet is not None:
  storage_packet_proof=json.loads((prepared_storage_packet/'proof.json').read_text());assert storage_packet_proof['complete']and storage_packet_proof['passed']
  if storage_packet_proof.get('shader_only'):assert storage_packet_proof['all_existing_position_normal_UV_rows_literal_from_v6']and storage_packet_proof['all_geometric_triangles_literal_from_v6']and storage_packet_proof['corrected_fragment_faces']==12
  else:assert storage_packet_proof['all_positions_normal_UV_arrays_literal_from_v5']
  assert storage_packet_proof['source_master_sha256']==H(CAT)and storage_packet_proof['producer_v4_sha256']==H(candidate/'cohort.npz')and storage_packet_proof['cohort_sha256']==H(prepared_storage_packet/'cohort.npz')and storage_packet_proof['maximum_surface_Hausdorff_bound_m']<=1e-8
  mesh=Mesh.load(prepared_storage_packet/'cohort.npz');owner=np.load(prepared_storage_packet/'ownership.npz',allow_pickle=False)['source_indices']
 assert len(owner)==len(mesh.faces)and set(owner)==set(selected)
 for i in selected:
  assert H(rows[i]['path'])==rows[i]['sha256'];frozen=Mesh.load(rows[i]['path']);actual=reader.read_part(source,records[rows[i]['id']],frame=frozen.frame);assert actual.geometry_sha==frozen.geometry_sha and actual.attributes_sha==frozen.attributes_sha
 assert not out.exists();out.mkdir(parents=True);write(out/'construction-state.json',dict(complete=False,passed=False,installed=False))
 vertices=mesh.vertices.copy();normal=mesh.attributes['normals'].copy();uv=mesh.attributes['uvs'].copy();faces=mesh.faces;prefix=a['source_position_prefix'];assert np.array_equal(vertices[:len(prefix)],prefix)
 # One coordinate table; only new intersection charts may receive a shared
 # normal. Original source word controls cannot silently change in packing.
 used=np.unique(faces);unique,inverse=np.unique(vertices,axis=0,return_inverse=True);bits=np.zeros(len(unique),np.uint8)
 for role,index in enumerate(selected):np.bitwise_or.at(bits,inverse[np.unique(faces[owner==index])],np.uint8(1<<role))
 shared_groups=np.flatnonzero((bits&(bits-1))!=0);canonicalized=0
 for group in ([]if storage_packet_proof and storage_packet_proof.get('shader_only')else shared_groups):
  indices=used[inverse[used]==group];words=np.unique(normal[indices],axis=0)
  if len(words)==1:continue
  word=words.astype(float).sum(0);assert np.linalg.norm(word)>0;word=np.rint(word/np.linalg.norm(word)*32767).astype(np.int16)
  original=indices[indices<len(prefix)];assert not len(original)or np.all(normal[original]==word),'Shared-normal canonicalization would change an original source word'
  normal[indices]=word;canonicalized+=1
 cleanup_proof=None
 if storage_failure_report is not None:
  failure=json.loads(storage_failure_report.read_text());assert not failure['passed']
  spec=importlib.util.spec_from_file_location('maxillary_storage_cleanup',ROOT/'server/maxillary-storage-cleanup-v5.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  vertices,faces,normal,uv,owner,cleanup_proof=module.cleanup(vertices,faces,normal,uv,owner,a,failure,mesh.frame)
  write(out/'storage-cleanup-proof.json',cleanup_proof)
 origin=np.array([0.,1.671875,0.]);stored=(vertices-origin).astype('<f4');actual=stored.astype(float)+origin
 assert np.array_equal(actual[:len(prefix)],prefix),'Original source controls lost Float32 identity'
 compiled=Mesh(actual,faces,mesh.frame,dict(normals=normal,uvs=uv));compiled.save(out/'compiled-cohort.npz')
 check=audit(compiled,require_closed=False,expected_components=2,max_pairs=600000,max_tests=3000000);write(out/'compiled-cohort-audit.json',check)
 if cleanup_proof is not None:cleanup_proof['actual_Float32_gate_passed']=check['passed'];cleanup_proof['actual_Float32_failures']=check['failures'];write(out/'storage-cleanup-proof.json',cleanup_proof)
 assert check['passed'],check['failures']
 loops=boundary_loops(compiled.faces)if len(np.unique(vertices,axis=0))==len(vertices)else None
 # UV/normal chart indices may repeat a geometric point. Classify actual
 # boundary cycles after exact-coordinate welding, never per-chart topology.
 gv,gi=np.unique(actual,axis=0,return_inverse=True);gf=gi[faces];loops=boundary_loops(gf);assert len(loops)==22 and all(len(x)==12 for x in loops)
 fields=np.load(candidate.parent/'source-course-arrays.npz',allow_pickle=False);expected=[]
 for port in p['external_ports']:
  descriptor=next(x for x in plan['fixed_terminal_ports']if all(x[k]==port[k]for k in ['source','component','end']));si=port['source'];frozen=Mesh.load(rows[si]['path']);ring=fields[descriptor['ring_array']][descriptor['ring_station']];points=frozen.vertices[ring];expected.append(frozenset(map(tuple,points)))
  for vi,point in zip(ring,points):
   matches=used[np.all(actual[used]==point,axis=1)];assert len(matches)
   assert all(np.array_equal(normal[j],frozen.attributes['normals'][vi])and np.array_equal(uv[j],frozen.attributes['uvs'][vi])for j in matches)
 assert len(set(expected))==22 and set(expected)=={frozenset(map(tuple,gv[loop]))for loop in loops}
 np.savez_compressed(out/'compiled-ownership.npz',vertices=actual,faces=faces,source_indices=owner,normal_words=normal,uvs=uv,shared_coordinate_ids=shared_groups)
 bank=bytearray()
 def append(x):
  x=np.ascontiguousarray(x)
  while len(bank)%x.dtype.itemsize:bank.append(0)
  offset=len(bank);bank.extend(x.tobytes());return offset
 posoff,noff,uoff=append(stored),append(normal.astype('<i2')),append(uv.astype('<f4'));after=copy.deepcopy(source);chunk=len(source['chunks'])
 limitation='Macroscopic native arterial exterior junction revision. Source c1 strips 0–3 are trimmed and its remaining provisional unnamed alveolar-side course joins the infraorbital artery through a constructed conical patch. No A12.2.05.079 or new child identity is assigned; no finite vessel wall or complete perfusion claim.'
 for i,pid in zip(selected,ids):
  record=next(x for x in after['parts']if x['id']==pid);ff=faces[owner==i].astype('<u4')
  for key in ['sourceOffset','geometryScale','normalSigns','geometryInstance']:record.pop(key,None)
  record.update(chunk=chunk,positions=posoff,normals=noff,uvs=uoff,indices=append(ff),geometryOrigin=origin.tolist(),vertexCount=len(actual),indexCount=int(ff.size),bounds=[actual[ff].min((0,1)).tolist(),actual[ff].max((0,1)).tolist()],sourceOuterSurfaceOnly=True,coverageLimitation=limitation,nativeMaxillaryJunctionRevision=dict(schema='native-local-junction-v4',sourceCatalogueSha256=H(CAT),sourceComponentOwner=i,finiteWall=False,newTA98TermCountClaim=0))
 bankpath=out/'maxillary-junction.bin';bankpath.write_bytes(bank);after['chunks'].append(dict(url='/local-models/'+str(bankpath.relative_to(ROOT/'.local-models')),bytes=len(bank)));after['triangles']=sum(x['indexCount']//3 for x in after['parts'])
 assert len(after['parts'])==len(source['parts'])and after['concepts']==source['concepts']and after['materials']==source['materials']and after['ta98TermBindings']==source['ta98TermBindings']
 assert all(x==records[x['id']]for x in after['parts']if x['id']not in ids)
 write(out/'proposed-atlas.json',after);(out/'registry-authority-inventory.json').write_bytes(invpath.read_bytes());write(out/'assembly-plan.json',dict(schema='native-maxillary-six-record-junction-replacement/v4',whole_organ='Bilateral native maxillary arterial exteriors and their existing inferior alveolar and infraorbital source courses.',source_master_sha256=H(CAT),changed_parts=ids,adjoining_interfaces=['Two side-specific exterior assemblies retain22literal exterior source mouths. Internal source junction classifiers are removed.','Only the declared local junction patch charts and c1 strips0–3 are changed; other source controls, originalnormalwords and UVs are retained.','The c1 alveolar-side branch is provisional and uncoded; no anterior superior alveolar .079 identity or extraTA98coverage is claimed.','Skeleton, brain, nerves, basilar artery, female catalogue and all other originalsource records stay fixed; actual packed context review remains required.'],contract=p['new_patch_contract'],finite_wall=False,new_TA98_term_count_claim=0,installed=False))
 result=dict(schema='native-maxillary-six-record-packed-candidate/v4',complete=True,passed=True,source_master_sha256=H(CAT),proposed_catalogue_sha256=H(out/'proposed-atlas.json'),candidate_proof_sha256=H(proofpath),ancestry_sha256=H(candidate/'ancestry.npz'),producer_cohort_sha256=H(candidate/'cohort.npz'),bank_sha256=H(bankpath),bank_bytes=len(bank),changed_parts=ids,unchanged_existing_records=len(source['parts'])-6,new_parts=0,new_concepts=0,new_TA98_term_count_claim=0,shared_coordinate_groups=len(shared_groups),new_shared_normal_groups_canonicalized=canonicalized,original_source_prefix_positions_literal=True,source_outer_surface_only=True,finite_wall=False,original_exterior_mouths=22,actual_context_passed=False,root_render_reviewed=False,installed=False,elapsed_seconds=time.monotonic()-t,storage_packet_proof_sha256=H(prepared_storage_packet/'proof.json')if prepared_storage_packet else None)
 write(out/'packing-proof.json',result);write(out/'construction-state.json',result);print(json.dumps(result))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--storage-failure-report',type=Path);p.add_argument('--prepared-storage-packet',type=Path);a=p.parse_args();run(a.candidate.resolve(),a.output.resolve(),a.storage_failure_report.resolve()if a.storage_failure_report else None,a.prepared_storage_packet.resolve()if a.prepared_storage_packet else None)
