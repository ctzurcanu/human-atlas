"""Bound a native MCA control proposal to authority, source geometry and context."""
from pathlib import Path
import json, sqlite3, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.ta98_authority import corrected_terms
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

B=(ROOT/sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v1'
FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();B.mkdir(exist_ok=False)
 master=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';raw=master.read_bytes();(B/'source-atlas.json').write_bytes(raw);atlas=json.loads(raw);parts={p['id']:p for p in atlas['parts']}
 resolved_path=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-blood-controls-v2/resolved-installed.json';resolved=json.loads(resolved_path.read_bytes())
 authority_path=ROOT/'.local-models/ta98-authority.sqlite';db=sqlite3.connect(authority_path);db.row_factory=sqlite3.Row;terms=corrected_terms({r['id']:dict(r)for r in db.execute('select * from ta98')});closure=[];pending=['A12.2.07.046']
 while pending:
  code=pending.pop(0)
  if code in closure:continue
  closure.append(code);pending.extend(k for k,v in terms.items()if v['parent_id']==code)
 assert len(closure)==24
 if resolved['sourceSha256']!=digest_file(B/'source-atlas.json'):
  previous_path=resolved_path.parent/'proposed-atlas.json';assert digest_file(previous_path)==resolved['sourceSha256'];previous=json.loads(previous_path.read_bytes());assert previous['chunks']==atlas['chunks']
  old_parts={p['id']:p for p in previous['parts']};assert old_parts.keys()==parts.keys();changed={pid for pid in parts if old_parts[pid]!=parts[pid]};allowed={'ZA:'+name+' pontine branches of basilar artery.'+side for name in ['Medial','Lateral']for side in ['l','r']};assert changed<=allowed
  assert all({k:v for k,v in old_parts[pid].items()if k!='ta98Term'}=={k:v for k,v in parts[pid].items()if k!='ta98Term'}for pid in parts)
  old_concepts={c['id']:c for c in previous['concepts']};new_concepts={c['id']:c for c in atlas['concepts']};scope={r['conceptId']for r in resolved['conceptBindings']if r['term']in closure};scope.update(c['id']for c in previous['concepts']if c.get('ta98Term')in closure)
  assert all(old_concepts[cid]==new_concepts[cid]for cid in scope)
  old_hierarchy=resolved_path.parent/'future-relations-workspace-v1/public/assets/ta98.json';old_nodes={n['id']:n for n in json.loads(old_hierarchy.read_bytes())['nodes']};new_nodes={n['id']:n for n in json.loads((ROOT/'public/assets/ta98.json').read_bytes())['nodes']};assert all(old_nodes[code]==new_nodes[code]for code in closure)
  atomic_json(B/'scoped-resolver-equivalence.json',dict(schema='literal-mca-snapshot-equivalence/v1',resolved_source_sha256=resolved['sourceSha256'],current_source_sha256=digest_file(B/'source-atlas.json'),all_geometry_records_literal=True,unrelated_identity_tag_changes=sorted(changed),mca_authority_nodes_literal=True,mca_concepts_literal=True,checked_concepts=len(scope),scope_terms=closure,global_census_reused=False))
 direct={n['id']:n['directPartIds']for n in resolved['resolved']};binding={code:[b for b in resolved['conceptBindings']if b['term']==code]for code in closure}
 inventory=dict(schema='ta98-geometry-inventory/v1',root_term='A12.2.07.046',term_count=24,applicable_term_count=24,model='local-ta98',model_sex='male',source_catalogue=str(B/'source-atlas.json'),source_catalogue_sha256=digest_file(B/'source-atlas.json'),authority_path=str(authority_path),authority_sha256=digest_file(authority_path),authority_corrections_path=str(ROOT/'geometry_tools/ta98_authority.py'),authority_corrections_sha256=digest_file(ROOT/'geometry_tools/ta98_authority.py'),terms=[])
 for code in closure:
  ids=sorted({pid for b in binding[code]for pid in b['partIds']if pid in parts});inventory['terms'].append(dict(authority=terms[code],representation_class='whole-organ',representation_provisional=True,distinct_representation_verified=False,missing_distinct_representation=not bool(direct.get(code)),source_meshes=[dict(part_id=pid,source_catalogue=str(B/'source-atlas.json'))for pid in ids],resolved_direct_part_ids=direct.get(code,[]),planned_representation=dict(part_ids=[],scope='Native outside-source identity/collection only; no finite wall or communicating lumen acceptance.')))
 atomic_json(B/'source-bound-inventory.json',inventory)
 reader=AtlasReader();cache={}
 def mesh(pid):
  if pid not in cache:cache[pid]=reader.read_part(atlas,parts[pid],frame=FRAME)
  return cache[pid]
 context=[];render_entries=[]
 for side in ['l','r']:
  pid='ZA:Postcentral arterial branch.'+side;m=mesh(pid);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];loops=boundary_loops(f);mouths=np.array([v[loop].mean(0)for loop in loops]);assert len(mouths)==4
  path=B/(side+'-postcentral.npz');m.save(path);render_entries.append(dict(mesh=str(path),name=pid,color=[.85,.16,.12,1.]))
  parent_ids=['ZA:Middle cerebral artery (M3 segment).l'if side=='l'else'ZA:Middle cerebral artery (M3-segment).r','ZA:Insular branches of middle cerebral artery (M2-segment).l'if side=='l'else'ZA:Insular branches of middle cerebral artery (M2).r','ZA:Artery of precentral sulcus.'+side,'ZA:Posterior parietal artery.'+side,'ZA:Anterior cerebral artery.'+side]
  peers=[]
  for peer in parent_ids:
   q=mesh(peer);near,_=SurfaceIndex(q).closest_points(mouths,certify=True);hits,proof=pair_contacts(m,q,0.);assert proof['complete_within_pair_budget'];shared=set(map(tuple,m.vertices)).intersection(map(tuple,q.vertices));peers.append(dict(id=peer,mouth_nearest_actual_surface_m=near['distance_m'].tolist(),literal_shared_coordinates=len(shared),actual_triangle_contacts=proof))
   if 'Middle cerebral artery (M3' in peer:
    p=B/(side+'-m3.npz');q.save(p);render_entries.append(dict(mesh=str(p),name=peer,color=[.58,.12,.10,1.]))
  recipients=[]
  for peer in ['ZA:Postcentral sulcus.'+side+':Parietal lobe','ZA:Postcentral sulcus.'+side+':Brain-Inner','ZA:Postcentral gyrus.'+side+':Parietal lobe','ZA:Central sulcus.'+side+':Interlobar sulci']:
   q=mesh(peer);near,_=SurfaceIndex(q).closest_points(m.vertices,certify=True);hits,proof=pair_contacts(m,q,0.);assert proof['complete_within_pair_budget'];recipients.append(dict(id=peer,minimum_actual_surface_to_source_vertex_m=float(near['distance_m'].min()),source_vertex_median_nearest_surface_m=float(np.median(near['distance_m'])),actual_triangle_contacts=proof))
   p=B/(side+'-'+peer.split(':')[-1].replace(' ','-')+'-'+str(len(recipients))+'.npz');q.save(p);render_entries.append(dict(mesh=str(p),name=peer,color=[.52,.48,.43,1.]))
  context.append(dict(id=pid,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,source_bank=m.source['bank_path'],source_bank_sha256=m.source['bank_sha256'],topology=topology(m),mouth_centers=mouths.tolist(),parent_candidates=peers,recipient_context=recipients))
 assembly_codes=['A12.2.07.046','A12.2.07.048','A12.2.07.055','A12.2.07.061'];assembly_rows=[]
 for code in assembly_codes:
  children=[k for k,v in terms.items()if v['parent_id']==code];assembly_rows.append(dict(term=code,name=terms[code]['name_en'],authority_direct_children=children,children_with_direct_control=[k for k in children if direct.get(k)],children_missing_direct_control=[k for k in children if not direct.get(k)],native_source_parts={k:[pid for pid in direct.get(k,[])if pid.startswith(('ZA:','TA98PART:ZA:'))]for k in children},complete_source_collection_candidate=(code=='A12.2.07.048'),qualification='Only .048 has both complete named child classes present; .046/.055/.061 remain incomplete or ambiguous and must not be promoted as exact complete collections.'))
 atomic_json(B/'context-and-collection-proof.json',dict(schema='native-mca-controls-context/v1',source_catalogue_sha256=digest_file(B/'source-atlas.json'),frame=FRAME,geometry_modified=False,skeleton_modified=False,anatomical_acceptance=False,postcentral=context,collection_candidates=assembly_rows,elapsed_seconds=time.monotonic()-start))
 atomic_json(B/'review-spec.json',dict(frame=FRAME,meshes=render_entries,resolution=900,wall_seconds=40,views=[dict(name='front',direction=[0,0,1]),dict(name='back',direction=[0,0,-1]),dict(name='left',direction=[1,0,0]),dict(name='right',direction=[-1,0,0])]))
 print(json.dumps(dict(source_sha256=digest_file(B/'source-atlas.json'),terms=24,proven_collection='A12.2.07.048',context=str(B/'context-and-collection-proof.json'),elapsed_seconds=time.monotonic()-start)))

if __name__=='__main__':main()
