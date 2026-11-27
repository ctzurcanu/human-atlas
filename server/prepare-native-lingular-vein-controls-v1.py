"""Bind three existing native lingular venous sources without geometry changes."""
from pathlib import Path
import copy,json,sqlite3,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_bytes,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology
from geometry_tools.ta98_authority import corrected_terms
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
SHA='66b4a76c63043a38639d483dd744faa3789b4e468685e2746f70f738b1c1d11b'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-next-venous-controls-v1'
PAIRS=[('A12.3.02.308','Lingular vein of left lung'),('A12.3.02.309','Superior lingular vein of left lung'),('A12.3.02.310','Inferior lingular vein of left lung')]

def main():
 start=time.monotonic();assert digest_file(SOURCE)==SHA
 before=json.loads(SOURCE.read_bytes());after=copy.deepcopy(before);parts={p['id']:p for p in after['parts']};concepts={c['id']:c for c in after['concepts']}
 oldcensus=json.loads((ROOT/'.local-models/geometry-workbench/circulatory-completion/native-pulmonary-controls-v1/census-installed-v1.json').read_bytes());census={r['id']:r for r in oldcensus['terms']}
 reader=AtlasReader();evidence=[];face_sets=[]
 for term,name in PAIRS:
  pid='ZA:'+name;p=parts[pid];c=concepts[pid]
  assert not p.get('suppressed') and p.get('indexCount',0)>0 and c['elements']==[pid] and not c.get('ta98Term')
  assert p['sourceId']==name and p['system']=='venous' and p['license']=='CC BY-SA 4.0'
  assert 'Pulmonary veins'in p['groups'] and "Left superior pulmonary vein'"in p['groups']
  assert not census[term]['resolved_direct_part_ids']
  mesh=reader.read_part(after,p,frame='model:local-ta98:assembled-metres');top=topology(mesh)
  face_set={tuple(sorted(tuple(point)for point in tri))for tri in mesh.triangles}
  assert all(not face_set.intersection(prev)for prev in face_sets);face_sets.append(face_set)
  limitation=('Existing named native lingular venous outside surface, independently selectable under TA98. '
   f"Source retains {len(top['components'])} component(s) and {top['open_edges']} open edges. "
   'Identity and hierarchy control only: geometry, source visibility and recipient interfaces are unchanged. '
   'Finite wall thickness, branch-junction ownership, lumen continuity and pulmonary drainage remain unverified.')
  c.update(ta98Term=term,ta98Kind='exact',ta98Method='reviewed-native-source-lingular-venous-identity',coverageLimitation=limitation)
  evidence.append(dict(term=term,concept_id=pid,part_id=pid,source_id=name,source_groups=p['groups'],
   geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,bank_path=mesh.source['bank_path'],bank_sha256=mesh.source['bank_sha256'],
   existing_topology=top,qualification=limitation,literal_faces_shared_with_other_batch_parts=False))
 after['nativeLingularVenousControlIntegration']=dict(schema='native-lingular-venous-source-controls/v1',sourceCatalogueSha256=SHA,
  conceptIds=[r['concept_id']for r in evidence],partIds=[r['part_id']for r in evidence],terms=[r['term']for r in evidence],
  geometryUnchanged=True,sourceVisibilityUnchanged=True,anatomyAccepted=False,finiteWallOrPatencyClaim=False)
 assert before['parts']==after['parts'] and before['chunks']==after['chunks'] and before['materials']==after['materials'] and before['defaultHidden']==after['defaultHidden']
 atomic_bytes(OUT/'source-atlas.json',SOURCE.read_bytes());atomic_bytes(OUT/'proposed-atlas.json',(json.dumps(after,separators=(',',':'))+'\n').encode())
 authority=ROOT/'.local-models/ta98-authority.sqlite';db=sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
 terms=corrected_terms({r['id']:dict(r)for r in db.execute('select * from ta98')});closure={'A12.3.02.001'}
 while True:
  enlarged=closure|{k for k,v in terms.items()if v.get('parent_id')in closure}
  if enlarged==closure:break
  closure=enlarged
 assert len(closure)==44,len(closure)
 selected={r['term']:r['part_id']for r in evidence};inventory=[]
 for code in sorted(closure):
  old=census[code];ids=sorted(set(old['active_exact_source_part_ids']+old['suppressed_exact_source_part_ids']+([selected[code]]if code in selected else [])))
  inventory.append(dict(authority=terms[code],representation_class='whole-organ',distinct_representation_verified=False,representation_provisional=True,
   missing_distinct_representation=not bool(old['resolved_direct_part_ids']),source_meshes=[dict(part_id=pid,source_catalogue=str(OUT/'source-atlas.json'))for pid in ids],
   planned_representation=dict(part_ids=[selected[code]]if code in selected else [],scope='Native outside-surface identity controls only; finite wall, communicating lumen and recipient fit remain unverified.')))
 inv=dict(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',root_term='A12.3.02.001',term_count=len(closure),applicable_term_count=len(closure),terms=inventory,
  source_catalogue=str(OUT/'source-atlas.json'),source_catalogue_sha256=SHA,authority_path=str(authority),authority_sha256=digest_file(authority),authority_corrections_path=str(ROOT/'geometry_tools/ta98_authority.py'),
  authority_corrections_sha256=digest_file(ROOT/'geometry_tools/ta98_authority.py'),checked_in_report=str(ROOT/'reports/mesh-audit/REVIEW.md'))
 atomic_json(OUT/'source-bound-inventory.json',inv)
 proof=dict(schema='native-lingular-venous-identity-control-proof/v1',source_catalogue_sha256=SHA,proposed_catalogue_sha256=digest_file(OUT/'proposed-atlas.json'),
  existing_part_records_literal=True,existing_banks_literal=True,materials_literal=True,visibility_literal=True,distinct_new_direct_reference_terms=3,
  new_mesh_records=0,new_triangles=0,new_geometry_bytes=0,anatomical_acceptance=False,installed=False,rows=evidence,elapsed_seconds=time.monotonic()-start)
 atomic_json(OUT/'identity-proof.json',proof);atomic_json(OUT/'concept-delta.json',dict(source_catalogue_sha256=SHA,concepts=[concepts[r['concept_id']]for r in evidence],integration_key='nativeLingularVenousControlIntegration',integration=after['nativeLingularVenousControlIntegration']))
 assert digest_file(SOURCE)==SHA
 print(json.dumps(dict(output=str(OUT),proposed_sha256=proof['proposed_catalogue_sha256'],references=3,authority_terms=len(closure),elapsed_seconds=proof['elapsed_seconds'])))

if __name__=='__main__':main()
