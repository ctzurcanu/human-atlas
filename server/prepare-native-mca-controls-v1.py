"""Prepare two qualified native control identities; do not publish geometry."""
from pathlib import Path
import copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3'
OUT=B/'proposal-v2';FRAME='model:local-ta98:assembled-metres'

def face_keys(mesh):
 t=mesh.triangles;order=np.lexsort((t[:,:,2],t[:,:,1],t[:,:,0]),axis=1)
 return np.ascontiguousarray(np.take_along_axis(t,order[:,:,None],axis=1)).view(np.dtype((np.void,t.dtype.itemsize*9))).reshape(-1)

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);a=json.loads((B/'source-atlas.json').read_bytes());before=copy.deepcopy(a);parts={p['id']:p for p in a['parts']};concepts={c['id']:c for c in a['concepts']};reader=AtlasReader();proof=[];render=[];groups=[]
 for side in ['l','r']:
  ids=['ZA:'+name+' lateral striate branches.'+side for name in ['Proximal','Distal']];meshes=[reader.read_part(a,parts[pid],frame=FRAME)for pid in ids]
  rows=[]
  for i,(pid,m) in enumerate(zip(ids,meshes)):
   p=OUT/(side+'-striate-'+str(i)+'.npz');m.save(p);render.append(dict(mesh=str(p),name=pid,color=[.88,.22,.13,1.]if i==0 else[.95,.46,.20,1.]));rows.append(dict(id=pid,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,bank=m.source['bank_path'],bank_sha256=m.source['bank_sha256'],source_outside_surface_unchanged=True))
  cid='TA98NATIVE:male:MCA:anterolateral-central:'+side;assert cid not in concepts
  limitation='Collection of the existing native proximal and distal lateral striate outside-surface networks, the two TA98 child classes of anterolateral central arteries. Every source triangle is drawn once through its original owner. Source casts are open/disconnected; this control creates no finite arterial wall, patent common junction, recipient perfusion or corrected vessel course.'
  c=dict(id=cid,name='Anterolateral central arteries ('+('left'if side=='l'else'right')+')',elements=ids,ta98Term='A12.2.07.048',ta98Kind='exact',ta98Method='complete-authority-child-class-native-source-collection',coverageLimitation=limitation)
  a['concepts'].append(c);groups.append(dict(term='A12.2.07.048',concept_id=cid,source_part_ids=ids,qualification=limitation));proof.append(dict(side=side,unique_existing_part_owners=len(set(ids)),source_geometry_not_copied_or_reconstructed=True,sources=rows))
  postcentral=concepts['ZA:Postcentral arterial branch.'+side];assert postcentral['elements']==['ZA:Postcentral arterial branch.'+side]and 'ta98Term'not in postcentral
  postcentral.update(ta98Term='A12.2.07.066',ta98Kind='exact',ta98Method='native-postcentral-name-and-mca-course-candidate',coverageLimitation='Existing native Postcentral arterial branch outside-surface source, proposed as artery of postcentral sulcus (TA98 A12.2.07.066; FMA50482). The source has two open components and retains its original MCA M3 root overlaps. The source follows the postcentral gyrus region; precise sulcal alignment, cortical penetration, finite wall, communicating lumen and perfusion remain unverified. Requires actual anatomical render review before publication.')
 a['nativeMcaControlIntegration']=dict(schema='native-mca-controls/v1',sourceCatalogueSha256=digest_file(B/'source-atlas.json'),terms=['A12.2.07.048','A12.2.07.066'],geometryUnchanged=True,sourceReferencesOnly=True,anatomicalAcceptance=False,proposedOnly=True)
 assert a['parts']==before['parts']and a['chunks']==before['chunks'];assert all(a.get(k)==before.get(k)for k in before if k!='concepts')
 proposed=OUT/'proposed-atlas.json';proposed.write_text(json.dumps(a,separators=(',',':'))+'\n')
 atomic_json(OUT/'proposal-proof.json',dict(schema='native-mca-control-proposal/v1',source_catalogue_sha256=digest_file(B/'source-atlas.json'),proposed_catalogue_sha256=digest_file(proposed),terms=['A12.2.07.048','A12.2.07.066'],new_distinct_source_meshes=0,new_triangles=0,new_geometry_bytes=0,all_source_geometry_material_visibility_and_ownership_literal=True,collections=groups,collection_geometry_checks=proof,deferred_incomplete_collections=['A12.2.07.046','A12.2.07.055','A12.2.07.061'],deferred_ambiguous_branches=['A12.2.07.057','A12.2.07.067'],postcentral_context_evidence=str(B/'context-and-collection-proof.json'),anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
 atomic_json(OUT/'striate-review-spec.json',dict(frame=FRAME,meshes=render,resolution=900,wall_seconds=40,views=[dict(name='front',direction=[0,0,1]),dict(name='back',direction=[0,0,-1]),dict(name='inferior',direction=[0,-1,0])]))
 inventory=json.loads((B/'source-bound-inventory.json').read_bytes())
 for row in inventory['terms']:
  code=row['authority']['id']
  if code=='A12.2.07.048':planned=[pid for group in groups for pid in group['source_part_ids']]
  elif code=='A12.2.07.066':planned=['ZA:Postcentral arterial branch.'+s for s in ['l','r']]
  else:continue
  row['planned_representation']['part_ids']=planned
  existing={r['part_id']for r in row['source_meshes']};row['source_meshes'].extend(dict(part_id=pid,source_catalogue=str(B/'source-atlas.json'))for pid in planned if pid not in existing)
 atomic_json(OUT/'source-bound-inventory.json',inventory)
 print(json.dumps(dict(proposed=str(proposed),sha256=digest_file(proposed),terms=2,geometry_changes=0,installed=False,elapsed_seconds=time.monotonic()-start)))

if __name__=='__main__':main()
