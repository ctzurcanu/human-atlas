"""Expose finite native chordae without moving or redrawing their source tissue."""
from pathlib import Path
import argparse,copy,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--bank',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args()
h=ROOT/'.local-models/geometry-workbench/heart';bank=h/args.bank;proof=json.loads((bank/'construction.json').read_text());assert proof['passed']and digest_file(bank/'cells.bin')==proof['bank_sha256'];sex=args.sex
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json';a=json.loads(source.read_text());before=copy.deepcopy(a);parts={p['id']:p for p in a['parts']};concepts={c['id']:c for c in a['concepts']};marker=a['finiteCardiacIntegration'];assert not marker.get('chordalPartIds');chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=proof['bytes']));reader=AtlasReader();new_ids=[];evidence=[];groups=marker['termGroups'];whole=concepts[marker['wholeConceptId']];parent_id=f'TA98HEART:{sex}:chordae'
papillary={1:'ZA:Inferior papillary muscle of right ventricle',5:'ZA:Inferior papillary muscle of left ventricle',8:'ZA:Septal papillary muscle of right ventricle',9:'TA98CARDIAC:FJ2418',10:'ZA:Anterior papillary muscle of right ventricle'}
def recipient(pid):
 if sex=='male':return pid
 matches=[p['id']for p in a['parts']if p.get('cardiacTemplateSourcePartId')==pid];assert len(matches)==1,(pid,matches);return matches[0]
for row in proof['parts']:
 old_id=recipient(row['source_part_id']);original=copy.deepcopy(parts[old_id]);chord_id=f'TA98HEART:{sex}:chordae:source-{row["source"]}';pm_id=recipient(papillary[row['source']]);mapping={row['source_part_id']:old_id,row['part_ids']['chordae']:chord_id};pair=[]
 for role,pid in [('leaflet',old_id),('chordae',chord_id)]:
  p=copy.deepcopy(original);layout=copy.deepcopy(row['parts'][role]);
  for key in ['positions','normals','indices','normalSigns','uvs','vertexCount','indexCount','triangles','tissueDrawGroups']:p.pop(key,None)
  p.update({k:v for k,v in layout.items()if k not in ['id','bounds','geometryOrigin','geometryScale','geometryInstance','sourceOffset']});p.update(id=pid,chunk=chunk,pairedTissue=True,sectionSolid=True,tissueWallAssembly=f'heart:{sex}',tissueWallLayer='chordae'if role=='chordae'else'valve',defaultHidden=False)
  for g in p['tissueDrawGroups']:
   if g.get('neighbor'):g['neighbor']=mapping[g['neighbor']]
  p['cardiacChordalPartition']=dict(schema='native-chordal-partition/v1',role=role,leafletPartId=old_id,chordalPartId=chord_id,papillaryPartId=pm_id,sharedLeafletInterface=True,papillaryInterfaceVerified=False,bankSha256=proof['bank_sha256'],sourcePartId=row['source_part_id'],sourceGeometrySha256=row['source_geometry_sha256'])
  p['coverageLimitation']='Native leaflet and branching chordal source tissue, separated at measured strand necks. The ownership boundary is illustrative; papillary ends remain adjacent and are not a verified shared tissue join.'
  if role=='chordae':
   side='left'if row['source']in[5,9]else'right';leaflet=('Septal'if row['source']==8 else'Anterior'if row['source']in[9,10]else'Posterior').lower();p.update(name=f'Chordae tendineae · {leaflet} leaflet ({side})',conceptId=chord_id,ta98Term='A12.1.00.023');p.pop('cardiacIdentityEvidence',None);p.pop('cardiacTemplateSourcePartId',None)
   c=dict(id=chord_id,name=p['name'],elements=[pid],ta98Term='A12.1.00.023',ta98Kind='parent',partOf=parent_id,coverageLimitation=p['coverageLimitation']);a['concepts'].append(c);groups.append(dict(key=f'chordae-source-{row["source"]}',name=p['name'],partIds=[pid],term='A12.1.00.023',conceptId=chord_id,parent='chordae',kind='parent'));a['parts'].append(p);new_ids.append(pid)
  else:parts[old_id].clear();parts[old_id].update(p);p=parts[old_id]
  mesh=reader.read_part(a,p,frame='local-ta98'if sex=='male'else'local-female-ta98');p['bounds']=mesh.bounds.tolist();pair.append(p)
  material=a['materials'][p['material']];material['pairedNormalSigns']=True
 evidence.append(dict(source=row['source'],leaflet=old_id,chordae=chord_id,papillary=pm_id,geometry_origin=original.get('geometryOrigin'),geometry_instance=original.get('geometryInstance'),world_bounds={p['id']:p['bounds']for p in pair}))
marker['chordalPartIds']=new_ids;marker['chordalBankSha256']=proof['bank_sha256'];marker['chordalPartitionEvidence']=str(bank/'construction.json');whole['elements']+=new_ids;next(g for g in groups if g['key']=='heart')['partIds']+=new_ids
a['concepts'].append(dict(id=parent_id,name='Chordae tendineae',elements=new_ids,ta98Term='A12.1.00.023',ta98Kind='exact',partOf=marker['wholeConceptId'],coverageLimitation='Collective control of the five source chordal fans. Papillary insertions and the complete finer chordal anatomy remain unverified.'));groups.append(dict(key='chordae',name='Chordae tendineae',partIds=new_ids,term='A12.1.00.023',conceptId=parent_id,parent='heart',kind='exact'))
# The recovered source is explicitly a left-ventricular papillary head.
for g in groups:
 if 'TA98CARDIAC:FJ2418'in[p.get('cardiacTemplateSourcePartId',p['id'])for p in a['parts']if p['id']in g['partIds']]and g['key'].startswith('papillary-'):g['parent']='left-ventricle'
remaining='Finite chamber walls, septal regions, endocardium, epicardium and native valves with five separately selectable chordal fans. Septal and strand-neck ownership boundaries and lining thickness are illustrative. Papillary insertions, complete fine valve anatomy, conduction tissues, pericardial detail and distal recipient vessel joins remain unfinished.';marker['remainingLimits']=remaining;whole['coverageLimitation']=remaining
out=h/args.stage;out.mkdir(exist_ok=False);atomic_json(out/'atlas.json',a);changed=[p['id']for p in a['parts']if p!=next((q for q in before['parts']if q['id']==p['id']),None)];atomic_json(out/'candidate.json',dict(sex=sex,source_catalogue=str(source),source_catalogue_sha256=digest_file(source),catalogue_sha256=digest_file(out/'atlas.json'),bank=str(bank/'cells.bin'),bank_sha256=proof['bank_sha256'],changed_parts=changed,physical_parts=len(whole['elements']),direct_terms=len({g['term']for g in groups if g['kind']=='exact'}),pairs=evidence,remaining_limits=remaining,installed=False));print(sex,len(changed),'changed parts;',len(whole['elements']),'heart tissues; one shared',proof['bytes'],'byte bank; no new textures',flush=True)
