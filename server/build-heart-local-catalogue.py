"""Prepare explicit heart controls from the shared wall and licensed AV bank."""
from pathlib import Path
import argparse,copy,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--wall',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);src=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(src.read_text());before=copy.deepcopy(a);parts={p['id']:p for p in a['parts']};concepts={c['id']:c for c in a['concepts']};wall=json.loads((b/args.wall/'emission.json').read_text());av=json.loads((b/'av-bank-v100/tissue.json').read_text());assert av['passed'];whole='TA98HEART:male:whole';prefix='TA98HEART:male:';native=['ZA:Right atrium','ZA:Right ventricle','ZA:Left atrium','ZA:Left ventricle'];retained=[p['id']for p in a['parts']if p.get('system')=='cardiac'and not p.get('suppressed')and p['id']not in set(a.get('defaultHidden',[]))and p['id']not in native];wallids=list(wall['part_ids'].values());newcusps=list(av['part_ids'].values());wholeids=wallids+retained+newcusps;groups=[];terms={};limits='Finite chamber walls and septal regions with separate endocardium and epicardium. Septal boundaries and relative lining thickness are illustrative. Fine valve anatomy, complete chordae, conduction tissues and pericardial detail remain unfinished.'
def concept(cid,name,ids,term,kind='exact',parent=whole):
 c=concepts.get(cid)
 if c is None:c=dict(id=cid);a['concepts'].append(c);concepts[cid]=c
 c.update(name=name,elements=list(dict.fromkeys(ids)),ta98Term=term,ta98Kind=kind,coverageLimitation=limits)
 if cid!=whole:c['partOf']=parent
 return c
def group(key,name,ids,term,parent='heart',kind='exact',cid=None):
 cid=cid or prefix+key;concept(cid,name,ids,term,kind,parent=next((g['conceptId']for g in groups if g['key']==parent),whole));row=dict(key=key,name=name,partIds=list(dict.fromkeys(ids)),term=term,conceptId=cid,parent=parent,kind=kind);groups.append(row);terms[term]=row;return row
layerterms={'myocardium':'A12.1.06.001','endocardium':'A12.1.05.001','epicardium':'A12.1.08.007','fibrous':'A12.1.00.015'};regionterms={'right-atrium':'A12.1.01.001','right-ventricle':'A12.1.02.001','left-atrium':'A12.1.03.001','left-ventricle':'A12.1.04.001','interatrial-septum':'A12.1.00.019','interventricular-muscular-septum':'A12.1.00.014','atrioventricular-septum':'A12.1.00.016','interventricular-membranous-septum':'A12.1.00.015'};colors={'myocardium':[155,69,63],'endocardium':[202,151,143],'epicardium':[188,108,101],'fibrous':[204,178,153]};reference=parts['ZA:Left ventricle']
inventory=json.loads((b/'source-v1-male/inventory.json').read_text());source_identity={}
for row in inventory['terms']:
 for match in row['source_meshes']:
  pid=match['part_id']
  if pid in retained and (pid not in source_identity or match['kind']=='exact'):source_identity[pid]=(row['authority']['id'],match['kind'])
for pid,(term,kind)in source_identity.items():
 c=concepts[parts[pid]['conceptId']];c.update(ta98Term=term,ta98Kind=kind)
for marker in [wall,av]:
 bank=Path(marker['bank']);assert digest_file(bank)==marker['bank_sha256'];ci=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size))
 for name,packed in marker['parts'].items():
  pid=marker['part_ids'][name];kind=name.split(':')[0];iswall=marker is wall;region=name.split(':')[1]if ':'in name else None;term=(regionterms[region]if region else layerterms[kind])if iswall else('A12.1.02.004'if name.startswith('tricuspid')else'A12.1.04.004');label=region.replace('-',' ').capitalize()if region else name.capitalize();color=colors[kind]if iswall else[212,174,155];appearance=dict(color=[v/255 for v in color],pattern='oblique'if kind=='myocardium'else'fibrous'if kind in ['epicardium','fibrous']else'solid',scale=.0007,relief=.000003 if kind!='myocardium'else.000006);material='HEART:'+kind if iswall else'HEART:valve';a['materials'][material]=dict(color=color,roughness=.63,metalness=0,pairedNormalSigns=True,tissueAppearance=appearance,geometryOrigin=marker.get('origin',[0,1.375,.03125]));donor=reference if iswall else parts['FJ2421'if name.startswith('tricuspid')else'FJ2420'];p={**packed,'chunk':ci,'conceptId':pid,'name':label if iswall else('Anterior leaflet of right atrioventricular valve'if name.startswith('tricuspid')else'Anterior leaflet of left atrioventricular valve'),'sourceId':'Heart:finite-wall'if iswall else'Heart:registered-AV-cusp','system':'cardiac','groups':['5: Cardiovascular system','Heart','Trunk'],'ta98Term':term,'pairedTissue':True,'sectionSolid':True,'material':material,'tissueWallAssembly':'heart:male','tissueWallLayer':('fibrous-septum'if kind=='fibrous'else kind)if iswall else'valve','tissueAppearance':appearance,'defaultHidden':False,'license':donor['license'],'provenance':donor['provenance'],'coverageLimitation':limits};a['parts'].append(p);parts[pid]=p;concept(pid,p['name'],[pid],term)
replacement={old:[prefix+'myocardium:'+region]for old,region in zip(native,['right-atrium','right-ventricle','left-atrium','left-ventricle'])}
for old,ids in replacement.items():parts[old].update(suppressed=True,defaultHidden=True,replacementParts=ids)
for c in a['concepts']:
 c['elements']=list(dict.fromkeys(pid for old in c.get('elements',[])for pid in replacement.get(old,[old])))
concept(whole,'Heart',wholeids,'A12.1.00.001');groups.append(dict(key='heart',name='Heart',partIds=wholeids,term='A12.1.00.001',conceptId=whole,parent=None,kind='exact'))
for old,region in zip(native,['right-atrium','right-ventricle','left-atrium','left-ventricle']):
 ids=[prefix+'myocardium:'+region];pap=[pid for pid in retained if 'papillary'in parts[pid]['name'].lower()and region.replace('-',' ')in parts[pid]['name'].lower()];group(region,region.replace('-',' ').capitalize(),ids+pap,regionterms[region],cid=old)
for region in ['interatrial-septum','atrioventricular-septum']:group(region,region.replace('-',' ').capitalize(),[prefix+'myocardium:'+region],regionterms[region])
group('interventricular-septum','Interventricular septum',[prefix+'myocardium:interventricular-muscular-septum',prefix+'fibrous:interventricular-membranous-septum'],'A12.1.00.013')
for kind,region in [('myocardium','interventricular-muscular-septum'),('fibrous','interventricular-membranous-septum')]:group(region,region.replace('-',' ').capitalize(),[prefix+kind+':'+region],regionterms[region],'interventricular-septum')
for kind in ['myocardium','endocardium','epicardium']:group(kind,kind.capitalize(),[pid for name,pid in wall['part_ids'].items()if name.split(':')[0]==kind],layerterms[kind])
valves={
 'tricuspid':('Tricuspid valve','A12.1.02.003',[prefix+'tricuspid-anterior-cusp','ZA:Inferior leaflet of right atrioventricular valve','ZA:Septal leaflet of right atrioventricular valve'],'right-ventricle'),
 'mitral':('Mitral valve','A12.1.04.003',[prefix+'mitral-anterior-cusp','ZA:Posterior leaflet of left atrioventricular valve'],'left-ventricle'),
 'aortic':('Aortic valve','A12.1.04.012',['ZA:Right coronary leaflet','ZA:Left coronary leaflet','ZA:Non-coronary leaflet'],'left-ventricle'),
 'pulmonary':('Pulmonary valve','A12.1.02.010',['ZA:Right semilunar leaflet of pulmonary valve','ZA:Left semilunar leaflet of pulmonary valve','ZA:Anterior semilunar leaflet of pulmonary valve'],'right-ventricle')}
for key,(name,term,ids,parent)in valves.items():
 group(key+'-valve',name,ids,term,parent)
 for pid in ids:
  c=concepts[parts[pid]['conceptId']];term=c.get('ta98Term')
  if term:group(key+'-leaflet-'+str(len(groups)),parts[pid]['name'],[pid],term,key+'-valve',c.get('ta98Kind','exact'),cid=c['id'])
for pid in retained:
 if 'papillary'not in parts[pid]['name'].lower():continue
 c=concepts[parts[pid]['conceptId']];term=c.get('ta98Term')
 if term:group('papillary-'+str(len(groups)),c['name'],[pid],term,'left-ventricle'if'left ventricle'in c['name'].lower()else'right-ventricle',c.get('ta98Kind','exact'),c['id'])
roots=['ZA:Ascending aorta','ZA:Pulmonary trunk','ZA:Superior vena cava','ZA:Inferior vena cava (thoracic part)','ZA:Left inferior pulmonary vein','ZA:Left superior pulmonary vein','ZA:Right inferior pulmonary vein','ZA:Right superior pulmonary vein'];roots=[('TA98HEPATIC:male:region:A12.3.09.001'if pid=='ZA:Inferior vena cava (thoracic part)'and parts[pid].get('suppressed')and 'TA98HEPATIC:male:region:A12.3.09.001'in parts and not parts['TA98HEPATIC:male:region:A12.3.09.001'].get('suppressed')else pid)for pid in roots];wholeids+=roots;concept(whole,'Heart',wholeids,'A12.1.00.001');next(g for g in groups if g['key']=='heart')['partIds']=wholeids;a['finiteCardiacIntegration']=dict(schema='cardiac-local-finite-v1',sex='male',wallPartIds=wallids,valvePartIds=[pid for pid in retained+newcusps if'leaflet'in parts[pid]['name'].lower()],papillaryPartIds=[pid for pid in retained if'leaflet'not in parts[pid]['name'].lower()],rootPartIds=roots,wholeConceptId=whole,termGroups=groups,sourceCatalogueSha256=digest_file(src),wallBankSha256=wall['bank_sha256'],remainingLimits=limits,anatomicalAcceptance=False);a['defaultHidden']=sorted((set(a.get('defaultHidden',[]))|set(native))-set(wholeids));a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts']);atomic_json(out/'atlas.json',a);atomic_json(out/'construction.json',dict(prepared=True,source_catalogue_sha256=digest_file(src),wall_bank_sha256=wall['bank_sha256'],av_bank_sha256=av['bank_sha256'],part_ids=wholeids,changed_parts=[p['id']for p in a['parts']if p!=next((q for q in before['parts']if q['id']==p['id']),None)],term_groups=len(groups),distinct_ta98_terms=len({g['term']for g in groups}),installed=False));print('Heart catalogue prepared:',len(wholeids),'parts;',len(groups),'controls;',len({g['term']for g in groups}),'TA98 terms',flush=True)
