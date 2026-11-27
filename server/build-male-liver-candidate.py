"""Prepare, but never install, the reviewed local male liver assembly.

The source catalogue and emitted geometry are hash-bound. Original variants
remain literal and hidden. Redistribution and clinical acceptance are outside
this local delivery.
"""
from pathlib import Path
import argparse,copy,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--workbench',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
 b=args.workbench.resolve();out=args.output.resolve();out.mkdir(exist_ok=False)
 source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
 original=json.loads(source.read_bytes());expected=json.loads((b/'source-v1/source.json').read_text())['catalogue_sha256'];assert digest_file(source)==expected,'Current catalogue changed'
 assembly=json.loads((b/'shared-lobes-v54/assembly.json').read_text());checks=json.loads((b/'assembly-checks-v56/evidence.json').read_text());assert checks['passed']and checks['prisms_strictly_inward'];assert checks['bank_sha256']==assembly['bank_sha256']==digest_file(assembly['bank'])
 for name in assembly['parts']:
  q=json.loads((b/'shared-audits-v55'/(name+'.json')).read_text());assert q['passed']and not q['self_intersections']['defect_pairs']
 a=copy.deepcopy(original);old={p['id']:p for p in a['parts']};native=old['ZA:Liver:Organ'];source_rows=json.loads((b/'source-v1/source.json').read_text())['sources'];retired=[r['id']for r in source_rows if r['id'].startswith('ZA:')and ('segment of liver ('in r['id']or r['id']in ['ZA:Liver:Organ','ZA:Liver:Ligament'])];ids=assembly['part_ids'];all_ids=list(ids.values());chunk=len(a['chunks']);bank=Path(assembly['bank']);url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'));a['chunks'].append(dict(url=url,bytes=bank.stat().st_size))
 for uid in retired:
  old[uid].update(suppressed=True,defaultHidden=True,reviewStatus='retained literal source; replaced by locally reviewed liver lobes and capsule')
  if uid in ['ZA:Liver:Organ','ZA:Liver:Ligament']:old[uid]['replacementParts']=all_ids
 a['defaultHidden']=sorted(set(a.get('defaultHidden',[]))|set(retired))
 terms={'right':'A05.8.01.026','left':'A05.8.01.027','quadrate':'A05.8.01.029','caudate':'A05.8.01.030','capsule':'A05.8.01.054'}
 names={'right':'Right lobe of liver','left':'Left lobe of liver','quadrate':'Quadrate lobe of liver','caudate':'Caudate lobe of liver','capsule':'Fibrous capsule of liver'}
 limitation='Local liver fit with four provisional gross lobes and a finite capsule. Porta openings are retained. Functional territories, microscopic structures and patent vascular/biliary wall interfaces remain unfinished.'
 for name,row in assembly['parts'].items():
  a.setdefault('materials',{})['LIVER:male:'+name]=dict(color=[107,61,43],roughness=.68,metalness=0)
  p={**copy.deepcopy(row),'chunk':chunk,'conceptId':ids[name],'name':names[name],'sourceId':'Liver:'+name,'system':'digestive','groups':['8: Visceral systems','Digestive system','Abdomen','Trunk','Liver'],'material':'LIVER:male:'+name,'license':native['license'],'provenance':copy.deepcopy(native['provenance']),'ta98Term':terms[name],'tissueWallAssembly':'liver:male:lobes-v1','tissueWallLayer':name,'pairedTissue':True,'sectionSolid':name=='capsule','coverageLimitation':limitation,'description':limitation,'reviewStatus':'local-liver-fit-and-four-lobe-preview; vascular-biliary-interfaces-unfinished','tissueAppearance':{'color':[.58,.25,.23]if name!='capsule'else[.61,.29,.26],'pattern':'solid'if name!='capsule'else'fibrous','relief':0 if name!='capsule'else .000002,'scale':.0004}}
  a['parts'].append(p)
  a['concepts'].append(dict(id=ids[name],name=names[name],elements=[ids[name]],ta98Term=terms[name],ta98Kind='exact',partOf='ZA:Liver'))
 for c in a['concepts']:
  if c['id']in ['ZA:Liver','ZA:collection:Liver']:
   c.update(elements=all_ids,name='Liver',ta98Term='A05.8.01.001',ta98Kind='exact')
 a['liverIntegration']=dict(schema='male-liver-lobes-capsule/v1',sourceCatalogueSha256=expected,sourceParts=retired,partIds=all_ids,wholeConceptId='ZA:Liver',termGroups=[dict(term='A05.8.01.001',conceptId='ZA:Liver',parts=all_ids)]+[dict(term=terms[n],conceptId=uid,parts=[uid])for n,uid in ids.items()],bankSha256=assembly['bank_sha256'],bankUrl=url,remainingLimits=limitation,anatomicallyAccepted=False)
 proposed=out/'atlas.json';atomic_json(proposed,a);changed=[p['id']for p in a['parts']if p['id']not in old or p!=next(q for q in original['parts']if q['id']==p['id'])]
 interfaces=[]
 for name,row in assembly['parts'].items():
  for group in row['tissueDrawGroups']:
   peer=group.get('neighbor')
   if not peer or ids[name]>=peer:continue
   other=next(q for n,q in assembly['parts'].items()if ids[n]==peer);opposite=next(g for g in other['tissueDrawGroups']if g.get('neighbor')==ids[name]);interfaces.append(dict(kind='shared-interface',first=ids[name],second=peer,parameters=dict(first_faces=list(range(group['start']//3,(group['start']+group['count'])//3)),second_faces=list(range(opposite['start']//3,(opposite['start']+opposite['count'])//3)))))
 spec=dict(model='local-ta98',proposed_catalogue=str(proposed),expected_catalogue_sha256=expected,changed_parts=changed,inventory=str(b/'source-v1/inventory.json'),require_closed=False,interface_checks=interfaces,assembly=dict(whole_organ='ZA:Liver',part_ids=all_ids,adjoining_interfaces=[dict(name='fixed diaphragm, skeleton, stomach, lungs and gallbladder/renal/adrenal',evidence=str(b/'assembly-checks-v56/evidence.json')),dict(name='porta, hepatic vessels and biliary ducts',evidence=str(b/'hilar-interfaces-v34/evidence.json'),status='source mismatches retained; no claimed patent geometric junction')],remaining_limits=limitation),evidence_files=[str(b/'surface-section-review-v57/evidence.json'),str(b/'assembly-checks-v56/evidence.json'),str(b/'shared-lobes-v53/evidence.json')])
 atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec));atomic_json(out/'candidate.json',dict(source_sha256=expected,proposed_sha256=digest_file(proposed),changed_parts=changed,bank_sha256=assembly['bank_sha256'],bytes=assembly['total_bytes'],external_deployment=False))
 print(json.dumps(dict(proposed=str(proposed),changed_parts=len(changed),bank_bytes=assembly['total_bytes'])))

if __name__=='__main__':main()
