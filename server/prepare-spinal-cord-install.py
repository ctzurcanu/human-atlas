"""Prepare a local native-frame cord core without changing unrelated anatomy."""
from pathlib import Path
import argparse,copy,json
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';before=json.loads(source.read_text());a=copy.deepcopy(before)
plan_path=base/'female-source-v1/inventory.json';plan=json.loads(plan_path.read_text());assert digest_file(source)==plan['source_catalogues']['female']['sha256']
bank=base/'female-region-bank-v2';manifest=json.loads((bank/'tissue.json').read_text());assert manifest['passed']
source_id='LOCAL:female:spinal_cord:864';old=next(p for p in a['parts']if p['id']==source_id);whole_id=old['conceptId']
old.update(suppressed=True,defaultHidden=True);a['defaultHidden']=list(dict.fromkeys([*a.get('defaultHidden',[]),source_id]))
index=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=(bank/'cells.bin').stat().st_size))
names={'grey-anterior':('Anterior horn of spinal cord','A14.1.02.103'),
    'grey-posterior':('Posterior horn of spinal cord','A14.1.02.115'),
    'grey-intermediate':('Intermediate column of spinal cord','A14.1.02.130'),
    'white-anterior':('Anterior funiculus','A14.1.02.202'),
    'white-lateral':('Lateral funiculus','A14.1.02.215'),
    'white-posterior':('Posterior funiculus','A14.1.02.239')}
part_ids=[];concepts=[]
for key,record in manifest['parts'].items():
    region,side=key.rsplit('-',1);label,term=names[region];pid=record['id'];part_ids.append(pid)
    color=[158,122,115]if key.startswith('grey')else[232,214,185];family=key.split('-')[0];material='spinal:female:'+family
    appearance=dict(color=[v/255 for v in color],pattern='solid',scale=.001,relief=0.)
    p=dict(record,conceptId=pid,name=label+' ('+side+')',system='nervous',chunk=index,material=material,
        sourceId='spinal_cord:'+key,ta98Term=term,pairedTissue=True,sectionSolid=True,
        sectionAssembly='spinal-cord:female',tissueWallAssembly='spinal-cord:female',tissueWallLayer=family,tissueAppearance=appearance,
        reviewStatus='Native-envelope macroscopic tissue; regional boundaries provisional. Nuclei and individual tracts remain unbuilt.')
    a['materials'][material]=dict(color=color,roughness=.68,metalness=0,tissueAppearance=appearance)
    a['parts'].append(p)
    concepts.append(dict(id=pid,name=p['name'],elements=[pid],ta98Term=term,ta98Kind='exact',partOf=whole_id,
        coverageLimitation='Reference-guided macroscopic horn/funiculus region inside the native female cord. This does not establish microscopic nuclei, laminae or individual tract boundaries.'))
whole=next(c for c in a['concepts']if c['id']==whole_id)
whole.update(elements=part_ids,componentElements=part_ids,ta98Term='A14.1.02.001',ta98Kind='exact')
a['concepts'].extend(concepts)
groups=[dict(term='A14.1.02.001',conceptId=whole_id,partIds=part_ids)]
def group(term,label,selected):
    cid='TA98CORD:female:term:'+term
    a['concepts'].append(dict(id=cid,name=label,elements=selected,componentElements=selected,ta98Term=term,ta98Kind='exact',partOf=whole_id,
        coverageLimitation='Collective control of the same stored macroscopic tissue; no additional overlapping surface or microscopic tissue is implied.'))
    groups.append(dict(term=term,conceptId=cid,partIds=selected))
grey=[p for p in part_ids if ':grey-'in p];white=[p for p in part_ids if ':white-'in p]
for code,label,selected in [
    ('020','Grey substance',grey),('101','Grey columns',grey),('024','White substance',white),('201','White substance',white),('013','Funiculi of spinal cord',white),
    ('021','Anterior horn',[p for p in grey if ':grey-anterior-'in p]),('102','Anterior column',[p for p in grey if ':grey-anterior-'in p]),('103','Anterior horn',[p for p in grey if ':grey-anterior-'in p]),
    ('023','Posterior horn',[p for p in grey if ':grey-posterior-'in p]),('114','Posterior column',[p for p in grey if ':grey-posterior-'in p]),('115','Posterior horn',[p for p in grey if ':grey-posterior-'in p]),
    ('130','Intermediate column',[p for p in grey if ':grey-intermediate-'in p]),
    ('202','Anterior funiculus',[p for p in white if ':white-anterior-'in p]),('215','Lateral funiculus',[p for p in white if ':white-lateral-'in p]),('239','Posterior funiculus',[p for p in white if ':white-posterior-'in p])]:group('A14.1.02.'+code,label,selected)
a['finiteSpinalIntegration']=dict(schema='spinal-cord-native-core/v1',sex='female',sourcePartId=source_id,wholeConceptId=whole_id,
    partIds=part_ids,bankSha256=manifest['bank_sha256'],sourcePlanSha256=digest_file(plan_path),termGroups=groups,
    representedCordTerms=len(groups),totalCordTerms=119,remainingCordTerms=119-len(groups),anatomicallyAccepted=False,
    remaining=['Segment regions, enlargements and conus controls','Complete finite meninges and root sleeves','Named nuclei, laminae and detailed tract territories'])
a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts']if not p.get('suppressed'))
atomic_json(out/'atlas.json',a)
inventory=dict(schema='ta98-geometry-inventory/v1',model='local-female-ta98',model_sex='female',root_term='A14.1.02.001',
    source_catalogue=str(source),source_catalogue_sha256=digest_file(source),authority_path=str(ROOT/plan['authority']['path']),authority_sha256=plan['authority']['sha256'],
    terms=[dict(r,source_meshes=r['source_meshes']['female'],applicable_to_model=True,
        candidate_part_ids=next((g['partIds']for g in groups if g['term']==r['authority']['id']),[]))for r in plan['terms']],
    term_count=len(plan['terms']),anatomical_acceptance=False)
atomic_json(out/'inventory.json',inventory)
layout=dict(manifest,parts={p['id']:p for p in manifest['parts'].values()},origin=[0.,1.35,-.05]);atomic_json(out/'audited-layout.json',layout)
interfaces=[]
for p in manifest['parts'].values():
    for g in p['tissueDrawGroups']:
        peer=g.get('neighbor')
        if not peer or p['id']>=peer:continue
        q=layout['parts'][peer]
        interfaces.append(dict(kind='shared-interface',first=p['id'],second=peer,parameters=dict(
            first_faces=list(range(g['start']//3,(g['start']+g['count'])//3)),
            second_faces=[i for gg in q['tissueDrawGroups']if gg.get('neighbor')==p['id']for i in range(gg['start']//3,(gg['start']+gg['count'])//3)])))
spec=dict(model='local-female-ta98',proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(source),
    changed_parts=[source_id,*part_ids],closed_parts=part_ids,inventory=str(out/'inventory.json'),max_pairs=250000,max_tests=2000000,
    tissue_audit_layouts=[str(out/'audited-layout.json')],tissue_audit_checkpoints=str(bank.with_name(bank.name+'-checks')),interface_checks=interfaces,
    assembly=dict(whole_organ=whole_id,part_ids=part_ids,fixed=['All skeleton','All other native female organs, nerves and vessels','Entire male model'],
        adjoining_interfaces=[dict(kind='shared-interface',scope='All stored grey/white and horn/funiculus tissue surfaces; reciprocal faces and normal words'),
            dict(kind='open-space',scope='Central canal remains a void; no solid plug'),
            dict(kind='unchanged-root-attachments',scope='45 native sources measured against the retained envelope',evidence=str(base/'female-surface-v1/native-root-attachments.json'))],
        remaining_limits=a['finiteSpinalIntegration']['remaining']),
    evidence_files=[str(bank/'tissue.json'),str(base/'female-surface-v1/surface.json'),str(base/'female-surface-v1/native-root-attachments.json')])
atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec))
print(json.dumps(dict(candidate=str(out/'atlas.json'),parts=12,cord_term_controls=len(groups),unrepresented_cord_terms=119-len(groups),bytes=manifest['total_bytes'],male_unchanged=True)),flush=True)
