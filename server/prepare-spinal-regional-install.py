"""Prepare reviewed regional cord tissues and pia for the ordinary local atlas."""
from pathlib import Path
import argparse,copy,json
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
before=json.loads(source.read_text());a=copy.deepcopy(before);old=a['finiteSpinalIntegration'];whole=old['wholeConceptId']
bank=base/'female-regional-pia-bank-v1';m=json.loads((bank/'tissue.json').read_text());assert m['passed']
meta=json.loads((base/'female-pia-v4/cells.json').read_text());attributes=meta['attributes']
old_ids=set(old['partIds'])
for p in a['parts']:
    if p['id']in old_ids:p.update(suppressed=True,defaultHidden=True)
a['defaultHidden']=list(dict.fromkeys([*a.get('defaultHidden',[]),*sorted(old_ids)]))
chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=m['total_bytes']))
concepts={c['id']:c for c in a['concepts']};core=[];covers=[];macro={};regional={}
for name,row in m['parts'].items():
    info=attributes[name];family=info['native_region'].split('-')[0];pid=row['id']
    if family=='pia':
        cid=pid;label='Spinal pia mater';term='A14.1.01.309';color=[206,171,164];covers.append(pid)
        concepts[cid]=dict(id=cid,name=label,elements=[pid],ta98Term=term,ta98Kind='exact',
            coverageLimitation='Finite pial tissue follows the native cord exterior. Root sleeves, filum attachment and microfenestrations remain unfinished.')
    else:
        native=info['native_region'];cid='TA98CORD:female:'+native;original=concepts[cid];label=original['name']+' · '+info['axial_region'].replace('-',' ')
        term=original['ta98Term'];core.append(pid);macro.setdefault(cid,[]).append(pid);regional.setdefault(info['axial_region'],[]).append(pid)
        color=[158,122,115]if family=='grey'else[232,214,185]
    material='spinal:female:'+family;appearance=dict(color=[x/255 for x in color],pattern='fibrous'if family=='pia'else'solid',scale=.001,relief=.01 if family=='pia'else 0.)
    a['materials'][material]=dict(color=color,roughness=.68,metalness=0,tissueAppearance=appearance)
    p=dict(row,conceptId=cid,name=label,system='nervous',chunk=chunk,material=material,sourceId='spinal_cord:'+name,
        ta98Term=term,pairedTissue=True,sectionSolid=True,sectionAssembly='spinal-cord:female',
        tissueWallAssembly='spinal-cord:female',tissueWallLayer=family,tissueAppearance=appearance,
        reviewStatus='Root-registered finite regional tissue; distal segment boundaries and microscopic territories remain provisional.')
    a['parts'].append(p)
for cid,ids in macro.items():concepts[cid].update(elements=ids,componentElements=ids)
groups=[]
for group in old['termGroups']:
    selected=[pid for prior in group['partIds']for pid in macro[prior]]
    concepts[group['conceptId']].update(elements=selected,componentElements=selected)
    groups.append(dict(group,partIds=selected))
def add(term,name,ids,limitation):
    cid='TA98CORD:female:term:'+term
    concepts[cid]=dict(id=cid,name=name,elements=ids,componentElements=ids,ta98Term=term,ta98Kind='exact',
        **({'partOf':whole}if term.startswith('A14.1.02.')else{}),coverageLimitation=limitation)
    groups.append(dict(term=term,conceptId=cid,partIds=ids))
def selected(*bands):return [pid for key in bands for pid in regional[key]]
boundary='Root-registered approximate region; caudal boundaries use published relative segment lengths and still need native root verification.'
for code,name,bands in [
    ('014','Cervical part',('cervical-upper','cervical-enlargement')),
    ('015','Thoracic part',('thoracic-t1','thoracic')),
    ('016','Lumbar part',('lumbar',)),('017','Sacral part',('sacral-proximal','sacral-distal')),
    ('018','Coccygeal part',('coccygeal',)),
    ('002','Cervical enlargement',('cervical-enlargement','thoracic-t1')),
    ('003','Lumbosacral enlargement',('lumbar','sacral-proximal')),
    ('004','Conus medullaris',('sacral-proximal','sacral-distal','coccygeal'))]:add('A14.1.02.'+code,name,selected(*bands),boundary)
for code,name in [('001','Meninges'),('003','Leptomeninx'),('301','Pia mater'),('309','Spinal pia mater')]:
    add('A14.1.01.'+code,name,covers,'Pial component only; dura, arachnoid, root sleeves and caudal continuations remain unfinished.'if code in['001','003']else'Finite pial coat; root continuations and microscopic fenestrations are not yet represented.')
a['concepts']=list(concepts.values());concepts[whole].update(elements=core,componentElements=core)
terms=sum(g['term'].startswith('A14.1.02.')for g in groups)
a['finiteSpinalIntegration']=dict(schema='spinal-cord-native-core/v2',sex='female',sourcePartId=old['sourcePartId'],wholeConceptId=whole,
    partIds=core,coveringPartIds=covers,retiredPartIds=sorted(old_ids),bankSha256=m['bank_sha256'],
    sourcePlanSha256=digest_file(base/'female-source-v1/inventory.json'),termGroups=groups,regionalGroups=regional,
    representedCordTerms=terms,totalCordTerms=119,remainingCordTerms=119-terms,anatomicallyAccepted=False,
    remaining=['Complete dural/arachnoid sac, root sleeves and filum','Verify caudal segment/root correspondence','Named nuclei, laminae and individual tract territories'])
a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts']if not p.get('suppressed'))
atomic_json(out/'atlas.json',a)
plan=json.loads((base/'female-source-v1/inventory.json').read_text())
inventory=dict(schema='ta98-geometry-inventory/v1',model='local-female-ta98',model_sex='female',root_term='A14.1.02.001',
    source_catalogue=str(source),source_catalogue_sha256=digest_file(source),authority_path=str(ROOT/plan['authority']['path']),authority_sha256=plan['authority']['sha256'],
    terms=[dict(r,source_meshes=r['source_meshes']['female'],applicable_to_model=True,
        candidate_part_ids=next((g['partIds']for g in groups if g['term']==r['authority']['id']),[]))for r in plan['terms']],term_count=len(plan['terms']),anatomical_acceptance=False)
atomic_json(out/'inventory.json',inventory)
layout=dict(m,parts={p['id']:p for p in m['parts'].values()},origin=[0.,1.35,-.05]);atomic_json(out/'audited-layout.json',layout)
interfaces=[]
for p in m['parts'].values():
    for g in p['tissueDrawGroups']:
        peer=g.get('neighbor')
        if not peer or p['id']>=peer:continue
        q=layout['parts'][peer]
        interfaces.append(dict(kind='shared-interface',first=p['id'],second=peer,parameters=dict(
            first_faces=list(range(g['start']//3,(g['start']+g['count'])//3)),
            second_faces=[i for gg in q['tissueDrawGroups']if gg.get('neighbor')==p['id']for i in range(gg['start']//3,(gg['start']+gg['count'])//3)])))
spec=dict(model='local-female-ta98',proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(source),
    changed_parts=[*sorted(old_ids),*core,*covers],closed_parts=[*core,*covers],inventory=str(out/'inventory.json'),max_pairs=250000,max_tests=2000000,
    tissue_audit_layouts=[str(out/'audited-layout.json')],tissue_audit_checkpoints=str(bank.with_name(bank.name+'-checks')),interface_checks=interfaces,
    assembly=dict(whole_organ=whole,part_ids=[*core,*covers],fixed=['Entire skeleton','All unrelated organs and peripheral nerves','Entire male model'],
        adjoining_interfaces=[dict(kind='shared-interface',scope='Reciprocal tissue charts and pial inner interface'),dict(kind='open-space',scope='Central canal and pial continuation rims remain unfilled')],
        remaining_limits=a['finiteSpinalIntegration']['remaining']),evidence_files=[str(bank/'tissue.json'),str(base/'female-root-regions-v2/root-regions.json'),str(base/'female-pia-v4/cells.json')])
atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec))
print(json.dumps(dict(cord_terms=terms,remaining_cord_terms=119-terms,finite_parts=len(core)+len(covers),bytes=m['total_bytes'],interfaces=len(interfaces))),flush=True)
