"""Stage the audited regional cord without changing any unrelated anatomy.

Native whole/segment controls remain stable. New microscopic territories are
finite owned reference regions, explicitly qualified rather than donor
histology. Existing identifiers are updated in place, never duplicated.
"""
from pathlib import Path
from collections import defaultdict
import argparse,copy,json
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser()
for option in ['bank','source','context','review','output']:ap.add_argument('--'+option,required=True)
ap.add_argument('--namespace',default='regional',help='New physical ownership namespace; anatomical concept identifiers remain stable')
ap.add_argument('--checkpoints',help='Source-bound completed audit stages for this bank')
args=ap.parse_args();base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
bank=Path(args.bank).resolve();source=Path(args.source).resolve();out=Path(args.output).resolve()
out.mkdir(parents=True,exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
a=json.loads(path.read_text());old=copy.deepcopy(a['finiteSpinalIntegration']);assert old['schema']in ['spinal-cord-native-core/v3','spinal-cord-native-core/v4']
assert args.namespace.replace('-','').isalnum()
m=json.loads((bank/'tissue.json').read_text());meta=json.loads((source/'cells.json').read_text())
assert m['passed']and meta['passed']and m['source_sha256']==digest_file(source/'cells.npz')
# A finer partition is a new physical ownership revision. Stable anatomical
# concepts select it; historical mesh identifiers retain their literal source
# geometry, including the old interface-forest evidence.
aliases={p['id']:'TA98CORD:female:'+args.namespace+':'+name for name,p in m['parts'].items()}
assert not(set(aliases.values())&set(old['partIds']+old['coveringPartIds'])),'A new ownership revision must retain the historical physical identifiers'
for name,part in m['parts'].items():
    part['id']=aliases[part['id']]
    for group in part['tissueDrawGroups']:
        if group.get('neighbor'):group['neighbor']=aliases[group['neighbor']]
m['part_ids']={name:part['id']for name,part in m['parts'].items()}
m['source_part_aliases']=aliases
context=Path(args.context).resolve();screen=json.loads(context.read_text())
assert screen['passed']and screen['actual_packed_float32']and screen['source_cell_sha256']==m['bank_sha256']
review=Path(args.review).resolve();visual=json.loads(review.read_text())
assert visual['source_sha256']==digest_file(source/'cells.npz')and len(visual['views'])>=4
plan=json.loads((base/'female-source-v1/inventory.json').read_text())
authority={r['authority']['id']:r['authority']for r in plan['terms']}
attributes=meta['attributes'];new_ids={p['id']for p in m['parts'].values()}
old_ids=set(old['partIds']+old['coveringPartIds']);retired=old_ids-new_ids
by_id={p['id']:p for p in a['parts']};assert len(by_id)==len(a['parts'])
concepts={c['id']:c for c in a['concepts']};macro=defaultdict(list);regions=defaultdict(list)
cell_regions=[];core=[];covers=[];new_terms=defaultdict(list);primary=defaultdict(list)
chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=m['total_bytes']))
for pid in retired:by_id[pid].update(suppressed=True,defaultHidden=True)
for name,row in m['parts'].items():
    info=attributes[name];native=info['native_region'];segment=info['axial_region'];family=native.split('-')[0]
    pid=row['id'];base_name=name.split('@')[0].split('~')[0];cid='TA98CORD:female:cell:'+base_name
    macro_cid='TA98CORD:female:'+native
    if family=='pia':term='A14.1.01.309';label='Spinal pia mater'
    elif info.get('ta98_terms'):
        term=info['ta98_terms'][0];label=authority[term]['name_en'].capitalize()
    else:
        term=concepts[macro_cid]['ta98Term'];label=authority[term]['name_en'].capitalize()
    if base_name.endswith('-left'):label+=' (left)'
    elif base_name.endswith('-right'):label+=' (right)'
    material='spinal:female:'+family;color={'white':[232,214,185],'grey':[158,122,115],'pia':[206,171,164]}[family]
    appearance=dict(color=[v/255 for v in color],pattern='fibrous'if family=='pia'else'solid',scale=.001,relief=.01 if family=='pia'else 0.)
    a['materials'][material]=dict(color=color,roughness=.68,metalness=0,tissueAppearance=appearance)
    part=dict(row,conceptId=cid,name=label+' · '+segment,system='nervous',chunk=chunk,material=material,
        sourceId='spinal_cord:'+name,ta98Term=term,pairedTissue=True,sectionSolid=True,
        sectionAssembly='spinal-cord:female',tissueWallAssembly='spinal-cord:female',tissueWallLayer=family,tissueAppearance=appearance,
        reviewStatus='Finite shared regional reference. Nuclear and tract borders are anatomical estimates; root interfaces remain unfinished.')
    if pid in by_id:by_id[pid].clear();by_id[pid].update(part)
    else:a['parts'].append(part);by_id[pid]=part
    (covers if family=='pia'else core).append(pid)
    primary[cid].append(pid);macro[macro_cid].append(pid)
    cell_regions.append(dict(partId=pid,segment=segment,nativeRegion=native))
    if family!='pia':regions[segment].append(pid)
    for t in info.get('ta98_terms',[]):
        ancestors=set()
        while t in authority and t not in ancestors:
            ancestors.add(t);new_terms[t].append(pid);t=authority[t]['parent_id']
    concepts[cid]=dict(id=cid,name=label,elements=[],componentElements=[],ta98Term=term,ta98Kind='exact',
        partOf=macro_cid if family!='pia'else old['wholeConceptId'],
        coverageLimitation='Reference tissue territory; microscopic borders and donor-specific cytoarchitecture are not verified.'if info.get('reference_estimate')else'Native-envelope regional tissue. Root sleeves and cranial/root continuity are unfinished.')
for cid,ids in primary.items():concepts[cid].update(elements=ids,componentElements=ids)
for cid,ids in macro.items():
    if cid in concepts:concepts[cid].update(elements=ids,componentElements=ids)
groups=[]
old_regions={r['partId']:(r['nativeRegion'],r['segment'])for r in old.get('cellRegions',[])}
old_roles=defaultdict(set)
for pid,role in old_regions.items():old_roles[role].add(pid)
for group in old['termGroups']:
    ids=[];selected=set(group['partIds'])
    complete_roles={role for role,members in old_roles.items()if members<=selected}
    old_logical={by_id[pid].get('sourceId','').removeprefix('spinal_cord:').split('~')[0].split('@')[0]+'@'+pid.split('@')[-1]for pid in selected}
    for name,row in m['parts'].items():
        info=attributes[name];prefix='TA98CORD:female:'+info['native_region']+'@'+info['axial_region']
        logical=name.split('@')[0].split('~')[0]+'@'+info['axial_region']
        if(old['schema'].endswith('/v3')and prefix in selected)or(info['native_region'],info['axial_region'])in complete_roles or logical in old_logical or group['term']in info.get('ta98_terms',[]):ids.append(row['id'])
    assert ids,group['term'];concepts[group['conceptId']].update(elements=ids,componentElements=ids)
    groups.append(dict(group,partIds=ids))
for term,ids in sorted(new_terms.items()):
    assert term in authority;cid='TA98CORD:female:term:'+term
    ids=list(dict.fromkeys(ids))
    existing=next((g for g in groups if g['term']==term),None)
    if existing:
        # A represented parent retains all of its native tissue, rather than
        # being reduced to only the newly named microscopic children.
        combined=list(dict.fromkeys([*existing['partIds'],*ids]));existing['partIds']=combined
        concepts[existing['conceptId']].update(elements=combined,componentElements=combined)
        continue
    concepts[cid]=dict(id=cid,name=authority[term]['name_en'].capitalize(),elements=ids,componentElements=ids,
        ta98Term=term,ta98Kind='exact',partOf=old['wholeConceptId'],
        coverageLimitation='Finite anatomical reference territory; microscopic borders are estimated. Eligibility follows the recorded level ranges, not donor histology.')
    groups.append(dict(term=term,conceptId=cid,partIds=ids))
segments=[]
for prior in old['segments']:
    ids=regions[prior['name']];concepts[prior['conceptId']].update(elements=ids,componentElements=ids)
    segments.append(dict(prior,partIds=ids,coveringPartIds=['TA98CORD:female:'+args.namespace+':pia@'+prior['name']]))
concepts[old['wholeConceptId']].update(elements=core,componentElements=core)
a['concepts']=list(concepts.values());a['defaultHidden']=[pid for pid in dict.fromkeys([*a.get('defaultHidden',[]),*sorted(retired)])if pid not in new_ids]
count=sum(g['term'].startswith('A14.1.02.')for g in groups)
a['finiteSpinalIntegration']=dict(old,schema='spinal-cord-native-core/v4',partIds=core,coveringPartIds=covers,
    bankSha256=m['bank_sha256'],termGroups=groups,segments=segments,cellRegions=cell_regions,regionalGroups=dict(regions),
    representedCordTerms=count,remainingCordTerms=119-count,anatomicallyAccepted=False,
    retiredPartIds=sorted((set(old.get('retiredPartIds',[]))|retired)-new_ids),
    courseFit=dict(source=digest_file(source/'cells.npz'),bank_sha256=m['bank_sha256'],
        fixed_skeleton_screen_sha256=digest_file(context),source_fit='Shared positive radial volume carrier; every tissue chart conforms to its affine volume pieces'),
    remaining=[str(119-count)+' cord terms remain unrepresented','Nuclear and tract borders are anatomical reference estimates',
        'Dura/arachnoid root ports and sleeves, filum and denticulate attachments','Verify native root/segment and cranial interfaces',
        'Several PAM50 tract fields end at S1 because the reference data stop there; this is not a physiological termination claim'])
a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts']if not p.get('suppressed'))
atomic_json(out/'atlas.json',a)
inventory=dict(schema='ta98-geometry-inventory/v1',model='local-female-ta98',model_sex='female',root_term='A14.1.02.001',
    source_catalogue=str(path),source_catalogue_sha256=digest_file(path),authority_path=str(ROOT/plan['authority']['path']),authority_sha256=plan['authority']['sha256'],
    terms=[dict(r,source_meshes=r['source_meshes']['female'],applicable_to_model=True,
        candidate_part_ids=next((g['partIds']for g in groups if g['term']==r['authority']['id']),[]))for r in plan['terms']],term_count=len(plan['terms']),anatomical_acceptance=False)
atomic_json(out/'inventory.json',inventory)
rows={p['id']:p for p in m['parts'].values()};layout=dict(m,parts=rows,origin=[0.,1.35,-.05]);atomic_json(out/'audited-layout.json',layout)
interfaces=[]
for pid,part in rows.items():
    for g in part['tissueDrawGroups']:
        peer=g.get('neighbor')
        if not peer or pid>=peer:continue
        interfaces.append(dict(kind='shared-interface',first=pid,second=peer,parameters=dict(
            first_faces=list(range(g['start']//3,(g['start']+g['count'])//3)),
            second_faces=[i for gg in rows[peer]['tissueDrawGroups']if gg.get('neighbor')==pid for i in range(gg['start']//3,(gg['start']+gg['count'])//3)])))
spec=dict(model='local-female-ta98',proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(path),
    changed_parts=sorted(retired|new_ids),closed_parts=sorted(new_ids),inventory=str(out/'inventory.json'),
    max_pairs=250000,max_tests=2000000,tissue_audit_layouts=[str(out/'audited-layout.json')],
    tissue_audit_checkpoints=str(Path(args.checkpoints).resolve()if args.checkpoints else base/'female-regional-carrier-bank-checks-v1'),interface_checks=interfaces,
    assembly=dict(whole_organ=old['wholeConceptId'],part_ids=sorted(new_ids),
        fixed=['Entire native skeleton','Every unrelated organ, peripheral nerve and vessel','Entire male model'],
        adjoining_interfaces=[dict(kind='shared-interface',scope='Reciprocal radial/axial tissue interfaces, stored once'),
            dict(kind='open-space',scope='Central canal remains a void; no independent lumen plug'),
            dict(kind='unfinished-native-attachment',scope='Root sleeves and cranial continuation are unfinished')],remaining_limits=a['finiteSpinalIntegration']['remaining']),
    evidence_files=[str(bank/'tissue.json'),str(source/'cells.json'),str(context),str(review)])
atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec))
print(json.dumps(dict(cord_terms=count,missing=119-count,finite_cells=len(new_ids),duplicate_part_ids=0,bytes=m['total_bytes'])),flush=True)
