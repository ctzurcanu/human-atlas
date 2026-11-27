"""Prepare 31 disjoint segment controls in both hierarchies, without term inflation."""
from pathlib import Path
import argparse,copy,json
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--bank',required=True);args=ap.parse_args()
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';before=json.loads(source.read_text());a=copy.deepcopy(before);old=a['finiteSpinalIntegration']
assert old['schema']=='spinal-cord-native-core/v2'
bank=Path(args.bank).resolve();m=json.loads((bank/'tissue.json').read_text());assert m['passed']
context=bank/'fixed-context-screen.json';screen=json.loads(context.read_text());assert screen['passed']and screen['actual_packed_float32']and screen['source_cell_sha256']==m['bank_sha256']
metadata=json.loads((base/'female-31-segments-v2/cells.json').read_text());oldmeta=json.loads((base/'female-cage-rounded-v3/cells.json').read_text())
old_ids=set(old['partIds']+old['coveringPartIds']);oldparts={p['id']:p for p in a['parts']if p['id']in old_ids}
donors={p['id']for p in a['parts']if p['id'].startswith('HRA:VH_F_')and'spinal_cord'in p['id'].lower()and not p.get('suppressed')}
retired=old_ids|donors
for p in a['parts']:
    if p['id']in retired:p.update(suppressed=True,defaultHidden=True)
a['defaultHidden']=list(dict.fromkeys([*a.get('defaultHidden',[]),*sorted(retired)]))
chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=m['total_bytes']))
concepts={c['id']:c for c in a['concepts']};core=[];covers=[];macros={};regions={};attributes=metadata['attributes']
for name,row in m['parts'].items():
    info=attributes[name];native=info['native_region'];segment=info['axial_region'];family=native.split('-')[0];pid=row['id'];cid='TA98CORD:female:'+native
    template=next(p for id,p in oldparts.items()if id.removeprefix('TA98CORD:female:').split('@')[0]==native)
    part=dict(template,**row,chunk=chunk,name=template['name'].split(' · ')[0]+' · '+segment,sourceId='spinal_cord:'+name,
        conceptId=cid,reviewStatus='Finite registered segment; root territories and microscopic boundaries remain approximate.')
    part.pop('suppressed',None);part.pop('defaultHidden',None);a['parts'].append(part)
    (covers if family=='pia'else core).append(pid);macros.setdefault(cid,[]).append(pid)
    if family!='pia':regions.setdefault(segment,[]).append(pid)
for cid,ids in macros.items():concepts[cid].update(elements=ids,componentElements=ids)
def replaces(new_name,prior):
    native=prior.removeprefix('TA98CORD:female:').split('@')[0];info=attributes[new_name]
    if info['native_region']!=native:return False
    if '@'not in prior:return True
    low,high=oldmeta['attributes'][prior.removeprefix('TA98CORD:female:')]['y_range_m'];middle=sum(info['y_range_m'])/2
    return low<=middle<=high
groups=[]
for g in old['termGroups']:
    ids=[row['id']for name,row in m['parts'].items()if any(replaces(name,prior)for prior in g['partIds'])]
    assert ids,g['term'];concepts[g['conceptId']].update(elements=ids,componentElements=ids)
    groups.append(dict(g,partIds=ids))
segments=[]
for segment,low,high in metadata['bands'][::-1]:
    code='018'if segment=='Co1'else{'C':'014','T':'015','L':'016','S':'017'}[segment[0]];term='A14.1.02.'+code
    parent=next(g['conceptId']for g in groups if g['term']==term);cid='TA98CORD:female:segment:'+segment;ids=regions[segment]
    concepts[cid]=dict(id=cid,name=segment+' spinal cord segment',elements=ids,componentElements=ids,partOf=parent,
        coverageLimitation='Supplemental segment control. Its transverse boundaries approximate overlapping root territories; caudal levels use reference relative lengths. This is not an additional TA98 term.')
    segments.append(dict(name=segment,conceptId=cid,parentTerm=term,partIds=ids,coveringPartIds=['TA98CORD:female:pia@'+segment],yRange=[low,high]))
concepts[old['wholeConceptId']].update(elements=core,componentElements=core);a['concepts']=list(concepts.values())
marker=dict(old,schema='spinal-cord-native-core/v3',partIds=core,coveringPartIds=covers,bankSha256=m['bank_sha256'],termGroups=groups,
    retiredPartIds=sorted(set(old.get('retiredPartIds',[]))|retired),regionalGroups=regions,segments=segments,segmentsAreAdditionalTA98Terms=False)
marker['courseFit']=dict(old['courseFit'],bank_sha256=m['bank_sha256'],segment_source_sha256=digest_file(base/'female-31-segments-v2/cells.npz'))
a['finiteSpinalIntegration']=marker;a['triangles']=sum(p['indexCount']//3 for p in a['parts']if not p.get('suppressed'));a['partsCount']=len(a['parts'])
atomic_json(out/'atlas.json',a)
plan=json.loads((base/'female-source-v1/inventory.json').read_text())
inventory=dict(schema='ta98-geometry-inventory/v1',model='local-female-ta98',model_sex='female',root_term='A14.1.02.001',source_catalogue=str(source),source_catalogue_sha256=digest_file(source),
    authority_path=str(ROOT/plan['authority']['path']),authority_sha256=plan['authority']['sha256'],terms=[dict(r,source_meshes=r['source_meshes']['female'],applicable_to_model=True,
        candidate_part_ids=next((g['partIds']for g in groups if g['term']==r['authority']['id']),[]))for r in plan['terms']],term_count=len(plan['terms']),anatomical_acceptance=False)
atomic_json(out/'inventory.json',inventory);rows={p['id']:p for p in m['parts'].values()};atomic_json(out/'audited-layout.json',dict(m,parts=rows,origin=[0.,1.35,-.05]))
interfaces=[]
for pid,p in rows.items():
    for g in p['tissueDrawGroups']:
        peer=g.get('neighbor')
        if not peer or pid>=peer:continue
        interfaces.append(dict(kind='shared-interface',first=pid,second=peer,parameters=dict(first_faces=list(range(g['start']//3,(g['start']+g['count'])//3)),
            second_faces=[i for gg in rows[peer]['tissueDrawGroups']if gg.get('neighbor')==pid for i in range(gg['start']//3,(gg['start']+gg['count'])//3)])))
spec=dict(model='local-female-ta98',proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(source),changed_parts=[*sorted(retired),*core,*covers],closed_parts=[*core,*covers],inventory=str(out/'inventory.json'),
    max_pairs=2000000,max_tests=6000000,tissue_audit_layouts=[str(out/'audited-layout.json')],tissue_audit_checkpoints=str(base/'female-31-segment-bank-v1-checks'),interface_checks=interfaces,
    assembly=dict(whole_organ=old['wholeConceptId'],part_ids=[*core,*covers],fixed=['Entire native skeleton','All unrelated anatomy','Entire male model'],
        adjoining_interfaces=[dict(kind='shared-interface',scope='Every radial and axial interface is stored once with reciprocal ownership'),dict(kind='open-space',scope='The central canal remains open, with annular pial caps')],
        remaining_limits=marker['remaining']),evidence_files=[str(bank/'tissue.json'),str(base/'female-31-segments-v2/cells.json'),str(context)])
atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec))
print(json.dumps(dict(segments=len(segments),finite_parts=len(rows),bytes=m['total_bytes'],interfaces=len(interfaces),cord_terms=marker['representedCordTerms'],hidden_donor_segments=len(donors))),flush=True)
