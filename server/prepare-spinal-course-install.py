"""Prepare the passed cord/pia placement correction without adding term claims.

This revision installs the 97 fully audited existing tissues only. The new
dura/arachnoid and unresolved root ports remain separate uninstalled work.
"""
from pathlib import Path
import argparse,copy,json
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
before=json.loads(source.read_text());a=copy.deepcopy(before);marker=a['finiteSpinalIntegration'];assert marker['schema']=='spinal-cord-native-core/v2'
bank=base/'female-cage-cord-bank-v1';m=json.loads((bank/'tissue.json').read_text());assert m['passed']
context=base/'female-cage-rounded-v3/fixed-context-screen.json';assert json.loads(context.read_text())['passed']
ids=[*marker['partIds'],*marker['coveringPartIds']];rows={p['id']:p for p in m['parts'].values()};assert set(ids)==set(rows)
chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=m['total_bytes']))
for part in a['parts']:
    if part['id']in rows:
        part.update(rows[part['id']],chunk=chunk,reviewStatus='Shared positive-cage cord/pia fit; regional anatomical boundaries and nerve-root interfaces remain provisional.')
marker.update(bankSha256=m['bank_sha256'],courseFit=dict(bank_sha256=m['bank_sha256'],source=digest_file(base/'female-common-course-v7/course.npz'),
    fixed_skeleton_screen_sha256=digest_file(context),source_fit='Positive piecewise-affine cage with conformed tissue charts and bounded Float32 repair',
    core_term_count_unchanged=True),anatomicallyAccepted=False)
marker['remaining']=['95 cord terms remain unrepresented','Dura/arachnoid root ports and sleeves, filum and denticulate attachments','Verify native root/segment and cranial interfaces','Named nuclei, laminae and individual tract territories']
a['triangles']=sum(p['indexCount']//3 for p in a['parts']if not p.get('suppressed'));a['partsCount']=len(a['parts'])
atomic_json(out/'atlas.json',a)
plan=json.loads((base/'female-source-v1/inventory.json').read_text());groups=marker['termGroups']
inventory=dict(schema='ta98-geometry-inventory/v1',model='local-female-ta98',model_sex='female',root_term='A14.1.02.001',source_catalogue=str(source),source_catalogue_sha256=digest_file(source),
    authority_path=str(ROOT/plan['authority']['path']),authority_sha256=plan['authority']['sha256'],terms=[dict(r,source_meshes=r['source_meshes']['female'],applicable_to_model=True,
        candidate_part_ids=next((g['partIds']for g in groups if g['term']==r['authority']['id']),[]))for r in plan['terms']],term_count=len(plan['terms']),anatomical_acceptance=False)
atomic_json(out/'inventory.json',inventory);layout=dict(m,parts=rows,origin=[0.,1.35,-.05]);atomic_json(out/'audited-layout.json',layout)
interfaces=[]
for pid,p in rows.items():
    for g in p['tissueDrawGroups']:
        peer=g.get('neighbor')
        if not peer or pid>=peer:continue
        interfaces.append(dict(kind='shared-interface',first=pid,second=peer,parameters=dict(first_faces=list(range(g['start']//3,(g['start']+g['count'])//3)),
            second_faces=[i for gg in rows[peer]['tissueDrawGroups']if gg.get('neighbor')==pid for i in range(gg['start']//3,(gg['start']+gg['count'])//3)])))
spec=dict(model='local-female-ta98',proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(source),changed_parts=ids,closed_parts=ids,inventory=str(out/'inventory.json'),
    max_pairs=2000000,max_tests=6000000,tissue_audit_layouts=[str(out/'audited-layout.json'),str(base/'female-cage-audit-reference-v1/audit-layout.json')],tissue_audit_checkpoints=str(bank.with_name(bank.name+'-checks')),interface_checks=interfaces,
    assembly=dict(whole_organ=marker['wholeConceptId'],part_ids=ids,fixed=['Entire skeleton','Every unrelated organ, vessel and peripheral nerve','Entire male model'],
        adjoining_interfaces=[dict(kind='shared-interface',scope='Every explicitly owned grey/white and pial interface'),dict(kind='open-space',scope='Central canal remains a void; cranial pial continuation remains open'),
            dict(kind='unfinished-native-attachment',scope='Root sleeves and verified cranial/root interfaces remain separate unfinished work')],remaining_limits=marker['remaining']),
    evidence_files=[str(bank/'tissue.json'),str(context),str(base/'female-cage-rounded-v3/rounding-repair.json'),str(base/'female-common-course-v7/course.json')])
atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec));print(json.dumps(dict(existing_finite_parts=len(ids),interfaces=len(interfaces),bytes=m['total_bytes'],cord_terms=marker['representedCordTerms'],male_unchanged=True)),flush=True)
