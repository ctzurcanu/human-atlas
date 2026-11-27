"""SMA-adapter pattern: source-bound preparation only, never installation."""
from pathlib import Path
import copy,json,sqlite3,sys,time,argparse
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh

args=argparse.ArgumentParser();args.add_argument('--side',choices=['left','right'],default='left');args.add_argument('--version',default='v5');args.add_argument('--layout',default='owned-bank-v2');args=args.parse_args();side=args.side;root_name='external-carotid-'+side
T=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/('male-'+side+'-lingual-finite-network-'+args.version);O=B/('male-'+side+'-lingual-install-proposal-v1');CURRENT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='9cf8b30456620f7b7bd6f52311ba818fc94f935beee24e6e07cb953d14106f64'if side=='left'else'978f5a0a25d505f1fcfcca9d2b4cd9fd07023ae22e1363f15131dcd7d0e750b1';assert digest_file(CURRENT)==EXPECTED;O.mkdir(exist_ok=False)
inputs=json.loads((B/'male-lingual-source-plan-v1/inputs.json').read_text());layout_path=D/args.layout/'tissue.json';layout=json.loads(layout_path.read_text());assert layout['passed']and layout['checks']['passed'];assert json.loads((D/'nearby-context-v1/evidence.json').read_text())['passed'];assert json.loads((D/'independent-check-v1/evidence.json').read_text())['numerical_route_gate_pass']
reader=AtlasReader();before=json.loads(CURRENT.read_bytes());after=copy.deepcopy(before);parts={p['id']:p for p in before['parts']};native_id='ZA:External carotid artery.'+('l'if side=='left'else'r');source_record=next(r for r in inputs['fixed_sources']if r['role']==root_name);native=reader.read_part(before,parts[native_id],frame=layout['frame'],catalogue_path=CURRENT,catalogue_sha=EXPECTED);assert parts[native_id]==source_record['source_record'];assert native.geometry_sha==source_record['geometry_sha256'];assert native.source['bank_sha256']==source_record['bank_sha256']
bank=Path(layout['bank']);assert digest_file(bank)==layout['bank_sha256'];chunk=len(after['chunks']);after['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size));owned=json.loads((D/'owned-bank-stage.json').read_text());term_by_name=dict(zip(owned['names'],owned['term_by_component']));term_names={r['id']:r['name_en']for r in inputs['terms']};parent_concept=parts[native_id].get('conceptId',native_id);main_id=layout['part_ids']['lingual'];newids=[pid for name,pid in layout['part_ids'].items()if name!=root_name];concept_groups={};changed=[]
limitation=side.capitalize()+' macroscopic source-calibre arterial family only. Native ECA source positions and far outer UV/normal corners are retained within declared Float32 encoding error. New lingual artery is genuinely open to the ECA; five terminal annuli enter tested source recipients. Wall thickness and small-branch taper are illustrative. Microscopic perfusion beds and pre-existing ECA branch mouths remain unfinished. The shared ostial support includes the first 8.35 mm of the new branch in the parent finite cell.'
for name,geometry in layout['parts'].items():
    p=copy.deepcopy(parts[native_id]);pid=geometry['id']
    for key in ['geometryScale','geometryInstance','sourceOffset','normalSigns','suppressed','hidden']:p.pop(key,None)
    p.update(copy.deepcopy(geometry),chunk=chunk,pairedTissue=True,sectionSolid=True,sectionAssembly=side+'-lingual-arterial-family:male',lingualArterialWallRevision='male-'+side+'-lingual-finite-network-'+args.version,coverageLimitation=limitation,uvMethod='Literal native ECA outer UV/normal words on retained source corners, barycentric at clipped source faces; generated walls use local functional material.')
    if name==root_name:after['parts'][after['parts'].index(parts[native_id])]=p
    else:
        term=term_by_name[name];cid='TA98ARTERIAL:male:dorsal-lingual:'+side+':whole'if term=='A12.2.05.017'else pid
        p.update(name=term_names[term].capitalize()+' ('+side+')'+(' · '+str(int(name[-1])+1)if name.startswith('dorsal-lingual-')else''),conceptId=cid,sourceId='independent-lingual-course/native-ECA-and-recipient-constrained',ta98Term=term,ta98Kind='exact')
        p['groups']=list(p.get('groups',[]))+['External carotid artery','Lingual artery'];after['parts'].append(p);concept_groups.setdefault(cid,dict(term=term,name=term_names[term].capitalize()+' ('+side+')',elements=[]))['elements'].append(pid)
    changed.append(pid)
for cid,item in concept_groups.items():
    after['concepts'].append(dict(id=cid,name=item['name'],elements=item['elements'],partOf=parent_concept if item['term']=='A12.2.05.015'else main_id,ta98Term=item['term'],ta98Kind='exact',coverageLimitation=limitation))
ap={p['id']:p for p in after['parts']};assert sorted(k for k in parts.keys()|ap.keys()if parts.get(k)!=ap.get(k))==sorted(changed);assert after['chunks'][:len(before['chunks'])]==before['chunks'];assert all(ap[p['id']]==p for p in before['parts']if p['id']not in changed)
proposal=O/'proposed-atlas.json';atomic_json(proposal,after)
# Independent literal source-corner replay against the actual packed binary.
emitted=reader.read_part(after,ap[native_id],frame=layout['frame'],catalogue_path=proposal);count=ap[native_id]['tissueDrawGroups'][0]['count']//3;tree=cKDTree(native.vertices);source_faces={tuple(sorted(map(int,f))):i for i,f in enumerate(native.faces)};matches=[];maxerror=0.
for fi,triangle in enumerate(emitted.triangles[:count]):
    distances,nearest=tree.query(triangle)
    if distances.max()>5e-9:continue
    original=source_faces.get(tuple(sorted(map(int,nearest))))
    if original is None:continue
    out=emitted.faces[fi]
    assert np.array_equal(emitted.attributes['uvs'][out],native.attributes['uvs'][nearest]);assert np.array_equal(emitted.attributes['normals'][out]*emitted.attributes['normalSigns'][out,None],native.attributes['normals'][nearest]);matches.append(original);maxerror=max(maxerror,float(distances.max()))
assert len(matches)==len(set(matches))and len(matches)>len(native.faces)-30
authority=ROOT/'.local-models/ta98-authority.sqlite';corrections=ROOT/'geometry_tools/ta98_authority.py';con=sqlite3.connect(authority);con.row_factory=sqlite3.Row;terms=[dict(r)for r in con.execute('WITH RECURSIVE d(id)AS(SELECT ? UNION ALL SELECT ta98.id FROM ta98 JOIN d ON ta98.parent_id=d.id)SELECT ta98.* FROM ta98 JOIN d USING(id)ORDER BY ta98.id',('A12.2.05.015',))];assert len(terms)==5
inventory=dict(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',root_term='A12.2.05.015',authority_path=str(authority),authority_sha256=digest_file(authority),authority_corrections_path=str(corrections),authority_corrections_sha256=digest_file(corrections),source_catalogue=str(CURRENT),source_catalogue_sha256=EXPECTED,term_count=5,applicable_term_count=5,checked_in_report=str(ROOT/'TA98-CIRCULATION.md'),automatic_mapping_is_not_representation_acceptance=True,anatomical_acceptance=False,terms=[])
for t in terms:inventory['terms'].append(dict(authority=t,applicable_to_model=True,representation_class='tissue-region',source_meshes=[dict(part_id=native_id,role='profile and parent mouth only, not alias of the new branch',geometry_sha256=native.geometry_sha)],missing_distinct_representation=False,distinct_representation_verified=False,planned_representation=dict(part_ids=[pid for name,pid in layout['part_ids'].items()if term_by_name[name]==t['id']],representation_class='finite-vessel-wall',anatomical_scope=limitation)))
inv=O/'source-bound-inventory.json';atomic_json(inv,inventory);interfaces=[]
for interface in owned['interfaces']:
    first=layout['part_ids'][interface['parent']];second=layout['part_ids'][interface['child']];params={}
    for role,pid,other in [('first',first,second),('second',second,first)]:
        params[role+'_faces']=[i for g in ap[pid]['tissueDrawGroups']if g.get('neighbor')==other for i in range(g['start']//3,(g['start']+g['count'])//3)];assert len(params[role+'_faces'])==24
    interfaces.append(dict(first=first,second=second,kind='shared-interface',parameters=params))
proof=dict(schema='left-lingual-source-bound-install-proposal/v1',source_catalogue_sha256=EXPECTED,proposed_catalogue_sha256=digest_file(proposal),changed_parts=changed,all_unrelated_records_and_banks_literal=True,native_parent_source_record_bound=True,native_outer_faces_with_literal_UV_normal_corners=len(matches),native_outer_faces=len(native.faces),maximum_native_position_encoding_error_m=maxerror,actual_new_arterial_owners=6,actual_unique_TA98_terms=5,zero_new_nearby_native_source_crossings=True,remaining_limits=[limitation],anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-T);atomic_json(O/'source-rebase-proof.json',proof)
encoding=D/'position-encoding-v1.json'
if not encoding.exists():
    encoding=O/'position-encoding-proof.json';atomic_json(encoding,dict(schema='direct-source-control-float32-serialization/v1',bank_sha256=layout['bank_sha256'],selected_origin_m=layout['parts']['lingual']['geometryOrigin'],maximum_position_encoding_error_m=layout['maximum_position_encoding_error_m'],all_actual_emitted_cells_and_interfaces_checked=True,passed=layout['passed'],anatomical_control_change=False))
evidence=dict(source_bound_inventory=str(inv),source_rebase_proof=str(O/'source-rebase-proof.json'),checked_in_inventory=str(ROOT/'TA98-CIRCULATION.md'),finite_network=str(D/'evidence.json'),owned_bank=str(layout_path),ownership=str(D/'owned-bank-stage.json'),patency_and_recipients=str(D/'independent-check-v1/evidence.json'),nearby_native_context=str(D/'nearby-context-v1/evidence.json'),float32_position_encoding=str(encoding),fixed_source_plan=str(B/'male-lingual-source-plan-v1/inputs.json'))
for name in ['independent-check-emitted-v1','nearby-context-emitted-v1']:
    path=D/name/'evidence.json'
    if path.exists():evidence[name]=str(path)
assembly=dict(whole_organ=dict(name='Native '+side+' ECA and five-term lingual arterial family',authority_root='A12.2.05.015',part_ids=changed,scope=limitation),adjoining_interfaces=owned['interfaces']+[dict(kind='five-open-terminal-mouths-inside-source-recipients',evidence=str(D/'independent-check-v1/evidence.json'),microscopic_perfusions_unmodeled=True),dict(kind='pre-existing-ECA-branch-mouths-pending',native_outer_courses_fixed=True)],fixed_skeleton=True,fixed_heart=True,fixed_all_unrelated_geometry=True,anatomical_acceptance=False)
parameters=dict(model='local-ta98',expected_catalogue_sha256=EXPECTED,proposed_catalogue=str(proposal),changed_parts=changed,closed_parts=changed,require_closed=True,inventory=str(inv),assembly=assembly,interface_checks=interfaces,evidence=evidence,remaining_limits=[limitation],tissue_audit_checkpoints=str(D/'owner-validation-checkpoints-v1'),tissue_audit_layouts=[str(layout_path)],max_pairs=250000,max_tests=2000000)
atomic_json(O/'catalogue-candidate-job.json',dict(kind='catalogue-candidate',parameters=parameters));assert digest_file(CURRENT)==EXPECTED;print(json.dumps(dict(proposed=str(proposal),sha256=digest_file(proposal),job=str(O/'catalogue-candidate-job.json'),changed_parts=len(changed),new_terms=5,bank_bytes=layout['total_bytes'],native_faces_with_literal_UV_normal_corners=len(matches),elapsed_seconds=time.monotonic()-T)))
