"""Prepare a source-bound open-reference partition job; never apply it."""
from pathlib import Path
import hashlib
import json
import sqlite3
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.ta98_authority import corrected_terms,term_sex

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-native-venous-regions-v1/literal-forearm-pack-v1'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'
CANDIDATE=BASE/'proposed-atlas.json'
PROPOSAL='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
FEMALE='8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'


def record_hash(records):
    return hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def main():
    began=time.monotonic();assert digest_file(SOURCE)==EXPECTED and digest_file(CANDIDATE)==PROPOSAL
    old=json.loads(SOURCE.read_bytes());new=json.loads(CANDIDATE.read_bytes());bp={p['id']:p for p in old['parts']};ap={p['id']:p for p in new['parts']}
    marker=new['nativeVenousRegionIntegration'];assert len(marker['groups'])==4
    changed=sorted(pid for pid in bp.keys()|ap.keys() if bp.get(pid)!=ap.get(pid));assert len(changed)==8
    original_ids={g['wholeConceptId'] for g in marker['groups']};new_ids={g['regionPartId'] for g in marker['groups']}
    assert set(changed)==original_ids|new_ids
    skeleton={pid:p for pid,p in bp.items() if p.get('system')=='skeletal'};cardiac={pid:p for pid,p in bp.items() if p.get('system')=='cardiac'}
    active_cardiac={pid:p for pid,p in cardiac.items() if not p.get('suppressed')}
    heart_path=ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/preview-proof.json'
    heart_proof=json.loads(heart_path.read_bytes());heart_ids={row['part'] for row in heart_proof['physical_tissues']};assert len(heart_ids)==33
    heart_ids.update(pid for name,values in old['finiteCardiacIntegration'].items() if name.endswith('PartIds') for pid in values)
    heart={pid:bp[pid] for pid in heart_ids};assert len(skeleton)==1044 and len(heart)==41 and len(cardiac)==48 and len(active_cardiac)==44
    assert all(ap[pid]==p for pid,p in bp.items() if pid not in original_ids)
    assert new['chunks'][:len(old['chunks'])]==old['chunks']
    female=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';assert digest_file(female)==FEMALE
    protected_files=[dict(path=str(female),sha256=FEMALE),dict(path=str(heart_path),sha256=digest_file(heart_path))]
    for sex,sha in [('male','375bed821531d0bb6c85f5b330362161c9e50dcef2af40f670292844a03e3040'),('female','f88dd44e18df5fe343d3a9724a68ef4d76d17803b7f0b388f9eb26d5772bded5')]:
        path=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/cardiac-label-sites.json';assert digest_file(path)==sha
        assert json.loads(path.read_bytes())['integration']['labeledTermCount']==113
        protected_files.append(dict(path=str(path),sha256=sha))
    reader=AtlasReader()
    for chunk in sorted({p['chunk'] for p in list(skeleton.values())+list(cardiac.values())+list(heart.values())}):
        path=reader.asset(old['chunks'][chunk]['url']);protected_files.append(dict(path=str(path),sha256=digest_file(path)))
    protection=dict(schema='native-forearm-protected-source-bindings/v1',source_catalogue_sha256=EXPECTED,proposal_sha256=PROPOSAL,
        unchanged_skeletal_records=len(skeleton),unchanged_skeletal_records_sha256=record_hash(skeleton),unchanged_active_cardiac_records=len(active_cardiac),
        unchanged_installed_heart_assembly_records=len(heart),unchanged_installed_heart_assembly_sha256=record_hash(heart),
        unchanged_all_cardiac_records=len(cardiac),unchanged_all_cardiac_records_sha256=record_hash(cardiac),
        both_cardiac_label_terms=113,protected_files=protected_files,all_existing_chunks_literal=True,all_unrelated_records_literal=True,
        female_catalogue_sha256=FEMALE,installed=False)
    protected=BASE/'protected-source-bindings.json';atomic_json(protected,protection)
    snapshot=BASE/'source-atlas.json';atomic_bytes(snapshot,SOURCE.read_bytes());assert digest_file(snapshot)==EXPECTED
    authority=ROOT/'.local-models/ta98-authority.sqlite';correction=ROOT/'geometry_tools/ta98_authority.py'
    con=sqlite3.connect(authority);con.row_factory=sqlite3.Row;terms=corrected_terms({r['id']:dict(r) for r in con.execute('SELECT * FROM ta98')});con.close()
    def descendants(root):
        selected={root};previous=-1
        while previous!=len(selected):
            previous=len(selected);selected|={t for t,r in terms.items() if r.get('parent_id') in selected}
        return selected
    full=descendants('A12.3.08.015')|descendants('A12.3.08.018');codes=full|{'A12.3.08.014','A12.3.08.021','A12.3.08.022'}
    contract=json.loads((ROOT/'reports/ta98-male/native-forearm-venous-region-plan-v1.json').read_bytes());retrieval={r['id']:r for r in contract['terms']}
    assert codes==set(retrieval) and all(terms[t]['parent_id']==retrieval[t]['parent_id'] for t in codes)
    source_plan=json.loads((BASE.parent/'source-plan.json').read_bytes());original_meshes={pid:source_plan['sources'][pid] for pid in original_ids}
    scope='Literal native outer-surface references partitioned at provisional radial-head/distal-radius ring stations. Four new side-qualified forearm controls provide two additional direct TA98 surface references. No coordinates, source triangles, normals or UVs change. No finite wall, blood-space patency, valve or repaired junction is claimed.'
    rows=[]
    for term in sorted(codes):
        applicable=term_sex(terms,term)!='female';active=retrieval[term].get('active_exact_source_part_ids',[])
        if term in ['A12.3.08.021','A12.3.08.022']:
            ids=[g['regionPartId'] for g in marker['groups'] if g['regionTerm']==term]
            source_ids=[g['wholeConceptId'] for g in marker['groups'] if g['regionTerm']==term]
        else:ids=active;source_ids=active
        sources=[]
        for pid in source_ids:
            p=bp[pid]
            m=reader.read_part(old,p,frame='model:local-ta98:assembled-metres',catalogue_path=SOURCE,catalogue_sha=EXPECTED)
            sources.append(dict(part_id=pid,source_catalogue=str(snapshot),geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,
                role='Unchanged contextual source reference; not copied, corrected or newly accepted' if pid not in original_ids else 'Literal licensed native whole outer-surface source before disjoint forearm partition',license=p.get('license')))
        rows.append(dict(authority=terms[term],applicable_to_model=applicable,representation_class='anatomical-collection' if term=='A12.3.08.014' else 'macroscopic-vessel',
            source_meshes=sources,missing_distinct_representation=not active,distinct_representation_verified=False,
            planned_representation=dict(part_ids=ids,representation_class='outer-surface-reference' if term in ['A12.3.08.015','A12.3.08.018','A12.3.08.021','A12.3.08.022'] else 'unchanged-context-or-missing',
                anatomical_scope=scope if term in ['A12.3.08.021','A12.3.08.022'] else 'Keep current representation or missing status unchanged; no parent alias substitutes for this distinct term.')))
    inventory=dict(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',root_term='A12.3.08.014',whole_source_roots=['A12.3.08.015','A12.3.08.018'],
        complete_applicable_whole_descendants=sorted(full),additional_explicit_regional_terms=['A12.3.08.021','A12.3.08.022'],context_root='A12.3.08.014',
        authority_path=str(authority),authority_sha256=digest_file(authority),authority_corrections_path=str(correction),authority_corrections_sha256=digest_file(correction),
        source_catalogue=str(snapshot),source_catalogue_sha256=EXPECTED,term_count=len(rows),applicable_term_count=sum(r['applicable_to_model'] for r in rows),
        checked_in_report=str(ROOT/'reports/ta98-male/native-forearm-venous-region-plan-v1.json'),automatic_mapping_is_not_representation_acceptance=True,anatomical_acceptance=False,terms=rows)
    invpath=BASE/'source-bound-inventory.json';atomic_json(invpath,inventory)
    interfaces=[]
    for group in marker['groups']:
        for interface in group['sharedInterfaces']:
            interfaces.append(dict(kind='literal-common-source-ring-not-a-new-finite-tissue-seam',first=group['remainderPartId'],second=group['regionPartId'],
                source_ring_station=interface['source_ring_station'],source_vertex_ids=interface['source_vertex_ids'],
                position_words_sha256=interface['position_words_sha256'],normal_words_sha256=interface['normal_words_sha256'],added_cap_triangles=0,
                blood_lumen_communication_unverified=True))
        interfaces.append(dict(kind='original-upstream-downstream-source-boundaries-retained',whole_concept=group['wholeConceptId'],
            complete_original_source_face_count=group['sourceFaceCount'],no_new_coordinate_or_triangle=True,no_new_patency_claim=True))
    assembly=dict(whole_organ=dict(name='Four native male cephalic/basilic source-surface assemblies',authority_roots=['A12.3.08.015','A12.3.08.018'],part_ids=changed,scope=scope),
        adjoining_interfaces=interfaces,fixed_skeleton=True,fixed_heart=True,fixed_female=True,fixed_all_unrelated_geometry=True,source_outer_surfaces_only=True,anatomical_acceptance=False)
    evidence=dict(checked_in_inventory=str(ROOT/'reports/ta98-male/native-forearm-venous-region-plan-v1.json'),source_bound_inventory=str(invpath),
        literal_partition_and_corner_words=str(BASE/'literal-partition-proof.json'),source_plan=str(BASE.parent/'source-plan.json'),actual_bank=str(BASE/'geometry.bin'),
        controls_and_no_duplicate_leaf_check=str(BASE/'controls-proof.json'),final_readiness=str(BASE/'readiness.json'),protected_source_bindings=str(protected),
        protected_source_files=[r['path'] for r in protected_files],source_snapshot=str(snapshot),native_review=str(BASE/'root-render-review-v1.json'),
        app_control_sources=[str(ROOT/'app/native-forearm-venous-controls.ts'),str(ROOT/'app/system-tree.tsx')],
        boundary_ancestry=[str(BASE/(g['regionTerm']+'-'+g['side']+'-ancestry.npz')) for g in marker['groups']])
    spec=dict(model='local-ta98',expected_catalogue_sha256=EXPECTED,proposed_catalogue=str(CANDIDATE),changed_parts=changed,
        require_closed=False,closed_parts=[],inventory=str(invpath),assembly=assembly,interface_checks=[],evidence=evidence,
        remaining_limits=[scope,'Thoraco-acromial vein A12.3.08.016 remains missing; contextual accessory-cephalic sources are not newly repaired or accepted.',
            'The source surfaces have no inner tissue wall; original upstream/downstream blood connectivity is unverified.',
            'Right diagnostic section removes most source surfaces and is not a lumen proof.'],max_pairs=250000,max_tests=2000000)
    atomic_json(BASE/'install-spec.json',spec);atomic_json(BASE/'catalogue-candidate-job.json',dict(kind='catalogue-candidate',parameters=spec))
    assert digest_file(SOURCE)==EXPECTED and digest_file(female)==FEMALE
    report=dict(candidate_sha256=PROPOSAL,source_sha256=EXPECTED,changed_parts=changed,new_direct_terms=2,applicable_inventory_terms=len(rows),whole_descendants=sorted(full),
        protected_skeleton=1044,protected_installed_heart_assembly=41,protected_active_cardiac=44,protected_total_cardiac=48,protected_banks=len(protected_files)-4,
        native_review_accepted_local_preview_only=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(BASE/'job-preparation-receipt.json',report);print(json.dumps(report))


if __name__=='__main__':main()
