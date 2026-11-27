"""Prepare the standard guarded job for the reviewed14-control female core."""
import copy,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'server')]
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from tongue_basal_contract import require
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
CANDIDATE=W/'female-core-only-emitted-v1/candidate-atlas.json'
FEMALE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
FEMALE_SHA='770f47f802c79cf6f3e4336bf2fa7c42cf710998069221f533a1ea3f8322e470'
OUT=W/'female-core-only-install-job-v1';FRAME='model:local-female-ta98:assembled-metres'
def read(p):return json.loads(Path(p).read_text())
def main():
    started=time.monotonic();require(not OUT.exists() and digest_file(FEMALE)==FEMALE_SHA,'Fresh job stage and unchanged female source required')
    after,before=read(CANDIDATE),read(FEMALE);packet_path=CANDIDATE.with_name('transfer-evidence.json');packet=read(packet_path)
    require(packet['candidate']['sha256']==digest_file(CANDIDATE) and packet['all_six_decoded_arrays_byte_identical'] is True,'Final emitted packet changed')
    marker=after['finiteTongueIntegration'];active=set(marker['partIds'])
    require(marker.get('partialMode')=='female-core-only' and len(active)==6 and len(marker['termGroups'])==14
            and marker.get('representedAdultControls')==14 and marker.get('missingAdultTerms')==21
            and marker.get('anatomicalAcceptance') is False and marker.get('tissueOwnershipClaim') is False,'Exact partial scope required')
    ap,bp=({p['id']:p for p in a['parts']} for a in (after,before));suppressed=set(packet['suppressed_source_parts'])
    require(set(ap)-set(bp)==active and not(set(bp)-set(ap)) and len(suppressed)==7,'Unexpected added/removed/suppressed anatomy')
    for pid,p in bp.items():require(ap[pid]==(dict(p,suppressed=True) if pid in suppressed else p),'Unrelated anatomy or skeleton changed: '+pid)
    require(all(after['materials'][k]==v for k,v in before['materials'].items()) and after['chunks'][:len(before['chunks'])]==before['chunks'],'Original materials/banks changed')
    aliases=set(packet['legacy_alias_elements']);bc={c['id']:c for c in before['concepts']};ac={c['id']:c for c in after['concepts']}
    require(all(ac[k]==v for k,v in bc.items() if k not in aliases),'An unrelated/native extrinsic concept changed')
    require(set(after['defaultHidden'])==set(before['defaultHidden'])|suppressed,'Unrelated visibility changed')
    audit_path=W/'female-core-only-independent-audit-v1/evidence.json';audit=read(audit_path)
    require(audit['candidate_sha256']==digest_file(CANDIDATE) and audit['packet_sha256']==digest_file(packet_path)
            and audit['passed'] is True and audit['actual_float32_banks_independently_read'] is True
            and audit['all_actual_jobs_complete'] is True and audit['self_defect_pairs']==0
            and audit['nonpositive_corner_normals']==0 and audit['proper_pair_crossings']==0
            and len(audit['all_fifteen_pair_constraints'])==15,'Actual six-part independent gates must pass')
    reader=AtlasReader();meshes={pid:reader.read_part(after,ap[pid],frame=FRAME) for pid in active};bindings={str(CANDIDATE):digest_file(CANDIDATE),str(FEMALE):FEMALE_SHA,str(packet_path):digest_file(packet_path),str(audit_path):digest_file(audit_path)}
    for row in packet['actual_world_meshes']:
        m=meshes[row['part_id']];require(digest_file(row['mesh']['path'])==row['mesh']['sha256'] and m.geometry_sha==row['geometry_sha256']
            and m.attributes_sha==Mesh.load(row['mesh']['path']).attributes_sha,'Final actual geometry/attributes changed')
        bindings[row['mesh']['path']]=row['mesh']['sha256']
    target_rows=[]
    for stage,count in [('bones',30),('soft',36)]:
        path=W/('female-partial-actual-'+stage+'-v1/evidence.json');d=read(path)
        require(d['passed'] is True and d['all_pairs_complete'] is True and d['proper_crossing_pairs']==0
                and d['coplanar_overlap_pairs']==0 and len(d['pairs'])==count,'Original literal-array target audit failed/incomplete')
        bindings[str(path)]=digest_file(path)
        for row in d['pairs']:
            require(row['partId'] in active and row['first_geometry_sha256']==meshes[row['partId']].geometry_sha,'Target audit binds another new mesh')
            target=reader.read_part(after,ap[row['targetPartId']],frame=FRAME);require(target.geometry_sha==row['second_geometry_sha256'],'Fixed target changed')
            meshes[row['targetPartId']]=target;target_rows.append(row)
            require(digest_file(row['witnesses'])==row['witnesses_sha256'],'Target witness changed');bindings[row['witnesses']]=row['witnesses_sha256']
    expected={(pid,r['part_id']) for pid in active for r in packet['fixed_context_meshes']}
    require({(r['partId'],r['targetPartId']) for r in target_rows}==expected and len(expected)==66,'All66 fixed target pairs required')
    fixed={second for _,second in expected};constraints=[dict(kind='fixed',first=pid,parameters=dict(geometry_sha256=meshes[pid].geometry_sha)) for pid in sorted(fixed)]
    topo={pid:topology(m) for pid,m in meshes.items()}
    for first,second in sorted(expected):
        params=dict(occupancy='solid') if topo[first]['closed_oriented_manifold'] and topo[second]['closed_oriented_manifold'] else {}
        constraints.append(dict(kind='disjoint',first=first,second=second,parameters=params))
    for row in audit['all_fifteen_pair_constraints']:
        first,second=row['first'],row['second'];require(row['passed'] is True and row['actual_occupied_intersection_m3']==0
            and row['first_geometry_sha256']==meshes[first].geometry_sha and row['second_geometry_sha256']==meshes[second].geometry_sha,'Actual within-core pair changed')
        params=dict(occupancy='solid')
        if row['literal_shared_interface']:params['allow_coplanar_contact']=True
        constraints.append(dict(kind='disjoint',first=first,second=second,parameters=params))
        if row['literal_shared_interface']:
            def faces(a,b):return [f for g in ap[a]['tissueDrawGroups'] if g.get('neighbor')==b for f in range(g['start']//3,(g['start']+g['count'])//3)]
            constraints.append(dict(kind='shared-interface',first=first,second=second,parameters=dict(first_faces=faces(first,second),second_faces=faces(second,first))))
    inventory_path=ROOT/'reports/ta98-tongue-inventory.json';inventory=read(inventory_path);authority=inventory['authority'];byterm={g['term']:g for g in marker['termGroups']}
    require(inventory['term_count']==len(inventory['terms'])==36,'Complete authority inventory required')
    inventory.update(schema='ta98-geometry-inventory/v1',model='local-female-ta98',model_sex='female',root_terms=['A05.1.04.001','A04.1.05.001'],
        source_catalogue=str(FEMALE),source_catalogue_sha256=FEMALE_SHA,authority_path=authority['path'],authority_sha256=authority['sha256'],
        authority_corrections_path=authority['corrections_path'],authority_corrections_sha256=authority['corrections_sha256'],controls=14,adult_terms_missing=21)
    for row in inventory['terms']:
        g=byterm.get(row['authority']['id']);row.update(source_meshes=row['source_meshes']['female'],applicable_to_model=row.get('required_in_typical_adult',True),
            candidate_control=copy.deepcopy(g),candidate_control_absent=g is None,candidate_part_ids=g['partIds'] if g else [],anatomical_acceptance=False)
    bindings.update(packet['source_bank_bindings']);bindings.update(audit['source_bindings'])
    for path in [Path(__file__),ROOT/'server/build-tongue-female-core-only.py',inventory_path,ROOT/'ATTRIBUTION.md',Path(authority['path']),Path(authority['corrections_path']),
                 CANDIDATE.with_name('frontend-controls-v1.json'),CANDIDATE.with_name('review-surface-right-v2.json'),CANDIDATE.with_name('review-sagittal-v2.json')]:bindings[str(path.resolve())]=digest_file(path)
    require(all(digest_file(p)==sha for p,sha in bindings.items()),'An explicit emitted/target/source input changed')
    OUT.mkdir();atomic_json(OUT/'inventory.json',inventory)
    spec=dict(model='local-female-ta98',proposed_catalogue=str(CANDIDATE),expected_catalogue_sha256=FEMALE_SHA,
        changed_parts=sorted(active|suppressed),closed_parts=sorted(active),require_closed=False,max_pairs=250000,max_tests=2000000,
        inventory=str(OUT/'inventory.json'),interface_checks=constraints,evidence_files=sorted(bindings),explicit_input_bindings=bindings,
        stage='partial-local-review',external_deployment=False,assembly=dict(whole_organ=marker['wholeConceptId'],part_ids=marker['partIds'],
            fixed=['Entire skeleton','Every unrelated female source','All native extrinsic geometry preserved and suppressed','Entire installed male model'],
            representation_mode=marker['representationMode'],tissue_ownership_claim=False,anatomical_acceptance=False,
            adjoining_interfaces=[dict(kind='literal-shared-tissue',scope='Actual62 M/APO and122 APO/SL reciprocal faces and signed words; all15 pair occupied volumes zero'),
                dict(kind='actual-fixed-context',scope='All66 new-part versus fixed-context triangle checks passed and are rechecked by this job')],
            remaining_limits=[marker['coverageLimitation'],'Coarse exterior and gaps between reference muscle envelopes remain; this is not a complete tongue or accepted histology.']))
    atomic_json(OUT/'install-spec.json',spec);atomic_json(OUT/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec))
    print(json.dumps(dict(job_spec=str(OUT/'job-spec.json'),constraints=len(constraints),controls=14,meshes=6,adult_missing=21,
                          job_started=False,installed=False,elapsed_seconds=time.monotonic()-started)))
if __name__=='__main__':main()
