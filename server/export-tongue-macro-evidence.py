"""Normalize saved ROI evidence for packing without altering its construction."""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file


def main(candidate,base_construction=None):
    began=time.monotonic();work=candidate.parent
    construction_path=candidate/'construction.json';construction_sha=digest_file(construction_path)
    record=json.loads(construction_path.read_text())
    assert record['numerical_and_resource_passed'] and not record['exclusive_tissue_partition_claimed']
    assert record['source_geometry_modified'] is False and record['anatomical_acceptance'] is False
    source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(source)==record['source_catalogue_sha256']
    native=Mesh.load(work/'source-v1/meshes/native-00.npz')
    assert native.source['catalogue_sha256']==record['source_catalogue_sha256']
    source_audit=work/'source-v1/source-audit.json'
    old_context_path=work/'intrinsics-v5/retained-context-audit/summary.json'
    old_context=json.loads(old_context_path.read_text())
    assert digest_file(source_audit)==old_context['source_audit_sha256']
    stylo=work/'styloglossus-v2/candidate.json'
    assert digest_file(stylo)==old_context['styloglossus_candidate_sha256']
    mode='macroscopic-muscle-extent-reference-envelope'
    common=dict(representation_mode=mode,tissue_ownership_claim=False,roi_boundaries='approximate',
        anatomical_acceptance=False,complete_muscle_volume_segmentation=False,installed=False,
        source_geometry_modified=False)
    normal_bindings={}
    normal_rows=None
    if 'presentation_normal_proof_sha256'in record:
        normal_path=candidate/'normal-only-proof.json';normal_proof=json.loads(normal_path.read_text())
        assert digest_file(normal_path)==record['presentation_normal_proof_sha256']
        assert normal_proof['passed'] and normal_proof['geometry_unchanged']
        assert base_construction is not None
        assert digest_file(base_construction)==record['base_construction_sha256']==normal_proof['base_construction_sha256']
        normal_rows={r['key']:r for r in normal_proof['meshes']}
        normal_bindings=dict(presentation_normal_proof=dict(path=str(normal_path),sha256=digest_file(normal_path)),
            base_construction=dict(path=str(base_construction),sha256=digest_file(base_construction)),
            geometry_checks_reused_from_exact_geometry=True,presentation_normal_side_proof_passed=True)
    metadata={r['key']:r for r in record['meshes']};checks={r['key']:r for r in record['checks']}
    assert len(metadata)==4 and set(metadata)==set(checks)
    rows=[]
    for key,row in metadata.items():
        assert digest_file(row['mesh'])==row['sha256']
        mesh=Mesh.load(row['mesh']);assert mesh.geometry_sha==row['geometry_sha256']
        if normal_rows is not None:
            proof=normal_rows[key]
            assert proof['variant_mesh_sha256']==row['sha256']
            assert proof['variant_geometry_sha256']==mesh.geometry_sha==proof['base_geometry_sha256']
            assert proof['side_orientation_proven']and proof['variant_nonpositive_face_corners']==0
        check=checks[key];assert check['passed'] and all(check['checks'].values())
        assert check['cpp_diagnostics']['geometry_sha256']==mesh.geometry_sha
        rows.append(dict(key=key,sha256=row['sha256'],geometry_sha256=mesh.geometry_sha,
            numerical_checks_passed=True,topology=check['topology'],cpp_diagnostics=check['cpp_diagnostics'],
            exterior_crossings=check['source_contacts'],minimum_absolute_winding=check['minimum_absolute_winding'],
            whole_tongue_contained=check['checks']['whole_tongue_contained'],
            source_crossing_check_complete=check['checks']['source_budget_complete'],
            actual_check_flags=check['checks'],anatomical_acceptance=False))
    intrinsic=[r for r in record['pair_diagnostics']if r['scope'].startswith('intrinsic')]
    assert len(intrinsic)==6 and all(r['complete_within_pair_budget']for r in intrinsic)
    defects=sum(r['crossing_face_pairs']+r['coplanar_overlap_face_pairs']for r in intrinsic)
    assert defects==0
    checked=dict(schema='tongue-intrinsic-numerical-audit/v1',**common,**normal_bindings,
        source_evidence='Existing actual emitted-mesh checks from construction.json; no invented clearance or query counts',
        passed=True,construction_sha256=construction_sha,source_catalogue_sha256=record['source_catalogue_sha256'],
        source_tongue_geometry_sha256=native.geometry_sha,source_bank_sha256=native.source['bank_sha256'],
        rows=rows,cross_muscle_contacts=intrinsic,cross_muscle_defect_pairs=defects,
        visual_review='pending',resource=record['resource'],storage_origin=record['storage_origin'],
        exported_by_script_sha256=digest_file(__file__))
    atomic_json(candidate/'checked-geometry.json',checked)
    context=[]
    for row in old_context['context']:
        assert digest_file(row['mesh'])==row['sha256']
        mesh=Mesh.load(row['mesh']);assert mesh.geometry_sha==row['geometry_sha256']
        context.append(row)
    assert len(context)==6
    external=[r for r in record['pair_diagnostics']if not r['scope'].startswith('intrinsic')]
    expected={(a,b['key'])for a in metadata for b in context}
    assert len(external)==24 and {(r['first'],r['second'])for r in external}==expected
    assert all(r['complete_within_pair_budget']for r in external)
    external_rows=[dict(intrinsic=r['first'],retained=r['second'],contact_summary=r,
        intrinsic_geometry_sha256=metadata[r['first']]['geometry_sha256'],
        retained_geometry_sha256=next(c['geometry_sha256']for c in context if c['key']==r['second']),
        interpretation='Measured contact between nonexclusive macroscopic extent reference envelopes; not an exclusive tissue-ownership test')
        for r in external]
    context_report=dict(schema='tongue-intrinsic-retained-context-audit/v1',**common,**normal_bindings,
        passed=True,passed_scope='Measurement completed with source/geometry hashes verified; not tissue ownership or anatomical acceptance',
        measurement_complete=True,assembly_exclusive_tissue_verified=False,
        source_catalogue_sha256=record['source_catalogue_sha256'],source_audit_sha256=digest_file(source_audit),
        construction_sha256=construction_sha,checked_geometry_sha256=digest_file(candidate/'checked-geometry.json'),
        candidate_check_sha256=digest_file(candidate/'checked-geometry.json'),styloglossus_candidate_sha256=digest_file(stylo),
        context=context,pairs=external_rows,pair_count=len(external_rows),expected_pair_count=24,
        complete_within_pair_budget=True,
        crossing_face_pairs=sum(r['crossing_face_pairs']for r in external),
        coplanar_overlap_face_pairs=sum(r['coplanar_overlap_face_pairs']for r in external),
        source_measurement_construction_sha256=construction_sha,
        source_measurement_elapsed_seconds=record['elapsed_seconds'],
        exported_by_script_sha256=digest_file(__file__))
    context_dir=candidate/'retained-context-audit';context_dir.mkdir(exist_ok=True)
    atomic_json(context_dir/'summary.json',context_report)
    coverage=('Original approximate macroscopic muscle-extent reference envelope containing interwoven fibres; '
        'not patient-specific histological segmentation or an exclusive tissue volume. T/V uses one shared carrier.')
    packing=dict(schema='tongue-macro-envelope-packing-input/v1',**common,**normal_bindings,
        construction_sha256=construction_sha,source_catalogue_sha256=record['source_catalogue_sha256'],
        source_audit_sha256=digest_file(source_audit),source_bank_sha256=native.source['bank_sha256'],
        source_tongue_geometry_sha256=native.geometry_sha,
        all_canonical_muscle_terms=['A05.1.04.106','A05.1.04.107','A05.1.04.108','A05.1.04.109'],
        checked_geometry_sha256=digest_file(candidate/'checked-geometry.json'),
        retained_context_summary_sha256=digest_file(context_dir/'summary.json'),storage_origin=record['storage_origin'],
        retained_context_policy=dict(allowed_intersection_part_ids=[r['key']for r in context],
            primary_reference_url=record['primary_reference'],interpretation='Nonexclusive intralingual muscle-extent reference envelopes'),
        intrinsic_mesh_overlays={key:dict(partId='TA98TONGUE:male:intrinsic:'+key,coverageLimitation=coverage,
            terms=row['terms'])for key,row in metadata.items()},
        resource=record['resource'],exported_by_script_sha256=digest_file(__file__))
    atomic_json(candidate/'packing-input.json',packing)
    assert digest_file(construction_path)==construction_sha
    print(json.dumps(dict(construction_sha256=construction_sha,construction_unchanged=True,
        checked_rows=len(rows),context_pairs=len(external_rows),actual_context_crossing_pairs=context_report['crossing_face_pairs'],
        representation_mode=mode,elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True,type=Path)
    ap.add_argument('--base-construction',type=Path);args=ap.parse_args()
    main(args.candidate.resolve(),args.base_construction.resolve()if args.base_construction else None)
