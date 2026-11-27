"""Bind reviewed curved intrinsics to the existing packer without altering them.

Adds a redundant literal face-corner normal field only. No physical triangles,
normal words/signs, fixed covering banks or ordinary atlas files are rewritten.
This establishes numerical compatibility, not anatomical or visual acceptance.
"""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'server')]
import numpy as np
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh
from tongue_basal_contract import FRAME, INTRINSIC_KEYS, MODE, corner_normal_mesh, require
from tongue_whole_input import whole_input

WORK = ROOT/'.local-models/geometry-workbench/tongue-completion'
ORIGIN = np.array([0., 1.625, 0.])


def read(path):
    return json.loads(Path(path).read_text())


def bound(path, **extra):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=digest_file(path), **extra)


def load_packer():
    spec = importlib.util.spec_from_file_location('curved_intrinsic_packer', ROOT/'server/build-tongue-coupled-candidate.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--construction', type=Path, default=WORK/'tv-central-core-dual-exclusion-v2-storage-patch/construction.json')
    parser.add_argument('--assessment', type=Path, default=WORK/'tv-central-core-dual-exclusion-v2-storage-patch/numerical-assessment.json')
    parser.add_argument('--selected', type=Path, default=WORK/'selected-coupling-audit-v2/coupled-inputs.selected.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    source, assessment, selected = read(args.construction), read(args.assessment), read(args.selected)
    require(not args.output.exists(), 'Use a new saved adapter stage')
    require(assessment['construction_sha256'] == digest_file(args.construction)
            and assessment['all_actual_geometric_and_normal_gates_passed'] is True
            and assessment['bone_pairs_complete'] == 32
            and assessment['covering_pairs_complete'] == 6
            and assessment['intrinsic_pairs_complete'] == 6,
            'Actual completed geometric evidence is required')
    geometric_construction_sha = digest_file(args.construction)
    evidence = []
    if assessment.get('geometric_evidence_predecessor'):
        predecessor = assessment['geometric_evidence_predecessor']
        for binding in predecessor.values():
            require(digest_file(binding['path']) == binding['sha256'], 'Literal normal-only predecessor changed')
        parent = read(predecessor['construction']['path'])
        parent_assessment = read(predecessor['assessment']['path'])
        require(parent_assessment['construction_sha256'] == predecessor['construction']['sha256']
                and parent_assessment['all_actual_geometric_and_normal_gates_passed'] is True
                and assessment['check_bindings'] == parent_assessment['check_bindings'],
                'Normal-only replay lacks complete original physical checks')
        previous = {r['key']:r for r in parent['meshes']}
        require(set(previous) == INTRINSIC_KEYS, 'Normal-only replay omitted a physical carrier')
        for row in source['meshes']:
            old = previous[row['key']]
            require(digest_file(old['mesh']) == old['sha256'] and digest_file(row['mesh']) == row['sha256'],
                    'Literal physical replay source changed')
            before, after = Mesh.load(old['mesh']), Mesh.load(row['mesh'])
            require(np.array_equal(before.vertices,after.vertices) and np.array_equal(before.faces,after.faces)
                    and before.geometry_sha == after.geometry_sha,
                    'Physical checks may only be reused for exactly identical indexed geometry')
            tri=after.triangles;normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
            corner=after.attributes['normals'][after.faces].astype(float)*after.attributes['normalSigns'][after.faces,None]
            require(np.all(np.einsum('fci,fi->fc',corner,normal)>0), 'Replayed normal-only field points through a face')
        correction=assessment['normal_field_correction']
        require(digest_file(correction['path']) == correction['sha256'] and read(correction['path'])['passed'] is True,
                'Normal-only correction proof is stale or failed')
        geometric_construction_sha=predecessor['construction']['sha256']
        evidence.extend([*predecessor.values(), correction])
    for path, sha in assessment['check_bindings'].items():
        p = Path(path) if Path(path).is_absolute() else ROOT/path
        doc = read(p)
        if doc.get('stage') == 'retained-SL-self':
            sl = next(r for r in source['meshes'] if r['key'] == 'superior-longitudinal')
            require(doc.get('source_mesh') == sl['mesh'] and doc.get('source_sha256') == sl['sha256']
                    and doc['checks'][0]['result']['geometry_sha256'] == sl['geometry_sha256'],
                    'Retained SL self check is not bound to the literal unchanged source')
        else:
            require(doc.get('construction_sha256') == geometric_construction_sha, 'Actual check binds another construction')
        require(digest_file(p) == sha and doc.get('passed') is True, 'Actual numerical check is stale or failed')
        evidence.append(bound(p))
    require({r['key'] for r in source['meshes']} == INTRINSIC_KEYS and len(source['meshes']) == 4,
            'Exactly four intrinsic carriers required')
    baseline_path = Path(selected['bindings']['intrinsics']['path'])
    covering_path = Path(selected['bindings']['coverings']['path'])
    require(digest_file(baseline_path) == selected['bindings']['intrinsics']['sha256']
            and digest_file(covering_path) == selected['bindings']['coverings']['sha256'],
            'Explicit selected semantics/coverings changed')
    baseline, coverings = read(baseline_path), read(covering_path)
    _, _, actual_whole = whole_input(Path(selected['bindings']['whole']['path']), Path(selected['bindings']['whole_fit']['path']))
    require(baseline['whole_input'] == coverings['whole_input'] == actual_whole, 'Different corrected-whole source')
    metadata = {r['key']:r for r in baseline['meshes']}
    args.output.mkdir(parents=True)
    packer = load_packer()
    rows, proofs, all_meshes = [], [], {}
    for original in source['meshes']:
        key, path = original['key'], Path(original['mesh'])
        require(digest_file(path) == original['sha256'], 'Original intrinsic source changed')
        mesh = Mesh.load(path)
        require(mesh.frame == FRAME and mesh.geometry_sha == original['geometry_sha256'], 'Intrinsic physical source changed')
        row = {name:copy.deepcopy(metadata[key][name]) for name in
               ('color', 'terms', 'representation', 'exclusive_tissue_volume', 'anatomical_acceptance')}
        row.update(key=key, triangles=len(mesh.faces), vertices=len(mesh.vertices), components=1)
        if key == 'superior-longitudinal':
            manifest = read(coverings['packed_manifest'])
            sl_id = manifest['part_ids'][key]
            sl_path = Path(coverings['packed_manifest']).with_name(sl_id.replace(':','_')+'.npz')
            literal = Mesh.load(sl_path)
            require(np.array_equal(mesh.triangles,literal.triangles)
                    and np.array_equal(mesh.attributes['normals'][mesh.faces],literal.attributes['normals'][literal.faces])
                    and np.array_equal(mesh.attributes['normalSigns'][mesh.faces],literal.attributes['normalSigns'][literal.faces]),
                    'Retained shared-bank SL differs from actual assessed geometry or raw normals')
            row.update(mesh=str(sl_path.resolve()),sha256=digest_file(sl_path),geometry_sha256=literal.geometry_sha,
                       unchanged_literal=True,shared_covering_bank_retained=True)
            proofs.append(dict(key=key, original=bound(path),selected=bound(sl_path),
                               literal_physical_triangles=True,literal_raw_normal_words=True,literal_normal_signs=True))
            all_meshes[key] = literal
            rows.append(row)
            continue
        words, signs = packer.old.normal_words(mesh)
        require(np.all(signs == 1), 'Current corner adapter must preserve +1 source signs literally')
        with np.load(path, allow_pickle=False) as raw:
            arrays = {name:raw[name].copy() for name in raw.files}
        require('corner_normals' not in arrays, 'Source unexpectedly has a second corner-normal field')
        arrays['corner_normals'] = words[mesh.faces].copy()
        target = args.output/(key+'.npz')
        atomic_npz(target, **arrays)
        with np.load(target, allow_pickle=False) as raw:
            require(set(raw.files) == set(arrays) and all(np.array_equal(raw[name],value) for name,value in arrays.items()),
                    'Adapter changed a literal source NPZ array')
        adapted = Mesh.load(target)
        require(np.array_equal(mesh.vertices,adapted.vertices) and np.array_equal(mesh.faces,adapted.faces)
                and mesh.geometry_sha == adapted.geometry_sha and mesh.attributes_sha == adapted.attributes_sha,
                'Adapter must retain indexed coordinates/faces and all attributes')
        expanded = corner_normal_mesh(target,mesh.geometry_sha)
        expanded_words, expanded_signs = packer.old.normal_words(expanded)
        require(np.array_equal(expanded.triangles,mesh.triangles)
                and np.array_equal(expanded_words[expanded.faces],words[mesh.faces])
                and np.array_equal(expanded_signs[expanded.faces],signs[mesh.faces]),
                'Real packing normal expansion changed triangles/encoded normals/signs')
        decoded = (expanded.vertices-ORIGIN).astype('<f4').astype(float)+ORIGIN
        legacy = (expanded.vertices-packer.old.ORIGIN).astype('<f4').astype(float)+packer.old.ORIGIN
        require(np.array_equal(decoded,expanded.vertices), 'Declared dyadic storage origin does not replay literal coordinates')
        origin_proof = dict(schema='tongue-literal-intrinsic-origin-replay/v1',
            mesh=bound(target,geometry_sha256=mesh.geometry_sha),geometryOrigin=ORIGIN.tolist(),
            literal_decoded_replay=True,source=bound(path),
            literal_original_NPZ_arrays=True,literal_physical_triangles=True,
            literal_raw_normal_words=True,literal_normal_signs=True,
            legacy_origin=packer.old.ORIGIN.tolist(),legacy_literal_decoded_replay=bool(np.array_equal(legacy,expanded.vertices)),
            legacy_changed_coordinate_count=int(np.count_nonzero(legacy != expanded.vertices)),
            legacy_maximum_position_error_m=float(np.max(np.linalg.norm(legacy-expanded.vertices,axis=1))))
        proof_path = args.output/(key+'-origin-proof.json')
        atomic_json(proof_path,origin_proof)
        row.update(mesh=str(target.resolve()),sha256=digest_file(target),geometry_sha256=mesh.geometry_sha,
                   geometryOrigin=ORIGIN.tolist(),storageOriginProof=bound(proof_path),
                   corner_normal_source=dict(**bound(target),array='corner_normals'),source=bound(path))
        actual_origin, strict = packer.intrinsic_storage_origin(row,expanded)
        require(strict and np.array_equal(actual_origin,ORIGIN), 'Actual packer failed explicit literal origin contract')
        bank = bytearray()
        packed = packer.pack_record(bank,expanded,dict(id='test:'+key),actual_origin,literal=strict)
        offset=packed['positions']
        positions=np.frombuffer(bank,dtype='<f4',count=expanded.vertices.size,offset=offset).reshape(-1,3)
        require(np.array_equal(positions.astype(float)+actual_origin,expanded.vertices), 'Actual packed buffer changed coordinates')
        origin_proof.update(key=key,origin_proof=bound(proof_path),actual_packed_bytes=len(bank),
                            actual_packer_literal_replay_passed=True,expanded_vertices=len(expanded.vertices))
        proofs.append(origin_proof)
        rows.append(row)
        all_meshes[key]=expanded
    changed_bytes=sum(r.get('actual_packed_bytes',0) for r in proofs)
    all_triangles=sum(len(m.faces) for m in all_meshes.values())
    require(changed_bytes <= 2000000 and all_triangles <= 30000, 'Intrinsic resources exceed the existing limits')
    report=dict(schema='tongue-curved-intrinsic-packing-adapter/v1',passed=True,
                source_construction=bound(args.construction),assessment=bound(args.assessment),
                selected_semantics=bound(baseline_path),coverings=bound(covering_path),
                checks=evidence,meshes=proofs,changed_intrinsic_packed_bytes=changed_bytes,
                total_intrinsic_triangles=all_triangles,
                coupled_covering_context_audits_must_be_regenerated=True,
                actual_surface_section_review_required=True,anatomical_acceptance=False,installed=False,
                elapsed_seconds=time.monotonic()-started)
    proof_path=args.output/'adapter-proof.json'
    atomic_json(proof_path,report)
    construction=dict(schema='tongue-macroscopic-muscle-reference-envelopes/v1',
        frame=FRAME,whole_input=actual_whole,source_catalogue_sha256=actual_whole['source_catalogue_sha256'],
        primary_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
        shared_tv_representation=copy.deepcopy(baseline['shared_tv_representation']),meshes=rows,
        numerical_and_resource_passed=True,intrinsic_pairs_crossing_free=True,
        fixed_bone_crossings_zero=True,fixed_bone_pairs_complete=True,normal_chart_variant_passed=True,
        exclusive_tissue_partition_claimed=False,exact_intrinsic_region_boundary_claim=False,
        anatomical_acceptance=False,installed=False,visual_review='pending',
        source_geometry_modified=False,fixed_bone_geometry_modified=False,
        numerical_evidence=bound(args.assessment),normal_and_storage_adapter=bound(proof_path),
        original_construction=bound(args.construction),
        resource=dict(changed_intrinsic_bank_bytes=changed_bytes,total_intrinsic_triangles=all_triangles,
                      whole_muscle_bank_limit_bytes=2000000,whole_muscle_bank_limit_triangles=30000,
                      complete_coupled_bank_still_requires_actual_packing=True))
    out=args.output/'construction.json'
    atomic_json(out,construction)
    packer.qualified_intrinsics(construction,out)
    print(json.dumps(dict(construction=str(out.resolve()),sha256=digest_file(out),proof=str(proof_path.resolve()),
                          elapsed_seconds=time.monotonic()-started,changed_bytes=changed_bytes,triangles=all_triangles)))


if __name__ == '__main__':
    main()
