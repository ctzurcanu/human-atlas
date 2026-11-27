#!/usr/bin/env python3
"""Comparison only: replay frozen lower-limb source faces, never refit."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs, exact_transverse_crossing


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def area(tri):
    return np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2


def distribution(a):
    return {'count': len(a), 'minimum': float(a.min()), 'maximum': float(a.max()), 'quantiles': dict(zip(['p01', 'p05', 'p50', 'p95', 'p99'], map(float, np.quantile(a, [.01, .05, .5, .95, .99])))), 'belowHalf': int(np.sum(a < .5)), 'aboveDouble': int(np.sum(a > 2)), 'outsideTwentyPercent': int(np.sum((a < .8) | (a > 1.2)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-stage', required=True)
    parser.add_argument('--output', default='trial-v1-audit')
    args = parser.parse_args()
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/biodigital-male-circulation-v1'
    reference = base / 'lower-limb-joints-v1'
    candidate = base / args.candidate_stage / 'candidate-v1'
    baseline = base / 'body-envelope-fit-v2/candidate-v1'
    output = base / args.candidate_stage / args.output
    output.mkdir(exist_ok=False)
    evidence = json.loads((reference / 'constraints.json').read_text())
    sample_path = reference / 'vessel-sample-ownership.npz'
    frozen = np.load(sample_path)
    # Cache every compressed ownership array once before indexing in loops.
    mesh_ids = frozen['candidate_mesh']
    candidate_faces = frozen['candidate_face']
    source_faces = frozen['source_face']
    source_instances = frozen['source_instance_index']
    source_ids = frozen['source_instance_id']
    old_triangles = frozen['triangles_m']
    triangles = np.empty_like(old_triangles)
    rows = json.loads((candidate / 'instance-ownership.json').read_text())['instances']
    rows_by_id = {row['sourceInstanceId']: row for row in rows}
    bindings = {str(reference / 'constraints.json'): sha(reference / 'constraints.json'), str(sample_path): sha(sample_path), str(candidate / 'candidate-report.json'): sha(candidate / 'candidate-report.json'), str(candidate / 'instance-ownership.json'): sha(candidate / 'instance-ownership.json')}
    strains = []
    for mesh_name in ['arterial', 'venous']:
        path, ancestry = candidate / (mesh_name + '.npz'), candidate / (mesh_name + '-ancestry.npz')
        bindings[str(path)], bindings[str(ancestry)] = sha(path), sha(ancestry)
        z, az = np.load(path), np.load(ancestry)
        v, f = z['vertices'], z['faces']
        sf, si = az['candidate_face_source_face'], az['candidate_face_source_instance']
        original_ids = az['source_instance_ids']
        use = np.flatnonzero(mesh_ids == mesh_name)
        ids = candidate_faces[use]
        assert np.array_equal(sf[ids], source_faces[use]) and np.array_equal(si[ids], source_instances[use])
        assert np.array_equal(original_ids[si[ids]], source_ids[use])
        triangles[use] = v[f[ids]]
        old_z, old_a = np.load(baseline / (mesh_name + '.npz')), np.load(baseline / (mesh_name + '-ancestry.npz'))
        ov, of = old_z['vertices'], old_z['faces']
        old_sf, old_si = old_a['candidate_face_source_face'], old_a['candidate_face_source_instance']
        for row in rows:
            if row['candidateMesh'] != mesh_name + '.npz' or not any(word in row['name'].lower() for word in ['popliteal', 'genicular', 'femoral artery', 'femoral vein', 'anterior tibial artery', 'anterior tibial veins', 'posterior tibial artery', 'posterior tibial veins', 'fibular artery', 'fibular veins']):
                continue
            start, count = row['candidateFaceStart'], row['candidateFaceCount']
            cut = slice(start, start + count)
            assert np.array_equal(old_sf[cut], sf[cut]) and np.array_equal(old_si[cut], si[cut])
            ratio = area(v[f[cut]]) / area(ov[of[cut]])
            strains.append({'sourceInstanceId': row['sourceInstanceId'], 'name': row['name'], 'retainedFaceCount': count, 'triangleAreaRatio': distribution(ratio)})
        del z, az, v, f, sf, si, original_ids, old_z, old_a, ov, of, old_sf, old_si
    context = base / 'regional-interface-v1/fixed-context'
    for path in [context / 'full.npz', context / 'context.json', base / 'global-envelope-v1/fixed-skin/skin.npz']:
        bindings[str(path)] = sha(path)
        assert bindings[str(path)] == evidence['inputBindings'][str(path)], 'Fixed bone/skin evidence changed'
    bone_z = np.load(context / 'full.npz')
    full_owner, source_bone_face = bone_z['face_owner_part_index'], bone_z['face_owner_source_face_id']
    bone_names = {row['partIndex']: row['id'] for row in evidence['fixedBoneParts']}
    bone_map = np.flatnonzero(np.isin(full_owner, list(bone_names)))
    frame = str(bone_z['frame'])
    bones = Mesh(bone_z['vertices'], bone_z['faces'][bone_map], frame)
    sample = Mesh(triangles.reshape((-1, 3)), np.arange(len(triangles) * 3).reshape((-1, 3)), frame)
    a, b = BVH(sample), BVH(bones)
    pair_chunks, strict_chunks, pieces = [], [], {}
    tested, complete = 0, False
    for pairs in a.candidate_batches(b, .001, batch_size=2048):
        if tested + len(pairs) > 350000 or time.monotonic() - began > 30:
            break
        result = triangle_pairs(triangles[pairs[:, 0]], b.triangles[pairs[:, 1]])
        tested += len(pairs)
        take = np.flatnonzero(result['distance_m'] <= .001)
        if len(take):
            pp = pairs[take]
            strict = np.zeros(len(pp), dtype=bool)
            for k in np.flatnonzero(result['relation'][take] == 'crossing'):
                strict[k] = exact_transverse_crossing(triangles[pp[k, 0]], b.triangles[pp[k, 1]])
            pair_chunks.append(pp)
            strict_chunks.append(strict)
            for key, value in result.items():
                pieces.setdefault(key, []).append(value[take])
        (output / 'checkpoint.json').write_text(json.dumps({'status': 'comparison-only strict query', 'testedPairs': tested, 'elapsedSeconds': time.monotonic() - began}) + '\n')
    else:
        complete = True
    pairs = np.concatenate(pair_chunks) if pair_chunks else np.empty((0, 2), dtype=np.int64)
    strict = np.concatenate(strict_chunks) if strict_chunks else np.empty(0, dtype=bool)
    result = {k: np.concatenate(v) for k, v in pieces.items()} if pieces else triangle_pairs(np.empty((0, 3, 3)), np.empty((0, 3, 3)))
    np.savez_compressed(output / 'actual-bone-triangle-witnesses.npz', sampled_triangles_m=triangles, pairs=pairs, strict_literal_crossing=strict, fixed_full_skeleton_face=bone_map[pairs[:, 1]], fixed_original_source_face=source_bone_face[bone_map[pairs[:, 1]]], fixed_source_part_index=full_owner[bone_map[pairs[:, 1]]], **result)
    now = Counter()
    for i in np.flatnonzero(strict):
        sample_id, bone_id = pairs[i]
        now[(str(source_ids[sample_id]), bone_names[int(full_owner[bone_map[bone_id]])])] += 1
    previous = Counter((row['sourceInstanceId'], row['bonePart']) for row in evidence['strictWitnessOwnership'])
    changes = [{'sourceInstanceId': key[0], 'name': rows_by_id[key[0]]['name'], 'bonePart': key[1], 'before': previous[key], 'after': now[key], 'delta': now[key] - previous[key]} for key in previous.keys() | now.keys()]
    changes.sort(key=lambda row: (-row['after'], row['name']))
    knee = lambda name: any(term in name.lower() for term in ['popliteal', 'genicular'])
    knee_before = sum(value for key, value in previous.items() if knee(rows_by_id[key[0]]['name']))
    knee_after = sum(value for key, value in now.items() if knee(rows_by_id[key[0]]['name']))
    report = {'schema': 'frozen-source-face-lower-limb-trial-comparison/v1', 'candidateStage': args.candidate_stage, 'inputBindings': bindings, 'sameExactSourceFaceIds': True, 'sampledVesselFaces': len(triangles), 'fixedJointBoneFaces': len(bone_map), 'queryCompleteForSelectedFacesWithin1mm': complete, 'candidatePairsTested': tested, 'strictPairsBefore': sum(previous.values()), 'strictPairsAfter': sum(now.values()), 'strictPairDelta': sum(now.values()) - sum(previous.values()), 'kneeFamilyStrictPairsBefore': knee_before, 'kneeFamilyStrictPairsAfter': knee_after, 'kneeFamilyStrictPairDelta': knee_after - knee_before, 'byVesselBone': changes, 'actualTriangleAreaStrain': strains, 'areaStrainQualification': 'Full retained named major trunk/genicular surface face ratios compared with saved body-envelope-v2 candidate; this is triangle surface strain, not lumen calibre, thickness or patency.', 'criterion': {'majorKneeSampledCrossingsReduced': complete and knee_after < knee_before, 'geometryAndInsideAcceptance': False, 'actualSkinInsideClearanceChecksStillRequired': True}, 'noRefit': True, 'noNewFieldsProposed': True, 'geometryModified': False, 'installed': False, 'elapsedSeconds': time.monotonic() - began}
    (output / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    assert all(sha(path) == value for path, value in bindings.items()), 'Bound input changed'
    (output / 'checkpoint.json').write_text(json.dumps({'status': 'complete' if complete else 'saved incomplete', 'elapsedSeconds': time.monotonic() - began}) + '\n')
    print(json.dumps({'output': str(output), 'complete': complete, 'strictBefore': sum(previous.values()), 'strictAfter': sum(now.values()), 'kneeBefore': knee_before, 'kneeAfter': knee_after, 'topChanges': changes[:12], 'elapsedSeconds': time.monotonic() - began}))


if __name__ == '__main__':
    main()
