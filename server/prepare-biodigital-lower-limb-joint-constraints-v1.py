#!/usr/bin/env python3
"""One bounded actual-triangle lower-limb constraint stage; no geometry edit."""
import json
import hashlib
from collections import Counter
from pathlib import Path
import re
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs, exact_transverse_crossing
from geometry_tools.algorithms.compact_positive_fields import compact_radial_map


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def section_points(triangles, y):
    points = []
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        a, b = triangles[:, i], triangles[:, j]
        take = ((a[:, 1] < y) & (b[:, 1] > y)) | ((a[:, 1] > y) & (b[:, 1] < y))
        a, b = a[take], b[take]
        if len(a):
            mu = (y - a[:, 1]) / (b[:, 1] - a[:, 1])
            points.extend(a + mu[:, None] * (b - a))
    if not points:
        return np.empty((0, 3))
    points = np.asarray(points)
    _, ids = np.unique(np.round(points / 1e-10).astype(np.int64), axis=0, return_index=True)
    return points[np.sort(ids)]


def z_intersections(tree, point):
    """All projected XY triangle hits; actual source face/barycentric witnesses."""
    x, y = point[:2]
    ids = tree.query_bounds(np.array([x - 1e-12, y - 1e-12, tree.low[0, 2]]), np.array([x + 1e-12, y + 1e-12, tree.high[0, 2]]))
    t = tree.triangles[ids]
    if not len(t):
        return []
    a, u, v = t[:, 0, :2], t[:, 1, :2] - t[:, 0, :2], t[:, 2, :2] - t[:, 0, :2]
    p = np.array([x, y]) - a
    det = u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
    s = np.divide(p[:, 0] * v[:, 1] - p[:, 1] * v[:, 0], det, out=np.full(len(t), np.inf), where=det != 0)
    q = np.divide(u[:, 0] * p[:, 1] - u[:, 1] * p[:, 0], det, out=np.full(len(t), np.inf), where=det != 0)
    w = np.c_[1 - s - q, s, q]
    take = np.isfinite(w).all(axis=1) & (w.min(axis=1) >= -1e-10)
    zs = np.einsum('ni,ni->n', w[take], t[take, :, 2])
    return sorted([{"face": int(f), "z": float(z), "barycentric": ww.tolist()} for f, z, ww in zip(ids[take], zs, w[take])], key=lambda r: r["z"])


def main():
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/biodigital-male-circulation-v1'
    out = base / 'lower-limb-joints-v1'
    out.mkdir(exist_ok=True)
    stage = base / 'body-envelope-fit-v2'
    candidate = stage / 'candidate-v1'
    context = base / 'regional-interface-v1/fixed-context'
    skin_path = base / 'global-envelope-v1/fixed-skin/skin.npz'
    paths = [stage / 'contract.json', candidate / 'candidate-report.json', candidate / 'instance-ownership.json', context / 'context.json', context / 'full.npz', skin_path]
    bindings = {str(p): sha(p) for p in paths}
    contract, candidate_report = json.loads(paths[0].read_text()), json.loads(paths[1].read_text())
    fixed_context = json.loads(paths[3].read_text())
    assert candidate_report['contractSha256'] == bindings[str(paths[0])]
    assert contract['globalBodyRegistration']['targetSkinSha256'] == bindings[str(skin_path)]
    assert contract['target_catalogue_sha256'] == fixed_context['targetCatalogueSha256'] == '46262d8ab07682fcae3ff0ee634f652abad766c1e010f6e2ec311755aea9abf5'
    rows = json.loads(paths[2].read_text())['instances']
    source_instances = json.loads((base / 'source-v2/instances.json').read_text())['instances']
    source_by_id = {r['id']: r for r in source_instances}
    terms = ['popliteal', 'genicular', 'femoral', 'tibial', 'fibular', 'peroneal', 'malleolar', 'dorsalis pedis', 'circumflex iliac']
    selected_rows = [r for r in rows if any(t in r['name'].lower() for t in terms)]
    main_route = re.compile(r'^(Left|Right) (popliteal (artery|vein)|femoral (artery|vein)|anterior tibial (artery|veins)|posterior tibial (artery|veins)|fibular (artery|veins))$', re.I)
    bone_rows = [r for r in fixed_context['partOwnership'] if any(t in r['name'].lower() for t in ['hip bone', 'sacrum', 'femur', 'tibia', 'fibula', 'patella', 'talus', 'calcaneus'])]
    bone_by_index = {r['partIndex']: r for r in bone_rows}
    z = np.load(paths[4])
    full_vertices, full_faces = z['vertices'], z['faces']
    full_owner, full_source_face = z['face_owner_part_index'], z['face_owner_source_face_id']
    bone_map = np.flatnonzero(np.isin(full_owner, list(bone_by_index)))
    frame = str(z['frame'])
    bones = Mesh(full_vertices, full_faces[bone_map], frame)
    bone_tree = BVH(bones)
    sz = np.load(skin_path)
    skin = Mesh(sz['vertices'], sz['faces'], frame)
    skin_tree = BVH(skin)
    triangles, sample_rows, popliteal, plot_triangles = [], [], {}, []
    for mesh_name in ['arterial', 'venous']:
        path, ancestry_path = candidate / (mesh_name + '.npz'), candidate / (mesh_name + '-ancestry.npz')
        bindings[str(path)], bindings[str(ancestry_path)] = sha(path), sha(ancestry_path)
        mz, az = np.load(path), np.load(ancestry_path)
        v, f = mz['vertices'], mz['faces']
        sf, si, ids = az['candidate_face_source_face'], az['candidate_face_source_instance'], az['source_instance_ids']
        for row in selected_rows:
            if row['candidateMesh'] != mesh_name + '.npz':
                continue
            start, count = row['candidateFaceStart'], row['candidateFaceCount']
            n = count if main_route.match(row['name']) else min(96 if 'genicular' in row['name'].lower() else 32, count)
            sample_ids = start + np.unique(np.linspace(0, count - 1, n).astype(np.int64))
            for face in sample_ids:
                assert str(ids[si[face]]) == row['sourceInstanceId']
                triangles.append(v[f[face]])
                sample_rows.append({"sourceInstanceId": row['sourceInstanceId'], "name": row['name'], "candidateMesh": mesh_name, "candidateFace": int(face), "sourceFace": int(sf[face]), "sourceInstanceIndex": int(si[face]), "role": 'nutrient opening candidate; verify actual opening' if 'nutrient' in row['name'].lower() else 'major trunk' if main_route.match(row['name']) else 'regional branch; no bony-opening exemption'})
            if re.match(r'^(Left|Right) popliteal (artery|vein)$', row['name'], re.I):
                popliteal[row['name']] = (v[f[start:start + count]], row)
            if any(t in row['name'].lower() for t in ['popliteal', 'genicular', 'femoral artery', 'femoral vein', 'tibial artery', 'tibial veins']):
                plot_triangles.extend(v[f[start:start + count]])
        del mz, az, v, f, sf, si, ids
    triangles = np.asarray(triangles)
    assert len(triangles) <= 32768, len(triangles)
    sample = Mesh(triangles.reshape((-1, 3)), np.arange(len(triangles) * 3).reshape((-1, 3)), frame)
    tree = BVH(sample)
    np.savez_compressed(out / 'vessel-sample-ownership.npz', triangles_m=triangles, candidate_mesh=np.asarray([r['candidateMesh'] for r in sample_rows]), candidate_face=np.asarray([r['candidateFace'] for r in sample_rows]), source_face=np.asarray([r['sourceFace'] for r in sample_rows]), source_instance_index=np.asarray([r['sourceInstanceIndex'] for r in sample_rows]), source_instance_id=np.asarray([r['sourceInstanceId'] for r in sample_rows]))
    checkpoint = out / 'checkpoint.json'
    pieces, pair_chunks, strict_chunks, tested, query_complete = {}, [], [], 0, False
    for pairs in tree.candidate_batches(bone_tree, .001, batch_size=2048):
        if tested + len(pairs) > 350000 or time.monotonic() - began > 27:
            break
        result = triangle_pairs(triangles[pairs[:, 0]], bone_tree.triangles[pairs[:, 1]])
        tested += len(pairs)
        take = np.flatnonzero(result['distance_m'] <= .001)
        if len(take):
            pp = pairs[take]
            strict = np.zeros(len(take), dtype=bool)
            for k in np.flatnonzero(result['relation'][take] == 'crossing'):
                strict[k] = exact_transverse_crossing(triangles[pp[k, 0]], bone_tree.triangles[pp[k, 1]])
            pair_chunks.append(pp)
            strict_chunks.append(strict)
            for key, value in result.items():
                pieces.setdefault(key, []).append(value[take])
        checkpoint.write_text(json.dumps({"status": "actual joint triangle query", "testedPairs": tested, "elapsedSeconds": time.monotonic() - began}) + '\n')
    else:
        query_complete = True
    pairs = np.concatenate(pair_chunks) if pair_chunks else np.empty((0, 2), dtype=np.int64)
    strict = np.concatenate(strict_chunks) if strict_chunks else np.empty(0, dtype=bool)
    result = {k: np.concatenate(v) for k, v in pieces.items()} if pieces else triangle_pairs(np.empty((0, 3, 3)), np.empty((0, 3, 3)))
    np.savez_compressed(out / 'actual-bone-triangle-witnesses.npz', pairs=pairs, strict_literal_crossing=strict, fixed_full_skeleton_face=bone_map[pairs[:, 1]], fixed_original_source_face=full_source_face[bone_map[pairs[:, 1]]], fixed_source_part_index=full_owner[bone_map[pairs[:, 1]]], **result)
    strict_rows = []
    for i in np.flatnonzero(strict):
        s, b = pairs[i]
        ff = bone_map[b]
        strict_rows.append({**sample_rows[s], "sampleIndex": int(s), "bonePart": bone_by_index[int(full_owner[ff])]['id'], "fixedContextFace": int(ff), "originalBoneSourceFace": int(full_source_face[ff]), "centroid_m": triangles[s].mean(axis=0).tolist()})
    guide_rows, guide_points, guide_meta = [], [], []
    for name, (t, row) in popliteal.items():
        ymin, ymax = t[:, :, 1].min(), t[:, :, 1].max()
        lo, hi = max(ymin + .001, .405), min(ymax - .001, .555)
        for y in np.linspace(lo, hi, 5) if hi > lo else []:
            points = section_points(t, float(y))
            if not len(points):
                continue
            indices = np.unique(np.linspace(0, len(points) - 1, min(24, len(points))).astype(int))
            for p in points[indices]:
                bh, sh = z_intersections(bone_tree, p), z_intersections(skin_tree, p)
                q = {"vessel": name, "sourceInstanceId": row['sourceInstanceId'], "point_m": p.tolist(), "sectionY_m": float(y), "actualSectionPointNotAnatomicalLandmark": True}
                if bh:
                    ff = bone_map[bh[0]['face']]
                    q.update(bonePosteriorZ_m=bh[0]['z'], bonePart=bone_by_index[int(full_owner[ff])]['id'], fixedContextFace=int(ff), originalBoneSourceFace=int(full_source_face[ff]), boneBarycentric=bh[0]['barycentric'], posteriorShiftFor1mmBoneZSeparation_m=bh[0]['z'] - .001 - float(p[2]))
                if sh:
                    q.update(skinPosteriorZ_m=sh[0]['z'], skinAnteriorZ_m=sh[-1]['z'], skinPosteriorFace=sh[0]['face'], minimumShiftFor2mmSkinZSeparation_m=sh[0]['z'] + .002 - float(p[2]))
                guide_rows.append(q)
                guide_points.append(p)
                guide_meta.append(name)
    guide_points = np.asarray(guide_points).reshape((-1, 3))
    fields, correction = [], []
    for side in ['Right', 'Left']:
        gs = [q for q in guide_rows if q['vessel'].startswith(side)]
        bone_constraints = [q['posteriorShiftFor1mmBoneZSeparation_m'] for q in gs if 'bonePosteriorZ_m' in q]
        skin_constraints = [q['minimumShiftFor2mmSkinZSeparation_m'] for q in gs if 'skinPosteriorZ_m' in q]
        if not bone_constraints:
            correction.append({"side": side, "status": "No projected bony corridor obstruction measured at these source section points"})
            continue
        required = min(bone_constraints)
        allowed = max(skin_constraints) if skin_constraints else None
        item = {"side": side, "measuredFullPosteriorShiftForSectionBoneZMargin_m": required, "mostRestrictiveSkinZShiftLowerBound_m": allowed, "guidePointCount": len(gs), "fullShiftIsNotAnAnatomicalTarget": True}
        if required >= 0:
            item['status'] = 'No posterior correction requested by these measured section guides'
        elif allowed is None or allowed >= 0:
            item['status'] = 'Posterior trial withheld; measured skin corridor leaves no verified posterior shift budget'
        else:
            d = max(-.008, .5 * required, allowed + .0005)
            if d >= 0:
                item['status'] = 'Posterior trial withheld by skin margin'
            else:
                centre = np.mean([q['point_m'] for q in gs], axis=0)
                field = dict(center=centre.tolist(), core_m=.055, support_m=.16, kind='translation', translation_m=[0., 0., float(d)])
                mapped, jac, proof = compact_radial_map(guide_points, **field)
                fields.append(field)
                item.update(status='Single damped common-field trial proposed; no geometry applied', fieldIndex=len(fields) - 1, translation_m=d, analyticProof=proof, sampledGuideMinimumJacobian=float(np.linalg.det(jac).min()), remainingFullSectionShiftResidual_m=required-d)
        correction.append(item)
    if fields:
        from geometry_tools.algorithms.compact_positive_fields import compose_compact_radial_maps
        trial, jacobian, proof = compose_compact_radial_maps(guide_points, fields)
        distance, distance_evidence = skin_tree.closest_points(trial, backend='bvh')
        guide_trial = {"pointsTested": len(trial), "minimumActualSkinDistance_m": float(distance['distance_m'].min()), "minimumGuideJacobian": float(np.linalg.det(jacobian).min()), "analyticProof": proof, "distanceEvidence": distance_evidence, "qualifications": 'Finite source cross-section points only, unsigned actual triangle distances. No full vessel-surface clearance, inside, patency or clinical margin acceptance.'}
        np.savez_compressed(out / 'popliteal-section-guides.npz', actual_source_section_points_m=guide_points, source_vessel_name=np.asarray(guide_meta), trial_mapped_section_points_m=trial, actual_skin_nearest_points_m=distance['points'], actual_skin_face=distance['face_ids'], actual_skin_distance_m=distance['distance_m'])
    else:
        guide_trial = {"status": "No numerical trial field warranted by measured bounds"}
        np.savez_compressed(out / 'popliteal-section-guides.npz', actual_source_section_points_m=guide_points, source_vessel_name=np.asarray(guide_meta))
    by_vessel = Counter((q['name'], q['bonePart']) for q in strict_rows)
    recipe = {"schema": "actual-lower-limb-joint-constraints/v1", "baselineCandidate": str(candidate), "inputBindings": bindings, "identityFrameBinding": {"candidate": 'model:local-ta98:assembled-metres', "fixed": frame, "evidence": 'Exact target catalogue and target skin hash binding in candidate contract; identity coordinates, no transform inferred from bounds'}, "fixedBoneParts": [{k: r[k] for k in ['partIndex', 'id', 'name', 'faceCount']} for r in bone_rows], "fixedBoneFaceCount": len(bone_map), "selectedVascularInstances": len(selected_rows), "sampledVesselFaces": len(triangles), "sampling": 'Main named femoral/popliteal/tibial/fibular trunks use all retained faces; genicular branch owners up to96 evenly spaced source faces; other local routes up to32. No whole-body scan.', "strictQueryCompleteForSelectedFacesWithin1mm": query_complete, "candidatePairsTested": tested, "relations": dict(Counter(result['relation'].tolist())), "strictLiteralCrossingPairs": len(strict_rows), "strictByVesselBone": [{"vessel": k[0], "bonePart": k[1], "sampledStrictPairs": n} for k, n in by_vessel.most_common()], "strictWitnessOwnership": strict_rows, "sourceClosure": str(out / 'source-closure.json'), "nutrientForaminalPolicy": 'Nutrient named owners are possible intended opening routes; accept only an actual source-bound opening/canal witness. Perforating, acetabular, genicular or recurrent labels alone never authorize crossing cortical triangles. No such opening was inferred from a centroid.', "poplitealActualSectionGuides": guide_rows, "firstCorrection": {"priority": 'Knee family', "basis": 'Measured popliteal surface/plane section points and actual posterior fixed-bone/skin triangle hits, not bone centroids or claimed foramina', "axisConvention": 'Y superior; negative Z posterior, inherited active atlas frame', "rule": 'One half of measured posterior bone-corridor residual, capped at8mm, further limited by actual posterior-skin ray margin; no parameter search', "boneMargin_m": .001, "boneMarginQualification": 'Operational collision trial margin only; no anatomical attachment gap is asserted', "skinZMargin_m": .002, "skinZMarginQualification": 'Directional Z separation only; actual closest triangle distance also measured at section points; no full-surface clearance claim', "sideMeasurements": correction, "fields": fields, "application": 'Compose these common fields AFTER the existing Gnew and body-envelope-v2 map for every source-v2 instance. Do not move any bone/skin or displace an isolated trunk. Preserve exact source vertex/bank ancestry and normal transport.', "trialGuideChecks": guide_trial}, "requiredFollowups": ['Apply only once as a named uninstalled candidate if the proposed fields are useful; no iterative centroid fitting', 'Replay identical joint strict triangle samples, major-trunk surface strain and actual skin/ray clearances before accepting the trial', 'Verify femoral/adductor and anterior-tibial interosseous route constraints and actual source connections before treating a semantic subtree as a communicating vessel closure', 'Hip and ankle evidence in this stage remains a measured constraint inventory, not a fitted route'], "geometryModified": False, "installed": False, "anatomicalAcceptance": False, "elapsedSeconds": time.monotonic() - began}
    (out / 'constraints.json').write_text(json.dumps(recipe, indent=2) + '\n')
    assert all(sha(p) == value for p, value in bindings.items()), 'Bound inputs changed during stage'
    checkpoint.write_text(json.dumps({"status": 'complete focused constraint stage', "elapsedSeconds": time.monotonic() - began, "geometryModified": False}) + '\n')
    print(json.dumps({"sampleFaces": len(triangles), "boneFaces": len(bone_map), "strictPairs": len(strict_rows), "completeSampleQuery": query_complete, "topStrict": recipe['strictByVesselBone'][:12], "firstCorrection": correction, "guideTrial": guide_trial, "elapsedSeconds": time.monotonic() - began}))


if __name__ == '__main__':
    main()
