"""Fail-closed input contracts for a coupled, local male tongue proposal.

This module neither packs geometry nor installs a catalogue.  In particular,
closed covering cells do not authorize coating muscle-entry boundaries.  Their
complete actual covering/muscle audit and a separate exposed-surface review are
required before a future packer may emit those cells.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np

from geometry_tools.common import digest_file
from geometry_tools.model import Mesh

FRAME = 'model:local-ta98:assembled-metres'
MODE = 'macroscopic-muscle-extent-reference-envelope'
MUCO = 'TA98TONGUE:male:mucosa'
APO = 'TA98TONGUE:male:aponeurosis'
INTRINSIC_KEYS = frozenset(('superior-longitudinal', 'inferior-longitudinal-left',
                            'inferior-longitudinal-right', 'transverse-vertical-shared'))
RETAINED_KEYS = frozenset(('ZA:Genioglossus muscle.l', 'ZA:Genioglossus muscle.r',
                           'ZA:Hyoglossus muscle.l', 'ZA:Hyoglossus muscle.r',
                           'styloglossus-v2-left', 'styloglossus-v2-right'))
EXTRINSIC_IDS = frozenset(pid for pid in RETAINED_KEYS if pid.startswith('ZA:'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_binding(binding, local_root):
    """Read bindings by explicit value, never infer the latest candidate."""
    require(isinstance(binding, dict), 'An explicit path/hash binding is required')
    path = Path(binding.get('path', '')).resolve()
    require(path.is_relative_to(local_root.resolve()) and path.is_file(),
            'Candidate input must be an existing explicit local file')
    require(digest_file(path) == binding.get('sha256'), 'Candidate input changed: '+str(path))
    return path


def selected_extrinsics_gate(manifest, source_sha, local_root):
    """Exactly four selected, actually passed GG/HG candidates, never fallbacks."""
    require(manifest.get('schema') == 'tongue-selected-extrinsic-candidates/v1'
            and manifest.get('passed') is True
            and manifest.get('sourceCatalogueSha256') == source_sha
            and manifest.get('fixedBonesModified') is False,
            'Select an actually passed, source-bound four-muscle insertion/origin candidate')
    rows = manifest.get('parts', [])
    require(len(rows) == 4 and {r.get('originalPartId') for r in rows} == EXTRINSIC_IDS,
            'Explicitly select both GG and both HG; old native HG is not an implicit fallback')
    result = {}
    for row in rows:
        path = checked_binding(dict(path=row['mesh'],sha256=row['sha256']),local_root)
        mesh = Mesh.load(path)
        require(mesh.frame == FRAME and mesh.geometry_sha == row['geometrySha256']
                and row.get('geometryOrigin') == [0.,1.625,0.],
                'Selected extrinsic geometry/frame or reviewed encoding origin changed')
        report_path = checked_binding(row.get('attachmentReport'),local_root)
        report = json.loads(report_path.read_text())
        matches = []
        def walk(value):
            if isinstance(value,dict):
                for name,sha_name in [('candidate','candidateSha256'),('mesh','sha256'),
                                      ('mesh','mesh_sha256'),('candidate','sha256')]:
                    if value.get(name) and value.get(sha_name) == row['sha256']:
                        if Path(value[name]).resolve() == path:
                            matches.append(value)
                for child in value.values():
                    walk(child)
            elif isinstance(value,list):
                for child in value:
                    walk(child)
        walk(report)
        gg_checked = (report.get('schema') == 'tongue-intralingual-mucosal-fit/v1'
                      and report.get('candidatePassedNumericalSurfaceGate') is True
                      and report.get('fixedBoneGeometry') is True
                      and any(m.get('checks',{}).get('passed') is True for m in matches))
        require(matches and (report.get('passed') is True or gg_checked or any(m.get('passed') is True
                or m.get('numerical_geometry_passed') is True for m in matches)),
                'Actual producer attachment report does not pass and bind this selected mesh')
        require(not report.get('source_bones_modified',False) and not report.get('fixedBonesModified',False),
                'Attachment evidence moved a fixed bone')
        if row.get('storedBank'):
            checked_binding(dict(path=row['storedBank'],sha256=row.get('storedBankSha256')),local_root)
        result[row['originalPartId']] = dict(row,attachmentReportPath=str(report_path),
            attachmentReportSha256=digest_file(report_path))
    return result


def complete_pair(summary, *, allow_coplanar=False):
    require(summary.get('complete_within_pair_budget') is True
            and summary.get('truncated') is False,
            'An actual triangle-pair measurement is absent or incomplete')
    require(summary.get('crossing_face_pairs') == 0,
            'Unresolved covering/tissue or fixed-bone crossings block packing')
    require(allow_coplanar or summary.get('coplanar_overlap_face_pairs') == 0,
            'Unresolved coplanar occupied surfaces block packing')
    require(summary.get('degenerate_pair_count', 0) == 0,
            'A degenerate pair cannot establish an accepted interface')


def covering_muscle_gate(audit, covers, muscles, *, check_files=True, shared_interfaces=None):
    """Require every 2 × 10 actual pair; no muscle-envelope exemption here.

    `covers` and `muscles` map each audit key to its actual mesh path/hash.
    This adapter accepts the existing measured audit schema, including a failed
    audit, and rejects unresolved intersections without normalizing them away.
    """
    require(audit.get('schema') == 'tongue-finite-covering-muscle-audit/v1',
            'Unknown covering/muscle audit schema')
    require(set(covers) == {'mucosa', 'aponeurosis'} and len(muscles) == 10,
            'Both finite coverings and all ten muscle counterparts are required')
    shared_interfaces = shared_interfaces or {}
    expected = {(cover, muscle) for cover in covers for muscle in muscles}
    pairs = audit.get('pairs', [])
    actual = {(p.get('cover'), p.get('muscle')) for p in pairs}
    require(len(pairs) == len(expected) and actual == expected,
            'Covering/muscle audit must measure all twenty identity-distinct pairs')
    failures = []
    for row in pairs:
        for role, record in [('cover', covers[row['cover']]),
                             ('muscle', muscles[row['muscle']])]:
            require(Path(row[role+'_path']).resolve() == Path(record['path']).resolve()
                    and row[role+'_sha256'] == record['sha256'],
                    'Covering audit binds another physical source: '+role)
            if check_files:
                require(digest_file(record['path']) == record['sha256'],
                        'Covering audit source has changed')
        witness = row.get('witnesses')
        require(witness and row.get('witnesses_sha256'),
                'Complete audit must retain its actual triangle witnesses')
        if check_files:
            require(digest_file(witness) == row['witnesses_sha256'],
                    'Covering audit triangle witnesses changed')
        try:
            seam = shared_interfaces.get((row['cover'], row['muscle']))
            if seam:
                proof = seam['actual_check']
                require(proof.get('passed') is True
                        and proof.get('literal_opposed_triangles') is True
                        and proof.get('literal_opposed_normal_words') is True,
                        'A declared seam needs actual literal opposed facet/normal proof')
                coplanar = row.get('coplanar_face_pairs', [])
                require(len(coplanar) == row.get('summary', {}).get('coplanar_overlap_face_pairs'),
                        'Every allowed coplanar pair needs explicit source-face witnesses')
                require(all(a in seam['cover_faces'] and b in seam['muscle_faces']
                            for a, b in coplanar),
                        'Coplanar occupied surfaces extend outside the literal declared seam')
            complete_pair(row.get('summary', {}), allow_coplanar=bool(seam))
        except ValueError as error:
            failures.append(f"{row['cover']} / {row['muscle']}: {error}")
    require(not failures, '; '.join(failures))
    require(audit.get('passed') is True, 'Covering/muscle audit did not pass')
    return dict(pairs=len(pairs), passed=True, muscle_envelope_exception_used=False,
                literal_shared_interfaces=len(shared_interfaces))


def retained_envelope_gate(audit, intrinsics, retained):
    """Allow only measured, declared muscle-envelope intersections.

    Repaired GG hashes must be supplied by the selected GG candidate.  Old
    native GG evidence therefore cannot silently pass this coupled contract.
    Coverings, the whole-organ exterior and bones are never exempted.
    """
    require(set(intrinsics) == INTRINSIC_KEYS and set(retained) == RETAINED_KEYS,
            'Four actual ROIs and six explicitly retained muscle sources are required')
    require(audit.get('representation_mode') == MODE
            and audit.get('tissue_ownership_claim') is False
            and audit.get('roi_boundaries') == 'approximate',
            'Muscle-envelope diagnostics must explicitly disclaim tissue ownership')
    expected = {(a, b) for a in intrinsics for b in retained}
    pairs = audit.get('pairs', [])
    require(len(pairs) == len(expected) == 24
            and {(p.get('intrinsic'), p.get('retained')) for p in pairs} == expected,
            'All twenty-four actual ROI/retained-muscle measurements are required')
    for row in pairs:
        first, second = intrinsics[row['intrinsic']], retained[row['retained']]
        require(row.get('intrinsic_geometry_sha256') == first['geometry_sha256']
                and row.get('retained_geometry_sha256') == second['geometry_sha256']
                and row.get('intrinsic_sha256') == first['sha256']
                and row.get('retained_sha256') == second['sha256'],
                'Macro-envelope report is stale, including actual repaired GG')
        measured = row.get('summary', row)
        require(measured.get('complete_within_pair_budget') is True
                and measured.get('truncated') is False,
                'Incomplete muscle-envelope measurement')
        require(all(isinstance(measured.get(k), int) and measured[k] >= 0
                    for k in ('crossing_face_pairs', 'coplanar_overlap_face_pairs')),
                'Missing actual muscle-envelope crossing counts')
        require(row.get('primary_reference_url', audit.get('primary_reference_url'))
                == 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
                'Interwoven-envelope allowance requires the selected primary reference')
    require(audit.get('complete_within_pair_budget') is True
            and audit.get('passed') is True,
            'Retained muscle measurements did not complete')
    return dict(pairs=24, passed=True, tissue_ownership_claim=False)


def corner_normal_mesh(path, expected_geometry_sha):
    """Expand storage charts only; preserve every literal physical triangle."""
    mesh = Mesh.load(path)
    require(mesh.frame == FRAME and mesh.geometry_sha == expected_geometry_sha,
            'Intrinsic physical geometry or frame changed')
    with np.load(path, allow_pickle=False) as raw:
        require('corner_normals' in raw, 'Selected intrinsic corner-normal words are absent')
        words = np.asarray(raw['corner_normals'])
    require(words.shape == (*mesh.faces.shape, 3)
            and words.dtype == np.dtype('int16')
            and np.all(np.linalg.norm(words.astype(float), axis=-1) > 0),
            'Intrinsic corner-normal words must be nonzero actual Int16 fields')
    keys = np.column_stack((mesh.faces.reshape(-1), words.reshape(-1, 3))).astype(np.int64)
    chart, inverse = np.unique(keys, axis=0, return_inverse=True)
    attributes={key:np.asarray(value)[chart[:,0]].copy()for key,value in mesh.attributes.items()
                if key not in ('normals','normalSigns')}
    attributes.update(normals=chart[:,1:].astype('<i2'),normalSigns=np.ones(len(chart),dtype='i1'))
    expanded = Mesh(mesh.vertices[chart[:, 0]], inverse.reshape(-1, 3), mesh.frame,
                    attributes, copy.deepcopy(mesh.source))
    require(np.array_equal(expanded.triangles, mesh.triangles)
            and expanded.geometry_sha == mesh.geometry_sha,
            'A normal chart may not move vertices or change physical face order')
    require(np.array_equal(expanded.attributes['normals'][expanded.faces], words),
            'Expanded storage did not preserve the actual selected corner words')
    normal = np.cross(np.diff(mesh.triangles[:, :2], axis=1)[:, 0],
                      mesh.triangles[:, 2]-mesh.triangles[:, 0])
    require(np.all(np.einsum('fci,fi->fc', words.astype(float), normal) > 0),
            'A selected corner normal points through its physical triangle')
    return expanded


def points_on_triangles(points, triangles, tolerance=2e-8):
    """Small bounded reference-control check; not an anatomical classifier."""
    triangles = np.asarray(triangles, float)
    require(len(triangles) > 0 and len(triangles) <= 30000,
            'Feature ancestry triangle budget exceeded or empty')
    a, b = triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]
    aa, ab, bb = np.einsum('ij,ij->i', a, a), np.einsum('ij,ij->i', a, b), np.einsum('ij,ij->i', b, b)
    den = aa*bb-ab*ab
    require(np.all(den > 0), 'Feature carrier has a degenerate ancestor triangle')
    normal = np.cross(a, b)
    lengths = np.linalg.norm(normal, axis=1)
    scale = np.maximum(np.minimum(np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)), 1e-12)
    for point in points:
        offset = np.asarray(point)-triangles[:, 0]
        au, bu = np.einsum('ij,ij->i', offset, a), np.einsum('ij,ij->i', offset, b)
        u, v = (au*bb-bu*ab)/den, (bu*aa-au*ab)/den
        eps = tolerance/scale
        distance = np.abs(np.einsum('ij,ij->i', offset, normal))/lengths
        require(np.any((distance <= tolerance) & (u >= -eps) & (v >= -eps)
                       & (u+v <= 1+eps)), 'Curve point is outside actual exposed source-face ancestry')


def bind_features(features, ancestry, source, carrier, bank_sha):
    """Map source features to selected exposed carrier faces; never drop facets.

    Ancestry rows are {sourceFace, carrierFace, cornerOrder:[0,1,2]}.  The
    corner order maps each carrier corner to that source triangle's corner.
    A later open-exposure *finite cell* may use any subset of source faces.
    Features which refer to a removed entry-boundary face must be regenerated,
    rather than being silently clipped or redirected to another tissue.
    """
    require(features.get('schema') == 'tongue-construction-features/v1'
            and features.get('hostPartId') == MUCO
            and features.get('additionalSolidMeshes') == 0,
            'Use typed source-bound features without extra solids')
    require(ancestry.get('schema') == 'tongue-covering-face-ancestry/v1'
            and ancestry.get('carrierPartId') == MUCO
            and ancestry.get('carrierBankSha256') == bank_sha
            and ancestry.get('sourceGeometrySha256') == source.geometry_sha
            and ancestry.get('carrierGeometrySha256') == carrier.geometry_sha,
            'Actual emitted carrier/source ancestry is absent or stale')
    require(source.frame == carrier.frame == FRAME, 'Feature coordinate frame mismatch')
    mapping, carrier_faces = {}, set()
    carrier_space = features.get('featureFaceSpace') == 'carrier'
    require(features.get('featureFaceSpace', 'source') in ('source', 'carrier'),
            'Unknown typed feature face-index space')
    if carrier_space:
        require(features.get('exposureRegenerated') is True
                and features.get('carrierGeometrySha256') == carrier.geometry_sha,
                'Carrier-index features must be explicitly regenerated on actual exposed tissue')
    for row in ancestry.get('faces', []):
        sf, cf, order = row['sourceFace'], row['carrierFace'], row.get('cornerOrder')
        require(isinstance(sf, int) and isinstance(cf, int)
                and 0 <= sf < len(source.faces) and 0 <= cf < len(carrier.faces),
                'Source ancestry face is outside its actual mesh')
        require(cf not in carrier_faces, 'Ambiguous duplicate carrier face ancestry')
        if order is not None:
            require(order in ([0, 1, 2], [1, 2, 0], [2, 0, 1])
                    and np.array_equal(carrier.triangles[cf], source.triangles[sf][order]),
                    'A literal exposed triangle changed coordinates or outward orientation')
        else:
            bary = np.asarray(row.get('sourceBarycentric'), float)
            require(bary.shape == (3, 3) and np.isfinite(bary).all()
                    and np.all(bary >= 0) and np.all(bary <= 1)
                    and np.all(np.abs(bary.sum(1)-1) <= 1e-12),
                    'Cut facets need actual three-corner source barycentric provenance')
            reconstructed = bary @ source.triangles[sf]
            require(np.linalg.norm(reconstructed-carrier.triangles[cf], axis=1).max() <= 2e-8,
                    'Cut carrier facet differs from actual source ancestry beyond storage tolerance')
            old_normal = np.cross(source.triangles[sf, 1]-source.triangles[sf, 0],
                                  source.triangles[sf, 2]-source.triangles[sf, 0])
            new_normal = np.cross(carrier.triangles[cf, 1]-carrier.triangles[cf, 0],
                                  carrier.triangles[cf, 2]-carrier.triangles[cf, 0])
            require(float(old_normal@new_normal) > 0, 'Cut carrier facet reversed source orientation')
            require(carrier_space, 'Cut subdivisions require explicitly regenerated carrier-index controls')
        mapping.setdefault(sf, []).append(row)
        carrier_faces.add(cf)
    require(mapping, 'At least one actual exposed source facet is required')
    result = copy.deepcopy(features['termGroups'])
    require(len(result) == features.get('terms') == 6, 'Expected six typed controls')
    for group in result:
        require(group.get('partIds') == [MUCO], 'Feature aliases must count only the actual carrier')
        if group.get('surface'):
            surface = group['surface']
            require(surface.get('partId') == MUCO, 'Wrong surface carrier identity')
            require(all(face in carrier_faces if carrier_space else face in mapping
                        for face in surface['faces']),
                    'Feature includes a removed mucosal entry-boundary facet; regenerate it')
            if not carrier_space:
                require(all(len(mapping[f]) == 1 and 'cornerOrder' in mapping[f]
                            for f in surface['faces']),
                        'Partial exposure needs deliberate feature regeneration, not silent clipping')
                surface['faces'] = [mapping[f][0]['carrierFace'] for f in surface['faces']]
            surface.pop('sourceMeshSha256', None)
            surface['sourceBankSha256'] = bank_sha
        if group.get('anchor'):
            anchor = group['anchor']
            require(anchor.get('partId') == MUCO
                    and (anchor['face'] in carrier_faces if carrier_space else anchor['face'] in mapping),
                    'Anchor has no exposed mucosal source ancestor')
            weights = np.asarray(anchor['barycentric'], float)
            require(weights.shape == (3,) and np.isfinite(weights).all()
                    and np.all(weights >= 0) and abs(float(weights.sum())-1) <= 1e-12,
                    'Invalid original anchor barycentric coordinates')
            if not carrier_space:
                require(len(mapping[anchor['face']]) == 1
                        and 'cornerOrder' in mapping[anchor['face']][0],
                        'A cut anchor needs deliberate regeneration on its actual exposed facet')
                row = mapping[anchor['face']][0]
                weights = weights[row['cornerOrder']]
                anchor['face'] = row['carrierFace']
            anchor['barycentric'] = weights.tolist()
            actual = weights @ carrier.triangles[anchor['face']]
            require(np.linalg.norm(actual-np.asarray(anchor['point'])) <= 2e-8,
                    'Feature anchor differs from actual emitted carrier')
            anchor['point'] = actual.tolist()
        for curve in group.get('curves', []):
            require(curve.get('partId') == MUCO, 'Wrong curve carrier identity')
            points_on_triangles(curve['points'], carrier.triangles[sorted(carrier_faces)])
            curve.pop('sourceMeshSha256', None)
            curve['sourceBankSha256'] = bank_sha
    return result
