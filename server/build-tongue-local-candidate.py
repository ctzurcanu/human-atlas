"""Pack an explicitly selected, numerically checked local male tongue proposal.

Never chooses the latest construction, installs anatomy or edits source files.
Run only after inspecting the selected intrinsic candidate's renders/sections.
The separate gated install job records that visual selection.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time
from urllib.parse import urlsplit

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from geometry_tools.constraints import check_constraint
from geometry_tools.contact import pair_contacts

WORK = ROOT / '.local-models/geometry-workbench/tongue-completion'
FRAME = 'model:local-ta98:assembled-metres'
ORIGIN = np.array([0., 1.5, 0.])
WHOLE = 'TA98TONGUE:male:whole'
MUSCLES = 'TA98TONGUE:male:muscles'
LEGACY_STYLO = ['LOCAL:male:l_styloglossus:1030', 'LOCAL:male:r_styloglossus:981']
NATIVE_MUSCLES = {
    'genioglossus': ('A05.1.04.101', ['ZA:Genioglossus muscle.l', 'ZA:Genioglossus muscle.r']),
    'hyoglossus': ('A05.1.04.102', ['ZA:Hyoglossus muscle.l', 'ZA:Hyoglossus muscle.r']),
}
INTRINSIC_TERMS = {'superior-longitudinal':'A05.1.04.106', 'inferior-longitudinal':'A05.1.04.107', 'transverse':'A05.1.04.108', 'vertical':'A05.1.04.109'}
NAMES = {'genioglossus':'Genioglossus', 'hyoglossus':'Hyoglossus', 'styloglossus':'Styloglossus', 'superior-longitudinal':'Superior longitudinal muscle', 'inferior-longitudinal':'Inferior longitudinal muscle', 'transverse':'Transverse muscle of tongue', 'vertical':'Vertical muscle of tongue'}
ENVELOPE_MODE = 'macroscopic-muscle-extent-reference-envelope'
STRICT_MODE = 'finite-exclusive-tissue'
ENVELOPE_COUNTERPARTS = {p for _, ids in NATIVE_MUSCLES.values() for p in ids} | {'styloglossus-v2-left', 'styloglossus-v2-right'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def checked_file(path, checksum):
    path = Path(path).resolve()
    require(path.is_relative_to(ROOT / '.local-models'), 'Geometry input must stay local')
    require(path.is_file() and digest_file(path) == checksum, f'Stale input: {path}')
    return path


def normal_words(mesh):
    normals = np.asarray(mesh.attributes['normals'])
    require(normals.shape == mesh.vertices.shape and np.isfinite(normals).all(), 'Invalid original normals')
    signs = np.asarray(mesh.attributes.get('normalSigns', np.ones(len(normals), dtype='i1')))
    require(signs.shape == (len(normals),) and np.isin(signs, [-1, 1]).all(), 'Invalid original normal signs')
    if np.issubdtype(normals.dtype, np.integer):
        require(np.abs(normals.astype(np.int64)).max() <= 32767, 'Normal words outside signed Int16 range')
        return normals.astype('<i2'), signs.astype('i1')
    lengths = np.linalg.norm(normals, axis=1)
    require(np.all(lengths > 0), 'Zero original normal')
    return np.rint(normals / lengths[:, None] * 32767).clip(-32767, 32767).astype('<i2'), signs.astype('i1')


def opposed_origin_chart(mesh, bone, origin):
    """Separate cap corners retain bone words/signs without a wall-normal conflict."""
    require(mesh.frame == bone.frame == FRAME, 'Attachment frame mismatch')
    source_faces = origin['source_faces']
    require(len(source_faces) == 3, 'Expected the reviewed three-face origin')
    words, signs = normal_words(mesh)
    bone_words, bone_signs = normal_words(bone)
    vertices, packed_words, packed_signs = list(mesh.vertices.copy()), list(words.copy()), list(signs.copy())
    faces = mesh.faces.copy()
    available = {tuple(sorted(map(tuple, bone.triangles[f]))): int(f) for f in source_faces}
    corners = {}
    for face_id in range(3):
        triangle = mesh.triangles[face_id]
        key = tuple(sorted(map(tuple, triangle)))
        require(key in available, 'Reviewed origin triangle is not literal native cortex')
        bone_face = bone.faces[available[key]]
        for corner, original_vertex in enumerate(mesh.faces[face_id]):
            xyz = mesh.vertices[original_vertex]
            matches = [int(v) for v in bone_face if np.array_equal(bone.vertices[v], xyz)]
            require(len(matches) == 1, 'Ambiguous native origin corner')
            native_vertex = matches[0]
            seam_key = (tuple(xyz), tuple(bone_words[native_vertex]), int(-bone_signs[native_vertex]))
            if seam_key not in corners:
                corners[seam_key] = len(vertices)
                vertices.append(xyz.copy())
                packed_words.append(bone_words[native_vertex].copy())
                packed_signs.append(-bone_signs[native_vertex])
            faces[face_id, corner] = corners[seam_key]
    result = Mesh(np.asarray(vertices), faces, FRAME, {'normals':np.asarray(packed_words, dtype='<i2'), 'normalSigns':np.asarray(packed_signs, dtype='i1')}, copy.deepcopy(mesh.source))
    require(result.geometry_sha == mesh.geometry_sha, 'Normal chart changed the physical styloglossus surface')
    spec = dict(kind='shared-interface', parameters=dict(first_faces=source_faces, second_faces=list(range(3))))
    proof = check_constraint(spec, bone, result)
    require(proof['passed'], 'Origin chart lacks literal opposed coordinates and normal words')
    return result, spec, dict(**proof, normal_corner_vertices_added=len(corners), wall_normal_words_retained=True)


def add_array(buffer, values):
    buffer.extend(b'\0' * ((-len(buffer)) % 4))
    offset = len(buffer)
    buffer.extend(values.tobytes())
    return offset


def intrinsic_terms(row):
    """One physical carrier may explicitly represent interwoven T/V territory."""
    terms = row.get('terms', [row['term']] if 'term' in row else [])
    require(isinstance(terms, list) and terms and len(set(terms)) == len(terms), 'Intrinsic row needs unique explicit term/terms metadata')
    require(set(terms) <= set(INTRINSIC_TERMS.values()), 'Unsupported intrinsic identity')
    if len(terms) > 1:
        require(set(terms) == {'A05.1.04.108', 'A05.1.04.109'}, 'Only the explicitly interwoven T/V territory may share a carrier')
        require(bool(row.get('coverageLimitation')), 'Shared fibre territory needs an explicit anatomical limitation')
    return terms


def normal_presentation_gate(path, construction, metadata):
    if not construction.get('presentation_normal_variant_passed'):
        return None
    binding, base_binding = metadata.get('presentation_normal_proof', {}), metadata.get('base_construction', {})
    proof_path = checked_file(binding.get('path', ''), binding.get('sha256'))
    base_path = checked_file(base_binding.get('path', ''), base_binding.get('sha256'))
    require(proof_path == path.with_name('normal-only-proof.json') and digest_file(proof_path) == construction.get('presentation_normal_proof_sha256'), 'Selected presentation-normal proof is stale')
    proof, base = read(proof_path), read(base_path)
    require(proof.get('schema') == 'tongue-normal-only-presentation-variant/v1' and proof.get('passed') is True and proof.get('geometry_unchanged') is True and proof.get('geometry_defects_repaired_or_hidden') is False, 'Normal-only changes cannot authorize a geometry change')
    require(base_binding['sha256'] == construction.get('base_construction_sha256') == proof.get('base_construction_sha256') and proof.get('source_catalogue_sha256') == construction.get('source_catalogue_sha256'), 'Normal variant base/source binding changed')
    source = read(WORK/'source-v1/source-audit.json')
    require(proof.get('all_original_source_normal_words_unchanged') is True and proof.get('styloglossus_bone_caps_untouched') is True and proof.get('source_tongue_mesh_sha256') == source['contextParts'][0]['meshSaveSha256'], 'Native normal/cap preservation evidence is missing')
    base_rows, proof_rows = {r['key']:r for r in base['meshes']}, {r['key']:r for r in proof['meshes']}
    require(set(base_rows) == set(proof_rows) == {r['key'] for r in construction['meshes']}, 'Normal variant carrier identities changed')
    records = []
    for row in construction['meshes']:
        old, record = base_rows[row['key']], proof_rows[row['key']]
        current_path, old_path = checked_file(row['mesh'], row['sha256']), checked_file(old['mesh'], old['sha256'])
        current, original = Mesh.load(current_path), Mesh.load(old_path)
        require(current.geometry_sha == original.geometry_sha == row['geometry_sha256'] == old['geometry_sha256'] and np.array_equal(current.vertices, original.vertices) and np.array_equal(current.faces, original.faces), 'Normal variant changed physical vertices or triangle order')
        require(record['variant_mesh_sha256'] == row['sha256'] and record['base_mesh_sha256'] == old['sha256'] and record['variant_geometry_sha256'] == current.geometry_sha, 'Normal proof does not bind the actual selected words')
        require(current.frame == original.frame and current.attributes.keys() == original.attributes.keys() and all(np.array_equal(value, original.attributes[key]) for key, value in current.attributes.items() if key != 'normals'), 'Normal variant changed a non-normal attribute or frame')
        words, signs = normal_words(current)
        effective = words.astype(float) * signs[:, None]
        lengths = np.linalg.norm(effective, axis=1)
        require(np.all(lengths > 0), 'Normal-only variant contains zero words')
        effective /= lengths[:, None]
        tri = current.triangles
        geometric = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        geometric /= np.linalg.norm(geometric, axis=1)[:, None]
        cosine = np.einsum('fci,fi->fc', effective[current.faces], geometric)
        require(np.isfinite(cosine).all() and np.all(cosine > 0), 'Actual selected Int16 normal words point through a face')
        records.append(dict(key=row['key'], mesh_path=str(current_path), mesh_sha256=row['sha256'], base_mesh_path=str(old_path), base_mesh_sha256=old['sha256'], geometry_sha256=current.geometry_sha, exact_vertices_and_indices=True, actual_minimum_face_corner_cosine=float(cosine.min()), actual_nonpositive_face_corners=0))
    return dict(path=str(proof_path), sha256=digest_file(proof_path), base_construction=dict(path=str(base_path), sha256=digest_file(base_path)), records=records, geometry_unchanged=True, anatomical_acceptance=False)


def packing_construction(path):
    """Normalize a source-bound envelope manifest without changing its evidence."""
    construction = read(path)
    if construction.get('schema') != 'tongue-macroscopic-muscle-reference-envelopes/v1':
        return construction, None
    metadata_path = path.with_name('packing-input.json')
    metadata = read(metadata_path)
    require(metadata.get('schema') == 'tongue-macro-envelope-packing-input/v1', 'Unknown envelope packing metadata schema')
    require(metadata.get('construction_sha256') == digest_file(path), 'Envelope packing metadata is stale')
    require(metadata.get('source_audit_sha256') == digest_file(WORK/'source-v1/source-audit.json'), 'Envelope source-audit binding is stale')
    require(metadata.get('source_catalogue_sha256') == construction.get('source_catalogue_sha256') and metadata.get('storage_origin') == ORIGIN.tolist(), 'Envelope packing source or local origin changed')
    require(metadata.get('checked_geometry_sha256') == digest_file(path.with_name('checked-geometry.json')) and metadata.get('retained_context_summary_sha256') == digest_file(path.parent/'retained-context-audit/summary.json'), 'Envelope numerical/context sidecar binding is stale')
    require(construction.get('exclusive_tissue_partition_claimed') is False and construction.get('anatomical_acceptance') is False and construction.get('source_geometry_modified') is False, 'Envelope construction must explicitly disclaim tissue partitioning and anatomical acceptance')
    for key in ('source_bank_sha256', 'source_tongue_geometry_sha256', 'all_canonical_muscle_terms', 'representation_mode', 'tissue_ownership_claim', 'roi_boundaries', 'retained_context_policy'):
        require(key in metadata, 'Envelope packing metadata is missing '+key)
        construction[key] = copy.deepcopy(metadata[key])
    require(construction['representation_mode'] == ENVELOPE_MODE and construction['tissue_ownership_claim'] is False, 'Envelope metadata cannot promote this construction to exclusive tissue')
    require(construction['retained_context_policy']['primary_reference_url'] == construction.get('primary_reference'), 'Envelope anatomical reference changed during normalization')
    overrides = metadata.get('mesh_overrides', metadata.get('intrinsic_mesh_overlays', {}))
    for row in construction['meshes']:
        require(row.get('representation') == ENVELOPE_MODE and row.get('exclusive_tissue_volume') is False and row.get('anatomical_acceptance') is False, 'Every envelope carrier must explicitly retain its reference-only identity')
        override = overrides.get(row['key'], {})
        row['partId'] = override.get('partId', 'TA98TONGUE:male:intrinsic:'+row['key'])
        row['coverageLimitation'] = override.get('coverageLimitation', 'Approximate macroscopic muscle-extent reference envelope, not exclusive histological tissue. Intersections with explicitly reported GG/HG/styloglossus muscle envelopes are not hidden. Exact individual anatomy and ROI boundaries remain unverified.')
    normal_proof = normal_presentation_gate(path, construction, metadata)
    return construction, dict(path=str(metadata_path), sha256=digest_file(metadata_path), presentation_normals=normal_proof)


def representation_mode(construction, checked):
    mode = construction.get('representation_mode', STRICT_MODE)
    require(mode in (STRICT_MODE, ENVELOPE_MODE), 'Unknown tongue representation mode')
    if mode == ENVELOPE_MODE:
        require(checked.get('representation_mode') == mode and checked.get('tissue_ownership_claim') is False and construction.get('tissue_ownership_claim') is False, 'Envelope mode requires explicit, matching non-ownership evidence')
        require(construction.get('roi_boundaries') == checked.get('roi_boundaries') == 'approximate', 'Envelope ROI boundaries must be marked approximate')
        policy = construction.get('retained_context_policy', {})
        allowed = policy.get('allowed_intersection_part_ids', [])
        require(allowed and len(set(allowed)) == len(allowed) and set(allowed) <= ENVELOPE_COUNTERPARTS, 'Only explicitly whitelisted native GG/HG and styloglossus muscle envelopes may overlap')
        url = urlsplit(policy.get('primary_reference_url', ''))
        require(url.scheme == 'https' and url.hostname and url.path not in ('', '/') and not url.username and not url.password, 'Envelope mode requires an explicit primary anatomical reference URL')
    return mode


def intrinsic_numerical_gate(row, mode):
    if mode != ENVELOPE_MODE:
        return
    cpp, topo, exterior = row.get('cpp_diagnostics', {}), row.get('topology', {}), row.get('exterior_crossings', {})
    require(cpp.get('self_intersections_checked') is True and cpp.get('self_intersection_pairs') == 0, 'Envelope self-intersections were not completely checked')
    require(topo.get('closed_oriented_manifold') is True and topo.get('components') and all(c.get('closed_orientable_genus') == 0 for c in topo['components']), 'Envelope closure/topology did not pass')
    require(exterior.get('complete_within_pair_budget') is True and exterior.get('truncated') is False and exterior.get('crossing_face_pairs') == exterior.get('coplanar_overlap_face_pairs') == 0, 'Envelope has incomplete or failing native-tongue containment crossings')
    query_evidence = row.get('independent_winding_queries', 0) > 0 or row.get('whole_tongue_contained') is True or row.get('checks', row.get('actual_check_flags', {})).get('whole_tongue_contained') is True
    require(query_evidence and row.get('minimum_absolute_winding', 0) > .99, 'Envelope native-tongue winding containment did not pass')


def retained_context_gate(directory, construction_path, checked_path, source_path, source_sha, stylo_path, construction=None, checked=None):
    construction = construction if construction is not None else packing_construction(construction_path)[0]
    checked = checked if checked is not None else read(checked_path)
    mode = representation_mode(construction, checked)
    audit_path = directory / 'summary.json'
    if not directory.exists():
        require(mode != ENVELOPE_MODE, 'Envelope mode requires a complete retained-context crossing report')
        return None
    require(audit_path.is_file(), 'Retained-context audit is incomplete')
    audit = read(audit_path)
    require(audit.get('schema') == 'tongue-intrinsic-retained-context-audit/v1', 'Unknown retained-context audit schema')
    require(audit.get('passed', audit.get('assembly_exclusive_tissue_verified')) is True, 'Retained-context checks did not pass')
    require(audit.get('complete_within_pair_budget') is True and audit.get('pair_count') == audit.get('expected_pair_count') and audit.get('pair_count', 0) > 0, 'Retained-context audit did not finish every pair')
    if mode == STRICT_MODE:
        require(all(audit.get(k) == 0 for k in ('crossing_face_pairs', 'coplanar_overlap_face_pairs', 'unresolved_ownership_pairs')), 'Retained-context tissue ownership remains unresolved')
    else:
        require(audit.get('representation_mode') == mode and audit.get('tissue_ownership_claim') is False and audit.get('roi_boundaries') == 'approximate', 'Retained-context report does not explicitly describe approximate, non-owning envelopes')
    require(audit.get('construction_sha256') == digest_file(construction_path) and audit.get('candidate_check_sha256') == digest_file(checked_path), 'Retained-context intrinsic evidence is stale')
    require(audit.get('source_catalogue_sha256') == source_sha and audit.get('source_audit_sha256') == digest_file(source_path) and audit.get('styloglossus_candidate_sha256') == digest_file(stylo_path), 'Retained-context native evidence is stale')
    context_rows = audit.get('context', [])
    for row in context_rows:
        checked_file(row['mesh'], row['sha256'])
    classifications = []
    if mode == ENVELOPE_MODE:
        source, stylo = read(source_path), read(stylo_path)
        expected = {r['partId']: (r['meshSaveSha256'], r['geometrySha256']) for r in source['contextParts'] if r['partId'] in ENVELOPE_COUNTERPARTS}
        expected.update({'styloglossus-v2-'+side: (r['sha256'], r['geometry_sha256']) for side, r in stylo['sides'].items()})
        require({r.get('key') for r in context_rows} == ENVELOPE_COUNTERPARTS and len(context_rows) == len(ENVELOPE_COUNTERPARTS), 'Envelope report must include exactly the six retained muscle counterparts')
        by_context = {r['key']:r for r in context_rows}
        require(all(expected[r['key']] == (r['sha256'], r['geometry_sha256']) for r in context_rows), 'Retained muscle report sources changed')
        by_intrinsic = {r['key']:r for r in construction['meshes']}
        pairs = audit.get('pairs', [])
        expected_pairs = {(key, retained) for key in by_intrinsic for retained in ENVELOPE_COUNTERPARTS}
        require(len(pairs) == len(expected_pairs) == audit['pair_count'] and {(p.get('intrinsic'), p.get('retained')) for p in pairs} == expected_pairs, 'Envelope report is missing measured muscle pairs')
        policy = construction['retained_context_policy']
        for pair in pairs:
            require(pair.get('intrinsic_geometry_sha256') == by_intrinsic[pair['intrinsic']]['geometry_sha256'] and pair.get('retained_geometry_sha256') == by_context[pair['retained']]['geometry_sha256'], 'Envelope crossing row geometry binding is stale')
            summary = pair.get('contact_summary', {})
            require(summary.get('complete_within_pair_budget') is True and summary.get('truncated') is False, 'Envelope crossing report contains truncated measurements')
            overlaps = bool(summary.get('crossing_face_pairs') or summary.get('coplanar_overlap_face_pairs') or pair.get('intrinsic_samples_inside_retained_envelope') or pair.get('retained_vertices_inside_intrinsic'))
            require(not overlaps or pair['retained'] in policy['allowed_intersection_part_ids'], 'Overlap with a non-whitelisted muscle envelope')
            classifications.append(dict(intrinsic=pair['intrinsic'], counterpart=pair['retained'], classification='allowed-macroscopic-envelope-overlap' if overlaps else 'no-measured-overlap', measurement_scope='Complete surface-crossing diagnostic; no absent containment measurement is inferred', sampled_containment_measured='intrinsic_samples_inside_retained_envelope' in pair, crossing_face_pairs=summary.get('crossing_face_pairs'), coplanar_overlap_face_pairs=summary.get('coplanar_overlap_face_pairs'), primary_reference_url=policy['primary_reference_url'] if overlaps else None))
        require(sum(p['contact_summary']['crossing_face_pairs'] for p in pairs) == audit['crossing_face_pairs'] and sum(p['contact_summary']['coplanar_overlap_face_pairs'] for p in pairs) == audit['coplanar_overlap_face_pairs'], 'Envelope crossing totals do not match the complete pair report')
    return dict(path=str(audit_path), sha256=digest_file(audit_path), passed=True, representation_mode=mode, tissue_ownership_claim=mode != ENVELOPE_MODE, classifications=classifications)


def bone_crossing_gate(mesh, source_rows):
    results = []
    for row in source_rows.values():
        if row['sourcePart'].get('system') != 'skeletal':
            continue
        path = checked_file(row['meshSavePath'], row['meshSaveSha256'])
        bone = Mesh.load(path)
        require(bone.geometry_sha == row['geometrySha256'], 'Bone check source geometry is stale')
        _, summary = pair_contacts(mesh, bone, 0., max_pairs=50000, max_tests=1000000)
        require(summary.get('complete_within_pair_budget') is True and summary.get('truncated') is False and summary.get('crossing_face_pairs') == summary.get('coplanar_overlap_face_pairs') == 0, 'Intrinsic envelope crosses native bone or the bone audit did not finish')
        results.append(dict(part_id=row['partId'], mesh_sha256=row['meshSaveSha256'], geometry_sha256=bone.geometry_sha, contact_summary=summary))
    require(len(results) >= 3, 'The mandible and bilateral temporal bone checks are required')
    return results


def append_mesh(buffer, mesh, *, part_id, concept_id, name, material, terms, chunk, limitation):
    words, signs = normal_words(mesh)
    relative = (mesh.vertices - ORIGIN).astype('<f4')
    indices = mesh.faces.astype('<u4').reshape(-1)
    decoded = relative.astype(float) + ORIGIN
    max_error = float(np.linalg.norm(decoded - mesh.vertices, axis=1).max())
    require(max_error <= 2e-8, 'Local Float32 decoding exceeds the declared 20 nm storage budget')
    part = dict(id=part_id, conceptId=concept_id, name=name, system='muscular', chunk=chunk,
        material=material, vertexCount=len(relative), indexCount=len(indices),
        positions=add_array(buffer, relative), normals=add_array(buffer, words), normalSigns=add_array(buffer, signs), indices=add_array(buffer, indices),
        geometryOrigin=ORIGIN.tolist(), bounds=[decoded[mesh.faces].reshape(-1, 3).min(0).tolist(), decoded[mesh.faces].reshape(-1, 3).max(0).tolist()],
        groups=['Tongue', 'Tongue muscles'], sectionSolid=True,
        license='CC BY-SA 4.0', provenance=dict(label='Human Atlas tongue reconstruction derived from Z-Anatomy',
            url='https://github.com/Z-Anatomy/Models-of-human-anatomy',
            detail='Derived tongue-region geometry or native styloid contact facets; original source preserved. Procedural materials; no external reference mesh or texture arrays copied.'),
        coverageLimitation=limitation, reviewStatus='local proposal; anatomical acceptance unestablished', reviewSource='local-original-construction')
    if len(terms) == 1:
        part['ta98Term'] = terms[0]
    else:
        part['ta98Terms'] = terms
        part['sharedFibreTerritory'] = True
    return part


def material(color):
    return dict(color=[round(float(c)*255) for c in color[:3]], roughness=.72, pairedNormalSigns=True,
        tissueAppearance=dict(color=list(color[:3]), pattern='fibrous', scale=.0004, relief=.000015))


def main(intrinsics, output):
    began = time.monotonic()
    deadline = began + 50
    def remaining():
        if time.monotonic() >= deadline:
            raise TimeoutError('Tongue packing exceeded its 50 second worker budget')
    construction_path = intrinsics.resolve()
    if construction_path.is_dir():
        construction_path = construction_path / 'construction.json'
    require(construction_path.is_relative_to(WORK) and construction_path.name == 'construction.json', 'Supply an explicit tongue construction.json or its directory')
    output = output.resolve()
    require(output.is_relative_to(WORK) and output != WORK and not output.exists(), 'Use a new local tongue output stage')
    construction, packing_metadata = packing_construction(construction_path)
    checked_path = construction_path.with_name('checked-geometry.json')
    checked = read(checked_path)
    mode = representation_mode(construction, checked)
    require(checked.get('passed') is True and checked.get('cross_muscle_defect_pairs') == 0, 'Selected intrinsic numerical checks did not pass')
    require(checked['construction_sha256'] == digest_file(construction_path), 'Intrinsic numerical check is stale')
    rows = construction['meshes']
    require(rows and len({r['key'] for r in rows}) == len(rows) and len({r['partId'] for r in rows}) == len(rows), 'Intrinsic mesh keys and physical part IDs must be unique')
    represented_terms = set(t for row in rows for t in intrinsic_terms(row))
    require(represented_terms == set(INTRINSIC_TERMS.values()) == set(construction['all_canonical_muscle_terms']), 'All four canonical intrinsic identities must be explicitly classified')
    require(construction.get('source_geometry_modified') is False, 'Construction changed native anatomy')
    source_path = WORK / 'source-v1/source-audit.json'
    source = read(source_path)
    catalogue_path = Path(source['sourceCatalogues']['male']['path'])
    source_sha = digest_file(catalogue_path)
    require(source_sha == source['sourceCatalogues']['male']['sha256'] == construction['source_catalogue_sha256'] == checked['source_catalogue_sha256'], 'Current catalogue and selected construction are not the same source revision')
    original = read(catalogue_path)
    by_id = {p['id']:p for p in original['parts']}
    reader = AtlasReader()
    source_rows = {r['partId']:r for r in source['contextParts']}
    for row in source_rows.values():
        remaining()
        require(digest_file(row['bankPath']) == row['bankSha256'], 'Native context bank is stale')
        require(by_id[row['partId']] == row['sourcePart'], 'Native part metadata changed after source audit')
    require(all(id in original.get('defaultHidden', []) for id in LEGACY_STYLO), 'Misaligned legacy styloglossus must remain hidden')
    inventory_path = ROOT / 'reports/ta98-tongue-inventory.json'
    inventory = read(inventory_path)
    require(inventory['term_count'] == 36, 'The complete 36-term inventory is required')
    stylo_path, context_path = WORK / 'styloglossus-v2/candidate.json', WORK / 'styloglossus-v2/context-checks.json'
    stylo, context = read(stylo_path), read(context_path)
    require(stylo['numerical_candidate_pass'] is True and context['passed'] is True and context['candidateSha256'] == digest_file(stylo_path), 'Styloglossus checks are missing, failed or stale')
    require(stylo['source_catalogue_sha256'] == source_sha and stylo['source_audit_sha256'] == digest_file(source_path) and stylo['inventory_sha256'] == digest_file(inventory_path), 'Styloglossus is based on another source or inventory')
    retained_context = retained_context_gate(construction_path.parent/'retained-context-audit', construction_path, checked_path, source_path, source_sha, stylo_path, construction, checked)
    features_path = WORK / 'surface-controls-v1/surface-controls.json'
    features = read(features_path)
    require(features['sourceCatalogueSha256'] == source_sha and features['inventorySha256'] == digest_file(inventory_path), 'Surface features are stale')
    require(features['hostPartId'] == 'ZA:Tongue' and features['terms'] == len(features['termGroups']) == 6 and features['additionalSolidMeshes'] == 0, 'Expected the six source-bound controls without extra solids')
    native = reader.read_part(original, by_id['ZA:Tongue'], frame=FRAME, catalogue_path=catalogue_path, catalogue_sha=source_sha)
    require(native.source['bank_sha256'] == construction['source_bank_sha256'] == features['carrierBankHashes']['ZA:Tongue'], 'Native tongue carrier hash does not match construction/features')
    require(native.geometry_sha == construction['source_tongue_geometry_sha256'], 'Native tongue geometry changed')
    candidate = copy.deepcopy(original)
    require(not candidate.get('finiteTongueIntegration'), 'Do not replace a previously installed tongue through this preparation script')
    bank = bytearray()
    chunk = len(candidate['chunks'])
    inputs, new_parts, materials, contact_inputs, shared_carriers, bone_checks = [], [], {}, {}, [], {}
    for row in rows:
        remaining()
        path = checked_file(row['mesh'], row['sha256'])
        mesh = Mesh.load(path)
        audit_row = next((r for r in checked['rows'] if r['key'] == row['key']), None)
        require(audit_row and audit_row.get('numerical_checks_passed') is True and audit_row['sha256'] == row['sha256'] and audit_row['geometry_sha256'] == mesh.geometry_sha == row['geometry_sha256'], 'Intrinsic NPZ lacks matching passed numerical evidence')
        intrinsic_numerical_gate(audit_row, mode)
        if mode == ENVELOPE_MODE:
            bone_checks[row['partId']] = bone_crossing_gate(mesh, source_rows)
        terms = intrinsic_terms(row)
        require(mesh.frame == FRAME, 'Intrinsic frame mismatch')
        base = next((key for key, term in INTRINSIC_TERMS.items() if terms == [term]), None)
        concept_id = 'TA98TONGUE:male:muscle:' + base if base else 'TA98TONGUE:male:carrier:' + row['key']
        material_id = 'TA98TONGUE:male:material:' + (base or row['key'])
        materials[material_id] = material(row['color'])
        if set(terms) == {'A05.1.04.108', 'A05.1.04.109'}:
            materials[material_id]['tissueAppearance']['fibres'] = dict(
                directions=[[1., 0., 0.], [0., 1., 0.]], weights=[.5, .5], control='tongue-tv')
        side = ' (left)' if row['key'].endswith('-left') else ' (right)' if row['key'].endswith('-right') else ''
        name = NAMES[base]+side if base else row.get('name', 'Shared transverse/vertical fibre territory')
        part = append_mesh(bank, mesh, part_id=row['partId'], concept_id=concept_id, name=name, material=material_id, terms=terms, chunk=chunk, limitation=row['coverageLimitation'])
        if mode == ENVELOPE_MODE:
            part.update(representationMode=mode, tissueOwnershipClaim=False, roiBoundaries='approximate', primaryReferenceUrl=construction['retained_context_policy']['primary_reference_url'])
        if not base:
            shared_carriers.append(dict(id=concept_id, name=name, elements=[part['id']], coverageLimitation=row['coverageLimitation'], sharedFibreTerritory=True, anatomicalAcceptance=False))
        new_parts.append(part); inputs.append((part, mesh, row, path))
    for side, row in stylo['sides'].items():
        remaining()
        require(side in ('left', 'right') and row['numerical']['passed'] is True and row['physical_origin_patch_faces_literal'] is True, 'Unverified styloglossus side')
        path = checked_file(row['mesh'], row['sha256'])
        source_mesh = Mesh.load(path)
        require(source_mesh.geometry_sha == row['geometry_sha256'], 'Styloglossus NPZ geometry changed')
        bone_id = row['origin']['host_part_id']
        bone = reader.read_part(original, by_id[bone_id], frame=FRAME, catalogue_path=catalogue_path, catalogue_sha=source_sha)
        require(bone.geometry_sha == row['origin']['host_geometry_sha256'], 'Styloid source changed')
        mesh, spec, chart = opposed_origin_chart(source_mesh, bone, row['origin'])
        part_id = 'TA98TONGUE:male:styloglossus:' + side
        material_id = 'TA98TONGUE:male:material:styloglossus'
        materials[material_id] = material([.74,.32,.34])
        limit = source_mesh.source.get('limitation', 'Original finite proposal; exact donor footprint and interdigitation remain unverified.')
        part = append_mesh(bank, mesh, part_id=part_id, concept_id='TA98TONGUE:male:muscle:styloglossus', name=f'Styloglossus ({side})', material=material_id, terms=['A05.1.04.105'], chunk=chunk, limitation=limit)
        new_parts.append(part); inputs.append((part, mesh, row, path)); contact_inputs[part_id] = (bone, spec, chart, bone_id)
    require(len(new_parts) == len(rows)+2 and len({p['id'] for p in new_parts}) == len(new_parts), 'Expected selected physical intrinsic records and two styloglossus records')
    require(all(p['id'] not in by_id for p in new_parts), 'Candidate mesh identity collides with installed anatomy')
    output.mkdir(parents=True)
    bank_path = output / 'geometry.bin'
    bank_path.write_bytes(bank)
    url = '/local-models/' + str(bank_path.relative_to(ROOT / '.local-models'))
    candidate['chunks'].append(dict(url=url, bytes=len(bank)))
    candidate['parts'].extend(new_parts)
    candidate.setdefault('materials', {}).update(materials)
    decoded_checks = []
    for part, mesh, row, path in inputs:
        remaining()
        decoded = reader.read_part(candidate, part, frame=FRAME)
        # Read the actual bank, accounting for compact indexed controls.
        vertex_ids = np.asarray(decoded.source['source_vertex_ids'])
        error = np.linalg.norm(decoded.vertices - mesh.vertices[vertex_ids], axis=1)
        require(float(error.max()) <= 2e-8, 'Emitted-bank decode error exceeds 20 nm')
        topo = topology(decoded)
        require(topo['closed_oriented_manifold'] and all(c['closed_orientable_genus'] == 0 for c in topo['components']), 'Emitted geometry is not closed genus-zero components')
        original_topology = topology(mesh)
        require(len(topo['components']) == len(original_topology['components']), 'Packing changed connected components')
        expected_words, expected_signs = normal_words(mesh)
        require(np.array_equal(decoded.attributes['normals'], expected_words[vertex_ids]) and np.array_equal(decoded.attributes['normalSigns'], expected_signs[vertex_ids]), 'Actual stored normal words/signs changed')
        evidence = dict(part_id=part['id'], input_npz=str(path), input_sha256=digest_file(path), source_geometry_sha256=mesh.geometry_sha, decoded_geometry_sha256=decoded.geometry_sha,
            actual_max_position_error_m=float(error.max()), normal_words_and_signs_literal=True, topology=topo, anatomical_acceptance=False)
        if part['id'] in contact_inputs:
            bone, spec, chart, bone_id = contact_inputs[part['id']]
            proof = check_constraint(spec, bone, decoded)
            require(proof['passed'], 'Emitted styloid interface lost literal opposed coordinates or normal words')
            require(float(error.max()) == 0., 'Styloglossus coordinates were not retained literally')
            evidence['styloid_interface'] = dict(host_part_id=bone_id, constraint=spec, prepared_chart=chart, actual_emitted_check=proof)
        decoded_checks.append(evidence)
    native_ids = ['ZA:Tongue'] + [id for _, ids in NATIVE_MUSCLES.values() for id in ids]
    part_ids = native_ids + [p['id'] for p in new_parts]
    muscle_ids = part_ids[1:]
    groups = [dict(key='tongue', name='Tongue', term='A05.1.04.001', conceptId=WHOLE, partIds=part_ids, kind='exact'),
        dict(key='muscles', name='Muscles of tongue', term='A04.1.05.001', conceptId=MUSCLES, parent='tongue', partIds=muscle_ids, kind='exact')]
    for key, (term, ids) in NATIVE_MUSCLES.items():
        groups.append(dict(key=key, name=NAMES[key], term=term, conceptId='TA98TONGUE:male:muscle:'+key, parent='muscles', partIds=ids, kind='exact', coverageLimitation='Existing native external-muscle envelopes, retained literally. Tissue interdigitation and full attachments remain unverified.'))
    for key, term in {**INTRINSIC_TERMS, 'styloglossus':'A05.1.04.105'}.items():
        matching = [p for p in new_parts if term in p.get('ta98Terms', [p.get('ta98Term')])]
        members = [p['id'] for p in matching]
        groups.append(dict(key=key, name=NAMES[key], term=term, conceptId='TA98TONGUE:male:muscle:'+key, parent='muscles', partIds=members, kind='exact', coverageLimitation=' '.join(dict.fromkeys(p['coverageLimitation'] for p in matching))))
    groups.extend(copy.deepcopy(features['termGroups']))
    key_by_term = {g['term']:g['key'] for g in groups}
    authoritative = {r['authority']['id']:r['authority'] for r in inventory['terms']}
    for g in groups:
        if g['key'] != 'tongue':
            g['parent'] = key_by_term.get(authoritative[g['term']]['parent_id'], g.get('parent', 'tongue'))
    existing = {c['id'] for c in candidate['concepts']}
    require(not any(g['conceptId'] in existing for g in groups) and not any(c['id'] in existing for c in shared_carriers), 'Do not overwrite native concept identities')
    candidate['concepts'].extend(shared_carriers)
    root_limit = 'Local tongue proposal: native outer boundary and external-muscle sources plus original finite intrinsic territory carriers and styloglossus. Exact term identity does not establish donor anatomy. A shared transverse/vertical carrier, when present, denotes interwoven fibre territory and is counted once. Full muscle segmentation, mucosal coating, body/root partitions, septum, aponeurosis, papillae and remaining terms are not completed.'
    group_by_key = {g['key']:g for g in groups}
    for g in groups:
        concept = dict(id=g['conceptId'], name=g['name'], elements=g['partIds'], ta98Term=g['term'], ta98Kind='exact', ta98Method='explicit-local-tongue-control; anatomy-not-accepted', coverageLimitation=g.get('coverageLimitation', root_limit))
        if g.get('parent'):
            concept['partOf'] = group_by_key[g['parent']]['conceptId']
        candidate['concepts'].append(concept)
    marker = dict(schema='tongue-local/v1', sex='male', wholeConceptId=WHOLE, partIds=part_ids, termGroups=groups,
        bankSha256=digest_file(bank_path), bankUrl=url, carrierBankHashes=features['carrierBankHashes'], sourcePartId='ZA:Tongue',
        sourceCatalogueSha256=source_sha, anatomicalAcceptance=False, completeMuscleVolumeSegmentation=False,
        coverageLimitation=root_limit, originalSourcePartsRetained=True)
    marker.update(representationMode=mode, tissueOwnershipClaim=mode != ENVELOPE_MODE, roiBoundaries='approximate', stage='partial-local-review')
    candidate['finiteTongueIntegration'] = marker
    candidate['triangles'] = sum(p['indexCount']//3 for p in candidate['parts'] if not p.get('suppressed'))
    require(candidate['parts'][:len(original['parts'])] == original['parts'], 'Original source records were modified')
    require(candidate['chunks'][:len(original['chunks'])] == original['chunks'], 'Original source bank references were modified')
    require(next(c for c in candidate['concepts'] if c['id']=='ZA:Tongue') == next(c for c in original['concepts'] if c['id']=='ZA:Tongue'), 'Native whole selection identity changed')
    require(candidate.get('defaultHidden') == original.get('defaultHidden') and all(id in candidate['defaultHidden'] for id in LEGACY_STYLO), 'Legacy visibility changed')
    require(digest_file(catalogue_path) == source_sha, 'Installed source changed during packing')
    present = {g['term']:g for g in groups}
    coverage = copy.deepcopy(inventory)
    coverage.update(schema='ta98-tongue-local-control-inventory/v1', baselineInventorySha256=digest_file(inventory_path), installed=False, anatomicalAcceptance=False, controls=len(groups), physicalMeshes=len(part_ids), newPhysicalMeshes=len(new_parts), geometrySourceCatalogueSha256=source_sha)
    for row in coverage['terms']:
        g = present.get(row['authority']['id'])
        row['candidate_control'] = None if not g else dict(conceptId=g['conceptId'], partIds=g['partIds'], representation=g.get('representation', 'assembly' if g['key'] in ('tongue','muscles') else mode if g['term'] in set(INTRINSIC_TERMS.values()) else 'finite-original-proposal' if g['term']=='A05.1.04.105' else 'native-source-envelope'), anatomicalAcceptance=False, roiBoundaries='approximate', tissueOwnershipClaim=False if mode==ENVELOPE_MODE and g['term'] in set(INTRINSIC_TERMS.values()) else None)
        row['candidate_control_absent'] = g is None
        row['distinct_representation_verified'] = False
    atomic_json(output/'inventory.json', coverage)
    evidence = dict(schema='tongue-local-packed-candidate/v1', source_catalogue=str(catalogue_path), source_catalogue_sha256=source_sha,
        source_audit_sha256=digest_file(source_path), baseline_inventory_sha256=digest_file(inventory_path), intrinsic_construction=str(construction_path), intrinsic_construction_sha256=digest_file(construction_path), intrinsic_checks_sha256=digest_file(checked_path),
        packing_metadata=packing_metadata,
        styloglossus_candidate_sha256=digest_file(stylo_path), styloglossus_checks_sha256=digest_file(context_path), surface_controls_sha256=digest_file(features_path), retained_context_audit=retained_context,
        output_bank_sha256=marker['bankSha256'], output_bank_bytes=len(bank), local_origin_m=ORIGIN.tolist(), records=decoded_checks,
        representation_mode=mode, tissue_ownership_claim=mode != ENVELOPE_MODE, roi_boundaries='approximate', stage='partial-local-review', intrinsic_bone_checks=bone_checks,
        current_source_records_literal=True, source_geometry_and_skeleton_changed=False, legacy_styloglossus_hidden=True, native_whole_selection_retained=True,
        new_image_assets=0, new_uv_assets=0, new_meshes=len(new_parts), assembly_meshes=len(part_ids), represented_controls=len(groups), unrepresented_terms=36-len(groups),
        all_gates_passed=True, visual_selection_record='Required separately by the gated install job; not inferred by this packer', installed=False, anatomical_acceptance=False, complete_muscle_volume_segmentation=False, limitations=[root_limit], elapsed_seconds=time.monotonic()-began)
    remaining()
    atomic_json(output/'packing-evidence.json', evidence)
    atomic_json(output/'atlas.json', candidate)
    print(json.dumps(dict(output=str(output), new_meshes=len(new_parts), bank_bytes=len(bank), controls=len(groups), unrepresented_terms=36-len(groups), all_gates_passed=True, installed=False, elapsed_seconds=evidence['elapsed_seconds'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--intrinsics', type=Path, required=True, help='Explicit reviewed construction.json or its directory; never auto-selected')
    parser.add_argument('--output', type=Path, required=True, help='New directory below the local tongue workbench')
    args = parser.parse_args()
    main(args.intrinsics, args.output)
