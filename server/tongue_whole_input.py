"""Hash-bound optional corrected whole-tongue inputs for local regeneration.

This module validates saved evidence; it never constructs, fits or installs a
mesh. A corrected input must preserve the native source face/control mapping so
the existing anatomical dorsum selection remains an explicit source selection.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from geometry_tools.common import digest_file
from geometry_tools.model import Mesh

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / '.local-models/geometry-workbench/tongue-completion'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
SOURCE_AUDIT = WORK / 'source-v1/source-audit.json'
NATIVE_WHOLE = WORK / 'source-v1/meshes/native-00.npz'


def _require(value, message):
    if not value:
        raise ValueError(message)


def _local_file(value):
    path = Path(value).resolve()
    _require(path.is_relative_to(WORK) and path.is_file(),
             f'Expected an existing file inside the local tongue workbench: {path}')
    return path


def fixed_bone_sources():
    """Return the actual eight immutable bone exports from the source audit."""
    audit = json.loads(SOURCE_AUDIT.read_text())
    rows = [r for r in audit['contextParts'] if r['sourcePart'].get('system') == 'skeletal']
    _require(len(rows) == 8, 'Expected the eight declared native fixed bones')
    for row in rows:
        path = _local_file(row['meshSavePath'])
        _require(digest_file(path) == row['meshSaveSha256'], 'Fixed bone export hash changed')
        _require(Mesh.load(path).geometry_sha == row['geometrySha256'], 'Fixed bone geometry changed')
        bank = Path(row['bankPath'])
        _require(digest_file(bank) == row['bankSha256'], 'Fixed bone source bank changed')
    return rows


def whole_input(whole=None, fit_result=None):
    """Load native defaults, or an exact successful source-bound fitted whole.

    Passing the two explicit paths selects a previously reviewed correction.
    It does not create a new anatomical-acceptance assertion. Derived cells must
    still pass their own topology, containment, contact and visual checks.
    """
    _require((whole is None) == (fit_result is None),
             '--whole-tongue and --fit-result must be supplied together')
    native = Mesh.load(NATIVE_WHOLE)
    audit = json.loads(SOURCE_AUDIT.read_text())
    source_row = next(r for r in audit['contextParts'] if r['partId'] == 'ZA:Tongue')
    catalogue_sha = digest_file(CATALOGUE)
    _require(catalogue_sha == native.source['catalogue_sha256'] ==
             audit['sourceCatalogues']['male']['sha256'], 'Parent catalogue binding changed')
    _require(digest_file(NATIVE_WHOLE) == source_row['meshSaveSha256'] and
             native.geometry_sha == source_row['geometrySha256'], 'Native whole export changed')
    _require(digest_file(Path(source_row['bankPath'])) == source_row['bankSha256'] ==
             native.source['bank_sha256'], 'Parent tongue bank changed')
    binding = dict(schema='tongue-reviewed-whole-input/v1', mode='native-source-default',
        source_catalogue_path=str(CATALOGUE), source_catalogue_sha256=catalogue_sha,
        source_audit_path=str(SOURCE_AUDIT), source_audit_sha256=digest_file(SOURCE_AUDIT),
        parent_whole_mesh=str(NATIVE_WHOLE), parent_whole_mesh_sha256=digest_file(NATIVE_WHOLE),
        parent_whole_geometry_sha256=native.geometry_sha, parent_source_bank_sha256=source_row['bankSha256'],
        whole_mesh=str(NATIVE_WHOLE), whole_mesh_sha256=digest_file(NATIVE_WHOLE),
        whole_geometry_sha256=native.geometry_sha, frame=native.frame,
        original_source_files_modified=False, fixed_bones_modified=False,
        anatomical_acceptance=False, derived_geometry_requires_independent_checks=True)
    if whole is None:
        return native, NATIVE_WHOLE, binding

    path, evidence = _local_file(whole), _local_file(fit_result)
    report = json.loads(evidence.read_text())
    _require(report.get('schema') == 'tongue-protected-basal-cage/v1', 'Unrecognized whole fit evidence schema')
    _require(report.get('passed') is True and report.get('audit', {}).get('passed') is True,
             'Corrected whole fit is not independently complete and passed')
    _require(Path(report['mesh']).resolve() == path and report['mesh_sha256'] == digest_file(path),
             'Fit result does not bind the exact selected whole NPZ')
    _require(Path(report['source']).resolve() == NATIVE_WHOLE and
             report['source_sha256'] == binding['parent_whole_mesh_sha256'] and
             report['source_geometry_sha256'] == native.geometry_sha and
             report['source_catalogue_sha256'] == catalogue_sha and
             report['source_bank_sha256'] == binding['parent_source_bank_sha256'],
             'Fit parent source/catalogue/bank binding differs from the audited source')
    _require(report.get('bone_geometry_modified') is False and report.get('source_file_modified') is False,
             'Fit must preserve literal fixed bones and original source files')
    _require(report.get('global_jacobian_lower_bound', 0) > 0,
             'Fit lacks a positive whole-field Jacobian bound')
    corrected = Mesh.load(path)
    _require(corrected.geometry_sha == report['geometry_sha256'] and corrected.frame == native.frame,
             'Selected corrected whole geometry/frame does not match its fit result')
    _require(corrected.vertices.shape == native.vertices.shape and np.array_equal(corrected.faces, native.faces),
             'Corrected whole must retain source face/control correspondence; regenerate anatomical masks explicitly otherwise')
    _require(corrected.source.get('catalogue_sha256') == catalogue_sha and
             corrected.source.get('bank_sha256') == binding['parent_source_bank_sha256'] and
             corrected.source.get('derived_from_geometry') == native.geometry_sha,
             'Corrected mesh provenance does not retain its parent authority bindings')
    normals = corrected.attributes.get('normals')
    _require(normals is not None and normals.shape == corrected.vertices.shape and
             np.isfinite(normals).all() and np.all(np.linalg.norm(normals, axis=1) > 0),
             'Corrected whole needs its actual finite nonzero transported normal field')
    bones = fixed_bone_sources()
    expected = {r['partId']: r for r in bones}
    contacts = report.get('bone_contacts', [])
    _require(len(contacts) == len(expected) and len({r['part_id'] for r in contacts}) == len(expected),
             'Fit must contain each fixed-bone pair exactly once')
    _require(set(report.get('fixed', [])) == set(expected), 'Fit fixed-bone identities differ from the source audit')
    for row in contacts:
        source = expected.get(row['part_id'])
        summary = row['summary']
        _require(source is not None and row['geometry_sha256'] == source['geometrySha256'] and
                 summary.get('complete_within_pair_budget') is True and
                 summary.get('crossing_face_pairs') == 0 and summary.get('coplanar_overlap_face_pairs') == 0,
                 f'Corrected whole has a failed, stale or incomplete fixed-bone pair: {row["part_id"]}')
    review_path = _local_file(evidence.parent / 'render-review.json')
    review = json.loads(review_path.read_text())
    _require(review.get('schema') == 'tongue-whole-prepartition-review/v1' and
             review.get('accepted_for_derived_construction') is True and
             review.get('fit_result_sha256') == digest_file(evidence) and
             review.get('mesh_sha256') == digest_file(path),
             'Corrected whole needs a source-bound actual surface/section review accepted for derived construction')
    views = review.get('views', [])
    _require(any(v.get('kind') == 'surface' for v in views) and
             any(v.get('kind') == 'section' for v in views), 'Whole review must include actual surfaces and sections')
    for view in views:
        image = Path(view['path']).resolve()
        _require(image.is_file() and (image.is_relative_to(WORK) or image.is_relative_to(WORK.parent / 'review')),
                 'Whole review image must be an existing local workbench render')
        _require(view.get('inspected') is True and digest_file(image) == view['sha256'],
                 'Whole review contains an uninspected or changed render')
    binding.update(mode='explicit-reviewed-corrected-whole', whole_mesh=str(path),
        whole_mesh_sha256=digest_file(path), whole_geometry_sha256=corrected.geometry_sha,
        fit_result_path=str(evidence), fit_result_sha256=digest_file(evidence),
        fit_schema=report['schema'], fit_numerical_and_fixed_bone_checks_passed=True,
        render_review_path=str(review_path), render_review_sha256=digest_file(review_path),
        source_face_correspondence_preserved=True,
        review_selection='Explicit caller-selected corrected whole; anatomical acceptance is not inferred from numerical evidence')
    return corrected, path, binding


def checked_inner_envelope(evidence_path, native, binding):
    """Reject a stale erosion of the pre-correction source."""
    path = _local_file(evidence_path)
    evidence = json.loads(path.read_text())
    _require(evidence.get('passed') is True and
             evidence.get('source_catalogue_sha256') == binding['source_catalogue_sha256'] and
             evidence.get('source_mesh_sha256') == binding['whole_mesh_sha256'] and
             evidence.get('source_geometry_sha256') == native.geometry_sha,
             'Inner erosion must be passed and derived from the exact selected whole')
    mesh_path = _local_file(evidence['mesh'])
    _require(digest_file(mesh_path) == evidence['mesh_sha256'], 'Inner erosion NPZ hash changed')
    inner = Mesh.load(mesh_path)
    _require(inner.geometry_sha == evidence['geometry_sha256'] and inner.frame == native.frame,
             'Inner erosion geometry/frame changed')
    if binding['mode'] != 'native-source-default':
        other = evidence.get('whole_input', {})
        _require(other.get('fit_result_sha256') == binding['fit_result_sha256'] and
                 other.get('whole_geometry_sha256') == native.geometry_sha,
                 'Inner erosion does not retain the corrected-whole fit evidence')
    return inner, path, evidence


def new_output(path):
    path = Path(path).resolve()
    _require(path.is_relative_to(WORK) and path != WORK and not path.exists(),
             'Use a new local tongue stage; existing geometry/evidence is immutable')
    return path
