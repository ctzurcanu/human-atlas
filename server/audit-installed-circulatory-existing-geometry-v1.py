"""Literal baseline/current field audit of reused installed circulatory records."""
from pathlib import Path
import gzip
import hashlib
import json
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, identity
from geometry_tools.model import AtlasReader

OUTPUT = ROOT / '.local-models/geometry-workbench/circulation-census-v1/installed-work-delta-v1'
BASELINE = OUTPUT / 'before-first-circulation-install.json'
CURRENT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
ORDINARY = ROOT / '.local-models/geometry-workbench/circulation-census-v1/ordinary-view-v1/resolved-male.json'
FIRST_MANIFEST = ROOT / '.local-models/ta98-review/male-circulation-renal-fork-v1/manifest.json'
FRAME = 'model:local-ta98:assembled-metres'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_corner_fields(mesh):
    # Canonical oriented triangle positions establish corner order across repacks.
    vertices, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
    faces = inverse[mesh.faces]
    rotations = np.stack([np.roll(faces, shift, axis=1) for shift in range(3)], axis=1)
    selected = np.argmin(rotations[:, :, 0], axis=1)
    canonical = rotations[np.arange(len(faces)), selected]
    order = np.lexsort(canonical.T[::-1])
    corners = (np.arange(3)[None, :] - selected[:, None]) % 3
    fields = []
    for name in ['normals', 'normalSigns', 'uvs', 'colors']:
        values = mesh.attributes.get(name)
        if values is None:
            if name == 'normalSigns':
                values = np.ones(len(mesh.vertices), dtype=np.int8)
            else:
                fields.append((name, None))
                continue
        array = np.asarray(values)[mesh.faces]
        if array.ndim == 2:
            array = array[:, :, None]
        array = np.take_along_axis(array, corners[:, :, None], axis=1)[order]
        fields.append((name, identity(array.astype('<f8'))))
    return identity('canonical-oriented-triangle-corner-fields/v1', *fields)


def main():
    began = time.monotonic()
    assert digest(BASELINE) == '1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39'
    assert digest(CURRENT) == '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
    before, after = json.loads(BASELINE.read_text()), json.loads(CURRENT.read_text())
    old_parts, new_parts = ({part['id']: part for part in atlas['parts']} for atlas in [before, after])
    assert len(old_parts) == len(before['parts']) == 6938
    assert len(new_parts) == len(after['parts']) == 7011
    new_ids = sorted(new_parts.keys() - old_parts.keys())
    removed_ids = sorted(old_parts.keys() - new_parts.keys())
    circ = [part for part in after['parts'] if '/circulatory-completion/' in after['chunks'][part['chunk']]['url']]
    reused = [part for part in circ if part['id'] in old_parts]
    blood = [part for part in reused if part['system'] in ['arterial', 'venous']]
    context = [part for part in reused if part['system'] not in ['arterial', 'venous']]
    assert len(new_ids) == 73 and not removed_ids and len(circ) == 146 and len(reused) == 73 and len(blood) == 71 and len(context) == 2
    assert all(part['id'] in {record['id'] for record in circ} for part in [new_parts[identity] for identity in new_ids])
    ordinary = json.loads(ORDINARY.read_text())
    assert ordinary['sourceSha256'] == digest(CURRENT)
    visibility = {row['id']: row for row in ordinary['effectivePartVisibility']}
    assert all(visibility[part['id']]['available'] for part in circ)
    finite_new = [identity for identity in new_ids if new_parts[identity].get('pairedTissue') and new_parts[identity].get('sectionSolid')]
    assert len(finite_new) == 34
    output = OUTPUT / 'independent-existing-blood-geometry-fields-v1.json'
    assert not output.exists()
    report = {
        'schema': 'installed-circulatory-reused-record-literal-geometry-and-fields/v1',
        'baseline_catalogue': str(BASELINE), 'baseline_sha256': digest(BASELINE),
        'current_catalogue': str(CURRENT), 'current_sha256': digest(CURRENT),
        'ordinary_resolver_sha256': digest(ORDINARY),
        'script_sha256': digest(__file__),
        'baseline_parts': 6938, 'current_parts': 7011, 'new_ids': new_ids, 'removed_ids': removed_ids,
        'current_circulatory_bank_records': 146, 'reused_record_ids': [part['id'] for part in reused],
        'reused_blood_records': 71, 'reused_context_records': context,
        'new_paired_finite_solid_ids': finite_new,
        'new_reference_or_regional_ids': [identity for identity in new_ids if identity not in finite_new],
        'all146_ordinary_available': True,
        'effective_default_hidden_ids': [part['id'] for part in circ if visibility[part['id']]['defaultHidden']],
        'availability_qualification': 'Available means installed selectable source, not initial visibility or anatomical acceptance. Existing18 defaultHidden hepatic owned views are reported separately, never counted as newly added controls.',
        'rows': [], 'complete': False, 'no_geometry_source_model_registry_or_DB_mutation': True,
    }
    atomic_json(output, report)
    reader = AtlasReader(ROOT)
    for part in blood:
        if time.monotonic() - began >= 33.:
            report['stop_reason'] = 'Declared33s internal stop; no automatic extension.'
            break
        old = reader.read_part(before, old_parts[part['id']], frame=FRAME, catalogue_path=BASELINE, catalogue_sha=report['baseline_sha256'])
        new = reader.read_part(after, part, frame=FRAME, catalogue_path=CURRENT, catalogue_sha=report['current_sha256'])
        old_geometry, new_geometry = old.geometry_sha, new.geometry_sha
        old_surface, new_surface = old.surface_sha, new.surface_sha
        old_fields, new_fields = canonical_corner_fields(old), canonical_corner_fields(new)
        geometry_equal = old_surface == new_surface
        fields_equal = old_fields == new_fields
        row = {'id': part['id'], 'name': part['name'], 'system': part['system'],
               'ordinary_visibility': visibility[part['id']],
               'baseline_bank_path': old.source['bank_path'], 'baseline_bank_sha256': old.source['bank_sha256'],
               'current_bank_path': new.source['bank_path'], 'current_bank_sha256': new.source['bank_sha256'],
               'baseline_indexed_geometry_sha256': old_geometry, 'current_indexed_geometry_sha256': new_geometry,
               'baseline_oriented_surface_sha256': old_surface, 'current_oriented_surface_sha256': new_surface,
               'baseline_raw_attributes_sha256': old.attributes_sha, 'current_raw_attributes_sha256': new.attributes_sha,
               'baseline_canonical_corner_fields_sha256': old_fields, 'current_canonical_corner_fields_sha256': new_fields,
               'same_actual_oriented_triangle_surface': geometry_equal,
               'same_canonical_corner_normalwords_signs_UV_colors': fields_equal,
               'indexed_geometry_order_unchanged': old_geometry == new_geometry,
               'baseline_triangles': len(old.faces), 'current_triangles': len(new.faces),
               'classification': 'identical_geometry_and_corner_fields_repacked' if geometry_equal and fields_equal else 'attributes_only_changed' if geometry_equal else 'geometry_changed',
               'change_is_not_anatomical_improvement_claim': True}
        if not geometry_equal:
            # Diagnostic coordinate-set distance, explicitly not a surface Hausdorff certificate.
            old_used, new_used = np.unique(old.vertices, axis=0), np.unique(new.vertices, axis=0)
            distance = max(float(cKDTree(old_used).query(new_used)[0].max()), float(cKDTree(new_used).query(old_used)[0].max()))
            row['bidirectional_nearest_USED_control_distance_max_m'] = distance
            row['coordinate_distance_qualification'] = 'Discrete used-coordinate sets only; not a continuous surface distance, correspondence or anatomical acceptance bound.'
        report['rows'].append(row)
        report['elapsed_seconds'] = time.monotonic() - began
        atomic_json(output, report, overwrite=True)
    else:
        report['complete'] = True
    counts = {label: sum(row['classification'] == label for row in report['rows']) for label in ['geometry_changed', 'attributes_only_changed', 'identical_geometry_and_corner_fields_repacked']}
    report['classification_counts'] = counts
    report['classified_reused_blood_records'] = len(report['rows'])
    first_manifest = json.loads(FIRST_MANIFEST.read_text())
    female_entry = next(row for row in first_manifest if row['path'].endswith('/female/atlas.json'))
    female_backup = ROOT / female_entry['backup']
    baseline_female = gzip.decompress(female_backup.read_bytes())
    current_female = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
    report['female_literal_check'] = {'first_install_backup': str(female_backup),
                                       'baseline_decompressed_sha256': hashlib.sha256(baseline_female).hexdigest(),
                                       'current_sha256': digest(current_female),
                                       'byte_identical': baseline_female == current_female.read_bytes()}
    report['elapsed_seconds'] = time.monotonic() - began
    assert digest(CURRENT) == report['current_sha256']
    atomic_json(output, report, overwrite=True)
    print(json.dumps({'proof': str(output), 'sha256': digest(output), 'complete': report['complete'],
                      'classified': len(report['rows']), 'counts': counts,
                      'female_unchanged': report['female_literal_check']['byte_identical'],
                      'elapsed_seconds': report['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
