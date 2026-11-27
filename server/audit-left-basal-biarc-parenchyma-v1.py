"""Read-only sampled placement in the actual five-cell left inferior lobe union."""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from geometry_tools.winding_backend import winding_numbers

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
CONSTRUCTION = BASE / 'biarc-field-finite-wall-v1'
WALL = CONSTRUCTION / 'candidate-stage-v2/finite-basal-wall.npz'
AXES = CONSTRUCTION / 'independent-final-v2/continuous-axis-patency-samples.npz'
FIELDS = BASE / 'whole-unbranched-root-stem-biarc-v1/proper-profile-transport-v1/course-fields.npz'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
OUTPUT = WALL.parent / 'independent-fixed-context-v1'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    began = time.monotonic()
    expected = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
    assert digest(CATALOGUE) == expected
    assert digest(WALL) == 'e5e65bc043e92de2a179a5ef65f8db4ea6b5ef5c7b0267cb019862cb2f8e810e'
    reg = json.loads((BASE / 'left-lower-source-recipient-registration-v2/registration-contract.json').read_text())
    atlas = json.loads(CATALOGUE.read_text())
    parts = {part['id']: part for part in atlas['parts']}
    identities = [identity for identity in reg['bindings'] if identity.startswith('TA98PULMONARY:male:left-inferior:')]
    assert len(identities) == 5
    reader = AtlasReader(ROOT)
    vertices, faces, offset, bindings = [], [], 0, []
    for identity in identities:
        part = parts[identity]
        mesh = reader.read_part(atlas, part, frame='model:local-ta98:assembled-metres', catalogue_path=CATALOGUE, catalogue_sha=expected)
        assert mesh.geometry_sha == reg['bindings'][identity]['geometry_sha256']
        keep = np.concatenate([np.arange(group['start'] // 3, (group['start'] + group['count']) // 3)
                               for group in part['tissueDrawGroups'] if group.get('neighbor') not in identities])
        vertices.append(mesh.vertices)
        faces.append(mesh.faces[keep] + offset)
        offset += len(mesh.vertices)
        bindings.append({'id': identity, 'geometry_sha256': mesh.geometry_sha, 'bank_sha256': mesh.source['bank_sha256']})
    vertices, inverse = np.unique(np.concatenate(vertices), axis=0, return_inverse=True)
    target = Mesh(vertices, inverse[np.concatenate(faces)], 'model:local-ta98:assembled-metres')
    target_topology = topology(target)
    assert target_topology['closed_oriented_manifold']
    wall = Mesh.load(WALL)
    with np.load(FIELDS, allow_pickle=False) as archive:
        fields = {key: archive[key].copy() for key in archive.files}
    with np.load(AXES, allow_pickle=False) as archive:
        axis_points = archive['midpoints_m'].copy()
        branch_points = archive['branch_anchor_points_m'].copy()
    slices = {int(curve): (int(first), int(last)) for curve, first, last in fields['curve_slices']}
    graph = json.loads((BASE / 'collector-left-common-basal-v1/contract-before-construction.json').read_text())['component_graph']
    tips = []
    for source in graph:
        begin, end = slices[source['curve']]
        tips.append(end - 1 if source['native_root_endpoint'] == 0 else begin)
    groups = {'wall_vertices': wall.vertices, 'all404_centres': fields['all_centres'],
              'continuous_patency_segment_midpoints': axis_points,
              'seven_branch_anchors': branch_points,
              'eight_distal_tip_centres': fields['all_centres'][tips],
              '96_outer_distal_controls': fields['all_outer'][tips].reshape(-1, 3),
              '96_inner_distal_controls': fields['all_inner'][tips].reshape(-1, 3)}
    numeric, arrays = {}, {}
    for label, points in groups.items():
        winding = np.abs(winding_numbers(target.vertices, target.faces, points))
        inside = winding > .5
        numeric[label] = {'total': len(points), 'inside': int(inside.sum()), 'outside': int((~inside).sum()),
                          'outside_indices': np.flatnonzero(~inside).tolist()}
        arrays[label + '_abs_winding'] = winding
        arrays[label + '_outside_points'] = points[~inside]
    modified = fields['biarc_modified_interior']
    center_inside = arrays['all404_centres_abs_winding'] > .5
    numeric['modified_biarc_interior_centres'] = {'total': int(modified.sum()), 'inside': int((center_inside & modified).sum()), 'outside': int((~center_inside & modified).sum())}
    witnesses = OUTPUT / 'parenchyma-winding-witnesses.npz'
    assert not witnesses.exists()
    np.savez_compressed(witnesses, **arrays)
    proof = {'schema': 'canonical-biarc-wall-current-left-inferior-union-sampled-placement/v1',
             'source_catalogue_sha256': expected, 'candidate_wall_sha256': digest(WALL),
             'course_fields_sha256': digest(FIELDS), 'patency_axis_samples_sha256': digest(AXES),
             'target_geometry_sha256': target.geometry_sha, 'target_bindings': bindings,
             'target_topology': target_topology, 'counts': numeric,
             'method': 'Remove only exact internal tissueDrawGroup interfaces among the five fixed HRA territory cells, exact-coordinate weld, require closed oriented union; abs fast-winding>0.5 at saved actual candidate points.',
             'scope': 'Numerical point containment only; no continuous whole-wall containment or pleural occupied-volume claim. Hilar root tissue may legitimately be outside parenchyma and needs anatomical review; no point movement follows this report.',
             'internal_territory_crossings_are_not_outer_pleura_crossings': True,
             'witnesses_sha256': digest(witnesses), 'no_source_or_geometry_edits': True,
             'anatomical_acceptance': False, 'elapsed_seconds': time.monotonic() - began}
    output = OUTPUT / 'parenchyma-placement.json'
    output.write_text(json.dumps(proof, indent=2) + '\n')
    assert digest(CATALOGUE) == expected
    print(json.dumps({'proof': str(output), 'sha256': digest(output),
                      'counts': {label: {key: value for key, value in row.items() if key != 'outside_indices'} for label, row in numeric.items()},
                      'elapsed_seconds': proof['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
