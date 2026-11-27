"""Split the existing finite splenic wall at one declared annular boundary.

The distal tissue includes every native daughter and distal main continuation.
No vessel course, source control or organ is fitted. One common cut chart has
opposed ownership; its central lumen is explicitly retained as a hole.
"""
import copy
import json
import signal
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.audit import audit
from geometry_tools.common import digest_file
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank, split_closed_cell
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion'
WORK = BASE / 'splenic-owned-branches-v1'
OUT = WORK / 'candidate-v1'
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
PLAN = BASE / 'native-foregut-source-partition-check-v1/two-cell-hilar-split-plan-v1/split-plan-and-section-proof.json'
PIDS = {'trunk': 'TA98PART:ZA:Splenic artery:trunk',
        'branches': 'TA98PART:ZA:Splenic artery:A12.2.12.051'}
FRAME = 'model:local-ta98:assembled-metres'


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def main():
    signal.alarm(50)
    started = time.monotonic()
    plan = json.loads(PLAN.read_text())
    assert plan['current_catalogue_sha256'] == digest_file(CAT)
    assert json.loads((WORK / 'source-audit.json').read_text())['passed']
    assert (WORK / 'visual-review.json').exists()
    assert not OUT.exists(), 'Retain candidate history'
    OUT.mkdir()
    assembly = dict(schema='source-bound-splenic-two-cell-construction/v1',
                    whole_organ='Existing finite splenic artery, proximal trunk and complete distal hilar source tree',
                    source_catalogue_sha256=digest_file(CAT), split_plan_sha256=digest_file(PLAN),
                    adjoining_interfaces=['One literal common annular tissue cut with opposed trunk/branch ownership; central lumen retained.',
                                         'All eight original arterial mouths retained; coeliac root and spleen perfusion remain under review.'],
                    boundary_status='Provisional trunk-to-terminal-tree boundary before first authored takeoff; no new anatomical branch courses inferred.',
                    all_outside_anatomy_fixed=True, skeleton_basilar_pons_fixed=True,
                    source_courses_and_calibres_fixed=True, anatomical_acceptance=False, installed=False)
    write(OUT / 'construction-plan.json', assembly)
    source = Mesh.load(WORK / 'current-finite-splenic-wall.npz')
    assert source.geometry_sha == plan['current_geometry_sha256']
    v, inverse = np.unique(source.vertices, axis=0, return_inverse=True)
    faces = inverse[source.faces]
    field = (source.attributes['normals'].astype(np.int32) * source.attributes['normalSigns'][:, None]).astype('<i2')[source.faces]
    bank = SharedCutBank(v)
    normal = np.array(plan['selected_plane']['distal_plane_normal'])
    offset = plan['selected_plane']['plane_equation_normal_dot_x_equals']
    split = split_closed_cell(bank, faces, normal, offset, 'splenic-native-station94', allow_holes=True)
    assert len(split.boundary_loops) == 2
    points = np.asarray(bank.positions)
    count = len(split.interface_faces)
    body_faces = {'trunk': split.negative_faces[:-count], 'branches': split.positive_faces[:-count]}
    parents = {'trunk': split.negative_ancestry[:-count], 'branches': split.positive_ancestry[:-count]}
    surfaces = dict(body_faces, interface=split.interface_faces)
    words = {}
    for role in body_faces:
        f, parent = body_faces[role], parents[role]
        original = faces[parent]
        match = f[:, :, None] == original[:, None, :]
        out = np.einsum('tki,tij->tkj', match.astype(np.int32), field[parent].astype(np.int32)).astype('<i2')
        for ti, corner in np.argwhere(~match.any(2)):
            tri = v[faces[parent[ti]]]
            a, b, q = tri[1] - tri[0], tri[2] - tri[0], points[f[ti, corner]] - tri[0]
            aa, ab, bb = a @ a, a @ b, b @ b
            u, w = np.linalg.solve([[aa, ab], [ab, bb]], [q @ a, q @ b])
            n = np.array([1 - u - w, u, w]) @ field[parent[ti]].astype(float)
            out[ti, corner] = np.rint(n / np.linalg.norm(n) * 32767).astype('<i2')
        words[role] = out
    capword = np.rint(normal / np.linalg.norm(normal) * 32767).astype('<i2')
    words['interface'] = np.broadcast_to(capword, (count, 3, 3)).copy()
    cells = {role: [dict(surface=role, side=1),
                    dict(surface='interface', side=1 if role == 'trunk' else -1,
                         neighbor='branches' if role == 'trunk' else 'trunk')]
             for role in PIDS}
    origin = np.array([0., 1.25, 0.])
    packed = pack_finite_cells(points, surfaces, words, cells, PIDS, origin)
    (OUT / 'cells.bin').write_bytes(packed.data)
    atlas = json.loads(CAT.read_text())
    after = copy.deepcopy(atlas)
    old = next(p for p in after['parts'] if p['id'] == 'ZA:Splenic artery')
    limits = ('Source-derived finite arterial region. Trunk/terminal-tree boundary is provisional; '
              'native faceting and angular collars remain. Coeliac attachment, hilar vessel/organ '
              'interfaces and distal perfusion remain unverified.')
    rows = []
    for role in PIDS:
        record = copy.deepcopy(old)
        for key in ['positions', 'normals', 'normalSigns', 'uvs', 'colors', 'indices',
                    'geometryOrigin', 'sourceOffset', 'geometryScale', 'geometryInstance',
                    'tissueDrawGroups', 'suppressed', 'defaultHidden', 'replacementParts']:
            record.pop(key, None)
        record.update(packed.cells[role])
        term = 'A12.2.12.040' if role == 'trunk' else 'A12.2.12.051'
        record.update(conceptId=PIDS[role], name='Splenic artery · main trunk' if role == 'trunk' else 'Splenic branches',
                      chunk=len(after['chunks']), parentId='ZA:Splenic artery',
                      anatomicalParentId='ZA:Splenic artery', ta98Term=term,
                      tissueWallAssembly='splenic-artery:male:owned-branches', tissueWallLayer='wall',
                      coverageLimitation=limits,
                      splenicBranchOwnership=dict(schema='source-bound-hilar-tree/v1', role=role,
                           sourceCatalogueSha256=assembly['source_catalogue_sha256'],
                           splitStation=94, boundaryProvisional=True, allSourceDaughtersRetained=True,
                           reciprocalAnnularInterface=True, organInterfaceVerified=False))
        after['parts'].append(record)
        after['concepts'].append(dict(id=PIDS[role], name=record['name'], elements=[PIDS[role]],
                                     partOf='ZA:Splenic artery', parentId='ZA:Splenic artery',
                                     anatomicalParentId='ZA:Splenic artery', ta98Term=term,
                                     ta98Kind='parent' if role == 'trunk' else 'exact', coverageLimitation=limits))
        rows.append(dict(role=role, part_id=PIDS[role], term=term))
    old.update(suppressed=True, defaultHidden=True, replacementParts=list(PIDS.values()))
    whole = next(c for c in after['concepts'] if c['id'] == 'ZA:Splenic artery')
    whole.update(elements=list(PIDS.values()), componentElements=list(PIDS.values()),
                 ta98Term='A12.2.12.040', ta98Kind='exact', coverageLimitation=limits)
    after['chunks'].append(dict(url='/local-models/' + str((OUT / 'cells.bin').relative_to(ROOT / '.local-models')),
                                bytes=len(packed.data)))
    checks = []
    reader = AtlasReader()
    for row in rows:
        record = next(p for p in after['parts'] if p['id'] == row['part_id'])
        mesh = reader.read_part(after, record, frame=FRAME)
        mesh.save(OUT / (row['role'] + '-packed.npz'))
        check = audit(mesh, require_closed=True, expected_components=1,
                      expected_genus=[1 if row['role'] == 'trunk' else 7],
                      max_pairs=300000, max_tests=1000000)
        write(OUT / (row['role'] + '-audit.json'), check)
        assert check['passed'], check['failures']
        checks.append(check)
    write(OUT / 'proposed-atlas.json', after)
    np.savez_compressed(OUT / 'partition-arrays.npz', positions=points,
                        encoded_positions=packed.positions.astype(float) + origin,
                        source_welded_vertices=v, source_welded_faces=faces,
                        trunk_faces=split.negative_faces, branches_faces=split.positive_faces,
                        interface_faces=split.interface_faces,
                        trunk_parents=split.negative_ancestry, branches_parents=split.positive_ancestry,
                        **{'chart_' + k: value for k, value in packed.surfaces.items()})
    proof = dict(complete=True, passed=True, source_catalogue_sha256=assembly['source_catalogue_sha256'],
                 proposal_sha256=digest_file(OUT / 'proposed-atlas.json'),
                 bank_sha256=digest_file(OUT / 'cells.bin'), bank_bytes=len(packed.data),
                 changed_parts=['ZA:Splenic artery', *PIDS.values()],
                 source_geometry_sha256=source.geometry_sha, rows=rows,
                 split_evidence=split.evidence, source_original_controls_literal=True,
                 all_six_daughters_and_distal_main_continuation_retained=True,
                 native_external_mouths=8, new_shared_interface_has_lumen_hole=True,
                 anatomical_acceptance=False, new_direct_TA98_terms=1,
                 remaining_limits=[limits], elapsed_seconds=time.monotonic() - started)
    write(OUT / 'packing-proof.json', proof)
    print(json.dumps(proof))


if __name__ == '__main__':
    main()
