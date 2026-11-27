"""Independent replay of actual packed five-part ACA/ACom candidate geometry.

Run on a completed candidate only. This verifier never installs or changes the
current atlas. The source catalogue must still match the construction input.
"""
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
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.audit import audit, topology
from geometry_tools.model import AtlasReader, Mesh

OLD = ['ZA:Anterior cerebral artery.l', 'ZA:Anterior cerebral artery.r',
       'ZA:Anterior communicating artery']
NEW = ['TA98PART:ZA:Anterior cerebral artery.l:A12.2.07.023',
       'TA98PART:ZA:Anterior cerebral artery.l:A12.2.07.033',
       'TA98PART:ZA:Anterior cerebral artery.r:A12.2.07.023',
       'TA98PART:ZA:Anterior cerebral artery.r:A12.2.07.033',
       'TA98PART:ZA:Anterior communicating artery:communicating-course']
TERMS = ['A12.2.07.023', 'A12.2.07.033', 'A12.2.07.023',
         'A12.2.07.033', 'A12.2.07.029']
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
H = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def canonical(p, n=None, uv=None):
    rows = []
    for k in range(3):
        value = np.roll(p, k, axis=0).astype('<f8').tobytes()
        if n is not None:
            value += np.roll(n, k, axis=0).astype('<i2').tobytes()
        if uv is not None:
            value += np.roll(uv, k, axis=0).astype('<f4').tobytes()
        rows.append(value)
    return min(rows)


def mesh_corners(mesh):
    return Counter(canonical(mesh.vertices[f], mesh.attributes['normals'][f],
                             mesh.attributes['uvs'][f]) for f in mesh.faces)


def run(out):
    started = time.monotonic()
    packing = json.loads((out / 'packing-proof.json').read_text())
    assert packing.get('complete') and packing.get('passed')
    assert H(CAT) == packing['source_catalogue_sha256'], 'Current source catalogue changed'
    proposed = out / 'proposed-atlas.json'
    assert H(proposed) == packing['proposed_catalogue_sha256']
    before, after = json.loads(CAT.read_text()), json.loads(proposed.read_text())
    old, new = {p['id']: p for p in before['parts']}, {p['id']: p for p in after['parts']}
    assert len(new) == len(after['parts']) and len(old) == len(before['parts'])
    assert new.keys() - old.keys() == set(NEW) and not old.keys() - new.keys()
    outside = [pid for pid in old if pid not in OLD]
    assert all(new[pid] == old[pid] for pid in outside), 'An unrelated source part record changed'
    for role, pid in enumerate(OLD):
        children = NEW[role * 2:role * 2 + 2] if role < 2 else NEW[4:]
        expected = dict(old[pid], suppressed=True, defaultHidden=True, replacementParts=children)
        assert new[pid] == expected and pid in after['defaultHidden']
    assert after['defaultHidden'][:len(before['defaultHidden'])] == before['defaultHidden']
    assert set(after['defaultHidden']) - set(before['defaultHidden']) <= set(OLD)
    assert after['materials'] == before['materials']
    assert after['chunks'][:-1] == before['chunks']
    assert H(out / 'aca-communicating-exterior.bin') == packing['bank_sha256']
    old_concepts = {c['id']: c for c in before['concepts']}
    new_concepts = {c['id']: c for c in after['concepts']}
    assert len(new_concepts) == len(after['concepts'])
    assert new_concepts.keys() - old_concepts.keys() == set(NEW)
    assert all(new_concepts[pid] == c for pid, c in old_concepts.items() if pid not in OLD)
    assert after['ta98TermBindings'][:-5] == before['ta98TermBindings']
    for i, pid in enumerate(NEW):
        parent = OLD[i // 2] if i < 4 else OLD[2]
        c = new_concepts[pid]
        assert c['elements'] == [pid] and c['partOf'] == c['parentId'] == parent
        assert c['ta98Term'] == TERMS[i] and c['ta98Kind'] == 'exact'
        assert not new[pid].get('suppressed') and not new[pid].get('defaultHidden')
        assert new[pid].get('sourceOuterSurfaceOnly') is True
        assert new[pid]['nativeACACommunicationAssembly']['finiteWall'] is False
    z = np.load(out / 'compiled-shared-geometry.npz', allow_pickle=False)
    v, f, owner = z['vertices'], z['faces'], z['face_owner']
    normal, uv = z['corner_normal_words'], z['corner_uvs']
    assert owner.shape == (len(f),) and set(owner) == set(range(5))
    plan = json.loads((out / 'assembly-plan-before-union.json').read_text())
    assert plan['finite_wall'] is False and plan['physiological_patency_claim'] is False
    assert '2.376 mm gap' in ' '.join(plan['adjoining_interfaces'])
    frame = str(np.load(out / 'composite-input.npz', allow_pickle=False)['frame'])
    reader, decoded = AtlasReader(), []
    for role, pid in enumerate(NEW):
        mesh = reader.read_part(after, new[pid], frame=frame)
        expected = Counter(canonical(v[tri], n, u) for tri, n, u in
                           zip(f[owner == role], normal[owner == role], uv[owner == role]))
        assert mesh_corners(mesh) == expected, 'Actual decoded bank differs from reviewed owned corners'
        decoded.append(mesh)
    # Independently reconstruct one exact-coordinate-welded actual bank mesh.
    offsets = np.cumsum([0] + [len(m.vertices) for m in decoded])
    actual_v = np.concatenate([m.vertices for m in decoded])
    actual_f = np.concatenate([m.faces + offsets[i] for i, m in enumerate(decoded)])
    vertices, inv = np.unique(actual_v, axis=0, return_inverse=True)
    faces = inv[actual_f]
    merged = Mesh(vertices, faces, frame)
    check = audit(merged, require_closed=False, expected_components=1,
                  max_pairs=600000, max_tests=2000000)
    assert check['passed'], check['failures']
    assert topology(merged)['components'][0]['euler'] == -2
    geometric_faces = Counter(canonical(vertices[tri]) for tri in faces)
    assert all(count == 1 for count in geometric_faces.values()), 'A face is owned more than once'
    coordinate_owners = {}
    coordinate_normals = {}
    for role, mesh in enumerate(decoded):
        for point, word in zip(mesh.vertices, mesh.attributes['normals']):
            key = tuple(point)
            coordinate_owners.setdefault(key, set()).add(role)
            coordinate_normals.setdefault(key, set()).add(tuple(word))
    shared = [p for p, owners in coordinate_owners.items() if len(owners) > 1]
    assert shared and all(len(coordinate_normals[p]) == 1 for p in shared), 'Shared controls have inconsistent packed normal words'
    original_mouths, source_port_controls = set(), []
    for role, pid in enumerate(OLD[:2]):
        source = reader.read_part(before, old[pid], frame=frame)
        for loop in boundary_loops(source.faces):
            original_mouths.add(frozenset(map(tuple, source.vertices[loop])))
            for index in loop:
                source_port_controls.append((tuple(source.vertices[index]),
                    tuple(source.attributes['normals'][index]), tuple(source.attributes['uvs'][index])))
    actual_loops = boundary_loops(merged.faces)
    actual_mouths = {frozenset(map(tuple, vertices[loop])) for loop in actual_loops}
    assert len(original_mouths) == len(actual_loops) == 4 and actual_mouths == original_mouths, 'Four original ACA mouths are not literal and open'
    actual_fields = {}
    for mesh in decoded:
        for p, n, u in zip(mesh.vertices, mesh.attributes['normals'], mesh.attributes['uvs']):
            actual_fields.setdefault(tuple(p), set()).add((tuple(n), tuple(u)))
    assert all(actual_fields.get(p) == {(n, u)} for p, n, u in source_port_controls), 'Original mouth normal words or UV fields changed'
    report = dict(schema='independent-ACA-ACom-packed-replay/v1', complete=True, passed=True,
        source_catalogue_sha256=H(CAT), proposed_catalogue_sha256=H(proposed),
        bank_sha256=packing['bank_sha256'], unchanged_existing_parts=len(outside),
        changed_existing_parts=OLD, new_parts=NEW, new_TA98_terms=['A12.2.07.023', 'A12.2.07.033'],
        exactly_five_active_owned_parts=True, every_decoded_face_owned_once=True,
        shared_controls=len(shared), shared_coordinates_and_normal_words_literal=True,
        literal_original_aca_mouths=4, literal_original_mouth_position_normal_UV_controls=len(source_port_controls),
        diagnostic_caps_not_emitted=True, source_outer_surface_only=True,
        finite_wall_or_complete_perfusion=False, proximal_ICA_source_gap_unrepaired_m=.002376,
        all_outside_part_records_literal=True, source_materials_and_unrelated_concepts_literal=True,
        merged_decoded_audit=check, anatomical_acceptance=False, installed=False,
        elapsed_seconds=time.monotonic() - started)
    (out / 'independent-compiled-replay.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: val for k, val in report.items() if k != 'merged_decoded_audit'}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    run(args.directory.resolve())
