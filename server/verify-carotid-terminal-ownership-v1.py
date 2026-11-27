"""Verify ownership-only ICA common-stem/M1 partition from actual bank views.

The source surfaces, normals, UVs, mouths and crossings must be unchanged.
No fitting, geometry construction, database job or install runs here.
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
from geometry_tools.audit import audit
from geometry_tools.model import AtlasReader, Mesh

CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
COHORT = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1/native-carotid-terminal-cohort-v7b'
SOURCE = ['ZA:Middle cerebral artery (M1-segment).l', 'ZA:Middle cerebral artery (M1-segment).r']
H = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def row(p, n, u):
    return min(np.roll(p, k, axis=0).astype('<f8').tobytes() +
               np.roll(n, k, axis=0).astype('<i2').tobytes() +
               np.roll(u, k, axis=0).astype('<f4').tobytes() for k in range(3))


def corners(mesh, face_ids=None):
    faces = mesh.faces if face_ids is None else mesh.faces[face_ids]
    return Counter(row(mesh.vertices[f], mesh.attributes['normals'][f], mesh.attributes['uvs'][f]) for f in faces)


def weld(mesh):
    v, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
    return Mesh(v, inverse[mesh.faces], mesh.frame)


def boundary_edges(mesh):
    mesh = weld(mesh)
    edges = np.vstack([mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]])
    unique, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    return {frozenset(map(tuple, mesh.vertices[e])) for e in unique[counts == 1]}


def run(out):
    started = time.monotonic()
    proposed = out / 'proposed-atlas.json'
    packing = json.loads((out / 'packing-proof.json').read_text())
    assert packing.get('complete') and packing.get('passed')
    before, after = json.loads(CAT.read_text()), json.loads(proposed.read_text())
    source_hash = H(CAT)
    asserted_source = packing.get('source_catalogue_sha256', packing.get('master_sha256'))
    assert asserted_source == source_hash
    old, new = {p['id']: p for p in before['parts']}, {p['id']: p for p in after['parts']}
    bc, ac = {c['id']: c for c in before['concepts']}, {c['id']: c for c in after['concepts']}
    assert len(new) == len(after['parts']) and old.keys() <= new.keys()
    assert bc.keys() == ac.keys(), 'Ownership update must not create duplicate child concepts'
    assert all(new[pid] == p for pid, p in old.items() if pid not in SOURCE)
    assert after['materials'] == before['materials']
    assert after['chunks'][:len(before['chunks'])] == before['chunks']
    new_ids, changed_concepts, side_rows = [], set(), []
    reader, frame = AtlasReader(), 'model:local-ta98:assembled-metres'
    for side, original_id in zip(['l', 'r'], SOURCE):
        ica = f'TA98PART:ZA:Internal carotid artery.{side}:A12.2.06.015'
        stem = f'TA98PART:ZA:Middle cerebral artery (M1-segment).{side}:terminal-ICA'
        m1 = f'TA98PART:ZA:Middle cerebral artery (M1-segment).{side}:A12.2.07.047'
        new_ids.extend([stem, m1]); changed_concepts.update([ica, original_id])
        assert stem in new and m1 in new
        assert new[stem]['conceptId'] == ica and new[m1]['conceptId'] == original_id
        assert not new[stem].get('suppressed') and not new[m1].get('suppressed')
        assert not new[stem].get('defaultHidden') and not new[m1].get('defaultHidden')
        assert ac[ica]['elements'] == bc[ica]['elements'] + [stem]
        assert ac[original_id]['elements'] == [m1]
        assert ac[original_id]['ta98Term'] == 'A12.2.07.047' and ac[original_id]['ta98Kind'] == 'exact'
        assert any(b.get('conceptId') == original_id and b.get('term') == 'A12.2.07.047' and b.get('kind') == 'exact' for b in after['ta98TermBindings'])
        hidden = new[original_id]
        assert hidden.get('suppressed') and hidden.get('defaultHidden') and original_id in after['defaultHidden']
        allowed_metadata = {'suppressed', 'defaultHidden', 'replacementParts', 'replacementPartIds'}
        assert {k: v for k, v in hidden.items() if k not in allowed_metadata} == {k: v for k, v in old[original_id].items() if k not in allowed_metadata}, 'Original source geometry record was altered'
        source = reader.read_part(before, old[original_id], frame=frame)
        frozen = Mesh.load(COHORT / f'mca-{side}.npz', frame=frame)
        assert source.geometry_sha == frozen.geometry_sha and source.attributes_sha == frozen.attributes_sha
        source_stations = np.load(COHORT / f'mca-{side}-stations.npz', allow_pickle=False)
        station = np.full(len(source.vertices), -1, dtype=np.int64)
        station[source_stations['rings']] = np.arange(len(source_stations['rings']))[:, None]
        assert np.all(station >= 0)
        face_station = station[source.faces].min(1)
        first_ids, second_ids = np.flatnonzero(face_station < 17), np.flatnonzero(face_station >= 17)
        prefix, suffix = reader.read_part(after, new[stem], frame=frame), reader.read_part(after, new[m1], frame=frame)
        assert corners(prefix) == corners(source, first_ids)
        assert corners(suffix) == corners(source, second_ids)
        assert corners(prefix) + corners(suffix) == corners(source), 'Assembled surface differs from unchanged native source'
        assert not (corners(prefix) & corners(suffix)), 'Shared faces duplicated between owners'
        seam_controls = source.vertices[source_stations['rings'][17]]
        expected_edges = {frozenset(map(tuple, seam_controls[[i, (i + 1) % len(seam_controls)]])) for i in range(len(seam_controls))}
        shared = boundary_edges(prefix) & boundary_edges(suffix)
        assert len(seam_controls) == len(shared) == 12 and shared == expected_edges
        for mesh in [prefix, suffix]:
            for point, source_index in zip(seam_controls, source_stations['rings'][17]):
                ids = np.flatnonzero(np.all(mesh.vertices == point, axis=1))
                assert len(ids)
                assert all(np.array_equal(mesh.attributes['normals'][i], source.attributes['normals'][source_index]) for i in ids)
                assert all(np.array_equal(mesh.attributes['uvs'][i], source.attributes['uvs'][source_index]) for i in ids)
        vertices = np.concatenate([prefix.vertices, suffix.vertices])
        faces = np.concatenate([prefix.faces, suffix.faces + len(prefix.vertices)])
        merged = weld(Mesh(vertices, faces, frame))
        original_ports = {frozenset(map(tuple, source.vertices[loop])) for loop in boundary_loops(source.faces)}
        actual_ports = {frozenset(map(tuple, merged.vertices[loop])) for loop in boundary_loops(merged.faces)}
        assert len(original_ports) == len(actual_ports) == 2 and original_ports == actual_ports
        check = audit(merged, require_closed=False, expected_components=1, max_pairs=600000, max_tests=2000000)
        assert check['passed'], check['failures']
        side_rows.append(dict(side=side, source_part=original_id, common_stem=stem, M1=m1,
            source_triangles=len(source.faces), common_stem_triangles=len(first_ids), M1_triangles=len(second_ids),
            shared_open_ring_controls=12, two_original_mouths_literal=True,
            all_original_oriented_position_normal_UV_corners_literal=True, decoded_audit=check))
    assert new.keys() - old.keys() == set(new_ids)
    assert len(after['parts']) == len(before['parts']) + 4
    assert all(ac[pid] == concept for pid, concept in bc.items() if pid not in changed_concepts)
    assert set(after['defaultHidden']) - set(before['defaultHidden']) <= set(SOURCE)
    report = dict(schema='native-terminal-source-ownership-independent-replay/v1', complete=True, passed=True,
        source_catalogue_sha256=source_hash, proposed_catalogue_sha256=H(proposed),
        changed_existing_parts=SOURCE, new_parts=new_ids, changed_parts=SOURCE + new_ids,
        unchanged_existing_parts=len(old) - 2, changed_concepts=sorted(changed_concepts),
        new_concepts=0, sides=side_rows, bank_hashes=[dict(url=c['url'], sha256=H(reader.asset(c['url']))) for c in after['chunks'][len(before['chunks']):]],
        all_assembled_native_source_corners_literal=True, no_new_surface_crossings_by_literal_identity=True,
        provisional_boundary_station=17, universal_clinical_boundary_claim=False,
        unchanged_limits=['Original approximately 193 micrometre ICA/MCA seam phase mismatch remains.', 'Original ACA root and ACom contacts/courses are unchanged.', 'No lumen or finite wall repair is claimed.'],
        global_new_TA98_reference_term_count_claim=0, installed=False, anatomical_acceptance=False,
        elapsed_seconds=time.monotonic() - started)
    path = out / 'independent-compiled-replay.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'sides'}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    run(args.directory.resolve())
