"""Finite local wall candidates on literal native foregut arterial courses.

No branch origin, calibre, vessel centreline or organ is fitted. The illustrative
inner envelope contracts each source-derived cross-section radially by 12%.
Directed annular rims close tissue thickness while both lumen mouths stay open.
"""
import argparse
import copy
import json
import signal
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.audit import audit, winding_contains
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.algorithms.tube_continuity import boundary_loops, annular_rim

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/foregut-arterial-batch-v4'
OUT = BASE / 'native-walls-v1'
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME = 'model:local-ta98:assembled-metres'
PARTS = [('ZA:Gastroduodenal artery', 'gda', 'A12.2.12.016'),
         ('ZA:Coeliac trunk', 'coeliac', 'A12.2.12.012')]


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def source_stage():
    OUT.mkdir(exist_ok=False)
    atlas = json.loads(CAT.read_text())
    contract = json.loads((BASE / 'whole-assembly-contract.json').read_text())
    assert digest_file(CAT) == contract['source_catalogue_sha256']
    plan = dict(schema='literal-native-foregut-wall-plan/v1',
                source_catalogue_sha256=digest_file(CAT),
                whole_organ='Coeliac trunk and gastroduodenal artery native courses',
                inventory=str(BASE / 'full-coeliac-authority-inventory.json'),
                source_inventory_sha256=digest_file(BASE / 'full-coeliac-authority-inventory.json'),
                changed_parts=[pid for pid, _, _ in PARTS],
                adjoining_interfaces=['Literal native proximal and distal outer profiles retained for each artery; no new branch ostia or communicating junction claimed.',
                                     'All surrounding organ and vessel geometry fixed. New tissue remains inside the native outer envelope.'],
                radial_wall_fraction=.12, histological_thickness_verified=False,
                skeleton_basilar_pons_fixed=True, all_outside_anatomy_fixed=True,
                anatomical_acceptance=False, installed=False)
    write(OUT / 'plan.json', plan)
    reader = AtlasReader()
    for pid, short, term in PARTS:
        part = next(p for p in atlas['parts'] if p['id'] == pid)
        mesh = reader.read_part(atlas, part, frame=FRAME, catalogue_path=CAT,
                                catalogue_sha=digest_file(CAT))
        courses = recover_station_rings(mesh.vertices, mesh.faces)
        assert len(courses) == 1
        rings = np.array(courses[0]['rings'])
        assert len(rings[0]) == 12
        check = audit(mesh, require_closed=False, expected_components=1)
        write(OUT / (short + '-source-audit.json'), check)
        assert check['passed'], check['failures']
        mesh.save(OUT / (short + '-source.npz'))
        np.savez_compressed(OUT / (short + '-rings.npz'), rings=rings)


def construct():
    began = time.monotonic()
    plan = json.loads((OUT / 'plan.json').read_text())
    assert digest_file(CAT) == plan['source_catalogue_sha256']
    assert (OUT / 'source-render-inspected.json').exists(), 'Inspect actual native source before construction'
    atlas = json.loads(CAT.read_text())
    after = copy.deepcopy(atlas)
    data = bytearray()
    results = []

    def append(a):
        a = np.ascontiguousarray(a)
        while len(data) % a.dtype.itemsize:
            data.append(0)
        offset = len(data)
        data.extend(a.tobytes())
        return offset

    for pid, short, term in PARTS:
        source = Mesh.load(OUT / (short + '-source.npz'))
        rings = np.load(OUT / (short + '-rings.npz'))['rings']
        v, f = source.vertices, source.faces
        centres = v[rings].mean(1)
        inner = v.copy()
        for ring, centre in zip(rings, centres):
            inner[ring] = centre + .88 * (v[ring] - centre)
        nv = len(v)
        wall_faces = [f, f[:, ::-1] + nv]
        loops = boundary_loops(f)
        assert len(loops) == 2
        for loop in loops:
            wall_faces.append(annular_rim(loop, loop + nv))
        faces = np.concatenate(wall_faces)
        # A zero origin preserves every literal native Float32 outer control.
        positions = np.r_[v, inner].astype('<f4')
        assert np.array_equal(positions[:nv].astype(float), v)
        normals = np.r_[source.attributes['normals'],
                        -source.attributes['normals'].astype(np.int32)].astype('<i2')
        uvs = np.r_[source.attributes['uvs'], source.attributes['uvs']].astype('<f4')
        # Retain source outer normals/UVs; rims have their own corner fields.
        cpos = positions[faces].reshape(-1, 3)
        cnorm = normals[faces].reshape(-1, 3)
        cuv = uvs[faces].reshape(-1, 2)
        rim = positions[faces[2 * len(f):]].astype(float)
        rn = np.cross(rim[:, 1] - rim[:, 0], rim[:, 2] - rim[:, 0])
        rn /= np.linalg.norm(rn, axis=1)[:, None]
        cnorm[6 * len(f):] = np.repeat(np.rint(rn * 32767).astype('<i2'), 3, axis=0)
        cuv[6 * len(f):] = 0
        dt = np.dtype([('p', '<f4', (3,)), ('n', '<i2', (3,)), ('uv', '<f4', (2,))])
        keys = np.empty(len(cpos), dtype=dt)
        keys['p'], keys['n'], keys['uv'] = cpos, cnorm, cuv
        _, first, indices = np.unique(keys, return_index=True, return_inverse=True)
        q, n, uv = cpos[first], cnorm[first], cuv[first]
        ff = indices.reshape(-1, 3).astype('<u4')
        layout = dict(positions=append(q), normals=append(n), uvs=append(uv), indices=append(ff))
        part = next(p for p in after['parts'] if p['id'] == pid)
        for key in ['geometryOrigin', 'sourceOffset', 'geometryScale', 'normalSigns', 'colors', 'geometryInstance']:
            part.pop(key, None)
        part.update(layout)
        limitation = ('Illustrative finite arterial wall following the unchanged native course and calibre. '
                      'Both lumen mouths remain open. Histological wall thickness, branch ostia and adjoining '
                      'parent/distal communicating junctions remain unverified; no new branches are implied.')
        part.update(chunk=len(atlas['chunks']), vertexCount=len(q), indexCount=int(ff.size),
                    pairedTissue=True, sectionSolid=True, tissueWallAssembly='foregut:male:' + short,
                    tissueWallLayer='wall', ta98Term=term,
                    tissueDrawGroups=[dict(start=0, count=len(f) * 3),
                                      dict(start=len(f) * 3, count=len(f) * 3),
                                      dict(start=len(f) * 6, count=(len(ff) - len(f) * 2) * 3)],
                    coverageLimitation=limitation,
                    nativeForegutWall=dict(schema='source-section-radial-wall/v1',
                                           sourceCatalogueSha256=plan['source_catalogue_sha256'],
                                           sourceGeometrySha256=source.geometry_sha,
                                           outerControlsAndCornerAppearanceLiteral=True,
                                           radialWallFraction=.12, histologyVerified=False,
                                           communicatingJunctionEstablished=False))
        Mesh(positions.astype(float), faces, FRAME).save(OUT / (short + '-geometric-wall.npz'))
        geometric = Mesh(positions.astype(float), faces, FRAME)
        check = audit(geometric, require_closed=True, expected_components=1, expected_genus=[1])
        write(OUT / (short + '-wall-audit.json'), check)
        assert check['passed'], check['failures']
        probes = []
        for loop in loops:
            centre = positions[nv + loop].astype(float).mean(0)
            hit = winding_contains(geometric, centre[None])
            assert not hit['inside'][0] and not hit['boundary_ambiguous'][0]
            probes.append(centre.tolist())
        assert np.array_equal(q[ff[:len(f)]], v[f].astype('<f4'))
        assert np.array_equal(n[ff[:len(f)]], source.attributes['normals'][f])
        assert np.array_equal(uv[ff[:len(f)]], source.attributes['uvs'][f])
        results.append(dict(part_id=pid, term=term, stations=len(rings),
                            source_outer_faces_literal=True, outer_corner_appearance_literal=True,
                            mouths_open=2, mouth_probe_centres_m=probes,
                            triangles=len(ff), attribute_vertices=len(q),
                            finite_wall_self_crossings=0, anatomical_acceptance=False,
                            remaining_limits=[limitation]))

    bank = OUT / 'native-walls.bin'
    bank.write_bytes(data)
    after['chunks'].append(dict(url='/local-models/' + str(bank.relative_to(ROOT / '.local-models')),
                                bytes=len(data)))
    proposal = OUT / 'proposed-atlas.json'
    write(proposal, after)
    for pid, short, term in PARTS:
        actual = AtlasReader().read_part(after, next(p for p in after['parts'] if p['id'] == pid), frame=FRAME)
        actual.save(OUT / (short + '-packed-wall.npz'))
        check = audit(actual, require_closed=True, expected_components=1, expected_genus=[1])
        write(OUT / (short + '-packed-audit.json'), check)
        assert check['passed'], check['failures']
    proof = dict(complete=True, passed=True, source_catalogue_sha256=plan['source_catalogue_sha256'],
                 proposed_catalogue_sha256=digest_file(proposal), bank_sha256=digest_file(bank),
                 bank_bytes=len(data), changed_parts=plan['changed_parts'], rows=results,
                 new_TA98_terms=0, all_other_parts_literal=True,
                 elapsed_seconds=time.monotonic() - began, installed=False)
    write(OUT / 'packing-proof.json', proof)
    print(json.dumps(proof))


if __name__ == '__main__':
    signal.alarm(40)
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['source', 'construct'])
    args = parser.parse_args()
    source_stage() if args.stage == 'source' else construct()
