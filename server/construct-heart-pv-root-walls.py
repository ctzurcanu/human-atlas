"""Finite source-course pulmonary roots with explicit Boolean provenance.

End disks close computational envelopes only. The installed candidate must
contain an annular tissue rim and an open fluid course, never a fluid plug.
"""
from pathlib import Path
from collections import deque
import argparse, json, sys, time
import numpy as np
import manifold3d as md

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.annular_pipe_envelopes import (
    face_components, boundary_loops, close_with_shared_disks, signed_volume,
)

ap = argparse.ArgumentParser()
ap.add_argument('--stage', required=True)
ap.add_argument('--separate-at-wall', action='store_true',
                help='Classify whole root pieces after wall subtraction; no fluid Boolean')
ap.add_argument('--extended-fluid', action='store_true',
                help='Put computational fluid caps beyond the complete end rim')
ap.add_argument('--outer-envelope', action='store_true',
                help='Clip against the single gated outer heart shell')
ap.add_argument('--root-ownership', action='store_true',
                help='Retain root tissue at host contact for reciprocal host partition')
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
path = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
atlas = json.loads(path.read_text())
reader = AtlasReader()
parts = {p['id']: p for p in atlas['parts']}
frame = 'model:local-ta98:assembled-metres'
origin = np.array([0., 1.375, .03125])
wall_ids = atlas['finiteCardiacIntegration']['wallPartIds']
walls = [reader.read_part(atlas, parts[pid], frame=frame, catalogue_path=path)
         for pid in wall_ids]
v = np.concatenate([p.vertices for p in walls])
offsets = np.cumsum([0] + [len(p.vertices) for p in walls[:-1]])
all_faces = np.concatenate([p.faces + k for p, k in zip(walls, offsets)])
owners = np.concatenate([np.full(len(p.faces), i) for i, p in enumerate(walls)])
ordinals = np.concatenate([np.arange(len(p.faces)) for p in walls])
v, inv = np.unique(v, axis=0, return_inverse=True)
all_faces = inv[all_faces]
_, first, count = np.unique(np.sort(all_faces, axis=1), axis=0,
                            return_index=True, return_counts=True)
assert count.max() == 2
selected = first[count == 1]
f = all_faces[selected]
atomic_npz(out / 'wall-parentage.npz', positions=v, faces=f,
           owner=owners[selected], ordinal=ordinals[selected])
next_id = md.Manifold.reserve_ids(64)
relations = {}


def solid(points, faces, name):
    global next_id
    identity = next_id
    next_id += 1
    faces = np.asarray(faces, dtype=np.uint64)
    result = md.Manifold(md.Mesh64(
        np.asarray(points) - origin, faces,
        run_index=np.array([0, 3 * len(faces)], dtype=np.uint64),
        run_original_id=np.array([identity], dtype=np.uint32),
        face_id=np.arange(len(faces), dtype=np.uint64), tolerance=0.,
    ))
    assert str(result.status()) == 'Error.NoError', (name, result.status())
    relations[identity] = name
    return result, identity


def record(body, name):
    z = body.to_mesh64()
    face_source = np.empty(len(z.tri_verts), int)
    for a, b, identity in zip(z.run_index, z.run_index[1:], z.run_original_id):
        face_source[a // 3:b // 3] = identity
    points = z.vert_properties[:, :3] + origin
    mesh = Mesh(points, z.tri_verts, frame)
    mesh.save(out / (name + '.npz'))
    atomic_npz(out / (name + '-provenance.npz'), original_id=face_source,
               source_face=np.asarray(z.face_id, np.int64))
    return mesh, face_source, np.asarray(z.face_id, np.int64)


tissue, tissue_id = solid(v, f, 'installed-wall')
envelope_mesh = Mesh.load(h / 'blood-spaces-v130/fluid-boundary-1.npz')
envelope, envelope_id = solid(envelope_mesh.vertices,
                             envelope_mesh.faces[:, [0, 2, 1]],
                             'gated-outer-heart-envelope')
fluid_mesh = Mesh.load(h / 'blood-spaces-v130/fluid-boundary-4.npz')
fluid, fluid_id = solid(fluid_mesh.vertices, fluid_mesh.faces,
                        'diagnostic-left-blood-space')
root_ids = ['ZA:Left inferior pulmonary vein', 'ZA:Left superior pulmonary vein',
            'ZA:Right inferior pulmonary vein', 'ZA:Right superior pulmonary vein']
rows = []
root_bodies = []
for root_index, pid in enumerate(root_ids):
    source = reader.read_part(atlas, parts[pid], frame=frame, catalogue_path=path)
    points, inverse = np.unique(source.vertices, axis=0, return_inverse=True)
    faces = inverse[source.faces]
    nc, labels = face_components(faces)
    outer_bodies, inner_bodies = [], []
    for k in range(nc):
        ids, inverse = np.unique(faces[labels == k], return_inverse=True)
        p, fs = points[ids], inverse.reshape(-1, 3)
        loops = boundary_loops(fs)
        assert len(loops) == 2
        adj = [set() for _ in p]
        for t in fs:
            for a, b in zip(t, np.roll(t, -1)):
                adj[a].add(int(b))
                adj[b].add(int(a))
        d = np.full(len(p), -1)
        d[loops[0]] = 0
        queue = deque(map(int, loops[0]))
        while queue:
            a = queue.popleft()
            for b in adj[a]:
                if d[b] < 0:
                    d[b] = d[a] + 1
                    queue.append(b)
        inner = p.copy()
        for j in range(d.max() + 1):
            ring = np.flatnonzero(d == j)
            assert len(ring) == len(loops[0])
            c = p[ring].mean(0)
            inner[ring] = c + .8 * (p[ring] - c)
        for name, coords, target in [('outer', p, outer_bodies),
                                     ('inner', inner, inner_bodies)]:
            vv, ff, disks = close_with_shared_disks(coords, fs)
            if name == 'inner' and args.extended_fluid:
                vv, ff, disks = list(inner), list(fs), []
                for loop in loops:
                    c = inner[loop].mean(0)
                    _, _, vh = np.linalg.svd(p[loop]-c, full_matrices=False)
                    direction = vh[-1]
                    end = int(d[loop[0]])
                    neighbor = inner[d == (1 if end == 0 else d.max()-1)].mean(0)
                    if np.dot(direction, neighbor-c) > 0:
                        direction = -direction
                    # Only computational cutter controls are extended. The
                    # finite source inner-side controls remain unchanged.
                    rim = np.r_[p[loop], inner[loop]]
                    distance = max(0., float(((rim-c)@direction).max()
                        - ((inner[loop]-c)@direction).min())) + .0005
                    new = np.arange(len(vv), len(vv)+len(loop))
                    vv.extend(inner[loop]+distance*direction)
                    tip = len(vv)
                    vv.append(c+distance*direction)
                    for a, b, aa, bb in zip(loop, np.roll(loop, -1),
                                            new, np.roll(new, -1)):
                        ff.extend([[b, a, aa], [b, aa, bb]])
                        disks.append([bb, aa, tip])
                vv, disks = np.asarray(vv), np.asarray(disks, int)
                ff = np.r_[np.asarray(ff, int), disks]
            if signed_volume(vv, ff) < 0:
                ff = ff[:, [0, 2, 1]]
            body, identity = solid(vv, ff, f'{pid}/{k}/{name}')
            atomic_npz(out / f'root-{root_index}-source-{k}-{name}.npz',
                       positions=vv, faces=ff, computational_disk_count=len(disks))
            target.append(body)
    outer = md.Manifold.batch_boolean(outer_bodies, md.OpType.Add)
    inner = md.Manifold.batch_boolean(inner_bodies, md.OpType.Add)
    finite = outer - inner
    outside_wall = finite - (envelope if args.outer_envelope else tissue)
    record(finite, f'root-{root_index}-finite')
    record(outside_wall, f'root-{root_index}-minus-wall')
    if args.root_ownership:
        result = finite - fluid
    elif args.outer_envelope:
        result = outside_wall
    elif args.separate_at_wall:
        import igl
        pieces = outside_wall.decompose()
        piece_records, kept = [], []
        for j, piece in enumerate(pieces):
            pm, _, _ = record(piece, f'root-{root_index}-piece-{j}')
            w = np.abs(igl.fast_winding_number(fluid_mesh.vertices,
                                              fluid_mesh.faces, pm.vertices))
            # Classify complete vertices; a straddling piece requires an
            # explicit local boundary solution and cannot be silently dropped.
            entirely_inside = bool(np.all(w > .5))
            entirely_outside = bool(np.all(w < .5))
            piece_records.append(dict(piece=j, triangles=len(pm.faces),
                volume_mL=piece.volume()*1e6, winding_min=float(w.min()),
                winding_max=float(w.max()), entirely_inside=entirely_inside,
                entirely_outside=entirely_outside))
            if not entirely_inside:
                kept.append(piece)
        atomic_json(out / f'root-{root_index}-pieces.json', piece_records)
        result = md.Manifold.batch_boolean(kept, md.OpType.Add)
    else:
        result = outside_wall - fluid
    assert str(result.status()) == 'Error.NoError', (pid, result.status())
    mesh, sources, parents = record(result, f'root-{root_index}')
    root_bodies.append(result)
    proof = audit(mesh, require_closed=True, max_pairs=1000000,
                  max_tests=7000000)
    atomic_json(out / f'root-{root_index}-audit.json', proof)
    contact = sources == tissue_id
    fluid_cut = sources == fluid_id
    atomic_npz(out / f'root-{root_index}-wall-interface.npz',
               faces=mesh.faces[contact], parent_faces=parents[contact])
    row = dict(part_id=pid, source_geometry_sha256=source.geometry_sha,
               source_components=nc, triangles=len(mesh.faces),
               geometry_sha256=mesh.geometry_sha, finite_volume_mL=finite.volume()*1e6,
               trimmed_volume_mL=result.volume()*1e6,
               wall_occupied_overlap_mL=(result ^ tissue).volume()*1e6,
               left_blood_occupied_overlap_mL=(result ^ fluid).volume()*1e6,
               wall_contact_triangles=int(contact.sum()),
               fluid_trim_triangles=int(fluid_cut.sum()), passed=proof['passed'],
               failures=proof['failures'])
    rows.append(row)
    atomic_json(out / 'progress.json', dict(roots=rows,
               elapsed_seconds=time.monotonic()-start), overwrite=True)
    print(json.dumps(row), flush=True)
pairs = []
for i in range(4):
    for j in range(i+1, 4):
        overlap = root_bodies[i] ^ root_bodies[j]
        pairs.append(dict(first=root_ids[i], second=root_ids[j],
                          occupied_overlap_mL=overlap.volume()*1e6))
atomic_json(out / 'construction.json', dict(
    passed=all(r['passed'] for r in rows), roots=rows, root_pairs=pairs,
    source_catalogue_sha256=digest_file(path), wall_ids=wall_ids,
    wall_original_id=tissue_id, fluid_original_id=fluid_id,
    original_names=relations, source_wall_fixed=True, inner_radial_fraction=.8,
    anatomical_thickness_verified=False, reciprocal_wall_retriangulation_pending=True,
    installed=False, elapsed_seconds=time.monotonic()-start))
