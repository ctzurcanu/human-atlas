"""Bounded storage rounding of new interior serosal controls only."""
from pathlib import Path
import argparse, itertools, json, sys, time
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.tissue import normal_words

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True)
ap.add_argument('--stage', required=True)
ap.add_argument('--joint-words', action='store_true')
ap.add_argument('--shared-new-boundary-words', action='store_true')
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
src, out = h / args.input, h / args.stage
out.mkdir(exist_ok=False)
meta = json.loads((src / 'ownership.json').read_text())
z = np.load(src / 'arrangement.npz')
arrays = {k: z[k].copy() for k in z.files}
p = arrays['positions']
q = p.astype('<f4').astype(float)
initial = q.copy()
row = next(r for r in meta['surfaces'] if r['id'].startswith('owned-interface:'))
f = arrays[row['faces']]
protected = np.unique(np.concatenate([arrays[r['faces']].ravel() for r in meta['surfaces'] if r != row]))
source = np.load(h / 'annular-crossing-round-v376/arrangement.npz')
prefix = len(source['positions'])
fixed = np.zeros(len(p), bool)
fixed[:prefix] = True
fixed[protected] = True
tri = p[f]
reference = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
def good(ids):
    t = q[f[ids]]
    n = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
    return np.einsum('fi,fi->f', n, reference[ids]) > 0
directions = sorted([v for v in itertools.product([-1, 0, 1], repeat=3) if any(v)],
                    key=lambda v: sum(x*x for x in v))
tree = cKDTree(initial)
moves = []
if args.joint_words:
    from geometry_tools.algorithms.quantized_orientation_repair import repair_orientation_words
    if args.shared_new_boundary_words:
        fixed[:] = False
        fixed[:prefix] = True
        all_faces, all_words = [], []
        for surface in meta['surfaces']:
            sf = arrays[surface['faces']]
            all_faces.append(sf)
            if surface['id'].startswith('owned-interface:'):
                tr = p[sf]
                fn = np.cross(tr[:, 1]-tr[:, 0], tr[:, 2]-tr[:, 0])
                all_words.append(np.repeat(fn[:, None], 3, axis=1))
            else:
                all_words.append(arrays[surface['normals']].astype(float))
        q, joint_proof = repair_orientation_words(p, np.concatenate(all_faces), fixed,
                    maximum_m=2e-8, seconds=25., corner_normals=np.concatenate(all_words))
    else:
        q, joint_proof = repair_orientation_words(p, f, fixed, maximum_m=2e-8, seconds=25.)
    moves = joint_proof['iterations']
else:
    joint_proof = None
for iteration in range(0 if args.joint_words else 16):
    bad = np.flatnonzero(~good(np.arange(len(f))))
    if not len(bad):
        break
    if time.monotonic()-start > 20:
        raise TimeoutError('Serosal storage-rounding search deadline')
    candidates = np.unique(f[bad])
    candidates = candidates[~fixed[candidates]]
    accepted = False
    for vertex in candidates:
        incident = np.flatnonzero(np.any(f == vertex, axis=1))
        before = int(np.sum(~good(incident)))
        old = q[vertex].copy()
        near = tree.query_ball_point(old, 4e-8)
        for steps in [1, 2, 4, 8]:
            for direction in directions:
                value = old.astype('<f4')
                for axis, sign in enumerate(direction):
                    for step in range(steps if sign else 0):
                        value[axis] = np.nextafter(value[axis], np.float32(sign*np.inf))
                value = value.astype(float)
                if np.linalg.norm(value-p[vertex]) > 2e-8:
                    continue
                if any(i != vertex and np.array_equal(q[i], value) for i in near):
                    continue
                q[vertex] = value
                after = int(np.sum(~good(incident)))
                if after < before and np.all(good(incident)[good(incident) | np.isin(incident, bad)]):
                    # No previously valid incident face may reverse.
                    q[vertex] = old
                    valid_before = good(incident)
                    q[vertex] = value
                    if not good(incident)[valid_before].all():
                        continue
                    moves.append(dict(vertex=int(vertex), from_m=old.tolist(), to_m=value.tolist(),
                                      previous_invalid_faces=before, remaining_invalid_faces=after))
                    accepted = True
                    break
                q[vertex] = old
            if accepted:
                break
        if accepted:
            break
        q[vertex] = old
    if not accepted:
        raise ValueError('No bounded internal Float32 correction for the remaining serosal facets')
else:
    if not args.joint_words:
        raise ValueError('Serosal storage correction exhausted its finite move budget')
assert np.array_equal(q[fixed], initial[fixed])
arrays['positions'] = q
arrays[row['normals']] = normal_words(q, {row['id']: f})[row['id']]
atomic_npz(out / 'arrangement.npz', **arrays)
meta['source_stage'] = args.stage
atomic_json(out / 'ownership.json', meta)
atomic_json(out / 'rounding.json', dict(
    passed=True, moves=moves, source_stage=args.input,
    moved_control_ids=np.flatnonzero(np.any(q != initial, axis=1)).tolist(),
    joint_word_proof=joint_proof,
    source_arrangement_sha256=digest_file(src / 'arrangement.npz'),
    maximum_total_rounding_m=float(np.linalg.norm(q-p, axis=1).max()),
    maximum_internal_word_correction_m=float(np.linalg.norm(q-initial, axis=1).max()),
    all_original_source_controls_literal=np.array_equal(q[:prefix], p[:prefix]),
    all_boundary_controls_fixed_after_encoding=not args.shared_new_boundary_words,
    new_shared_cut_boundary_words_may_move=args.shared_new_boundary_words,
    face_indices_literal=True, full_self_topology_and_occupied_volume_audits_pending=True,
    elapsed_seconds=time.monotonic()-start, anatomical_acceptance=False, installed=False))
print('Serosal storage repaired:', int(np.sum(np.any(q != initial, axis=1))), 'new controls;',
      np.linalg.norm(q-initial, axis=1).max()*1e9, 'nm maximum word correction;',
      time.monotonic()-start, 's; full audits pending', flush=True)
