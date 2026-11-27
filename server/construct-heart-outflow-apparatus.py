"""Keep each arterial root and its three cusps together during a bounded fit."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import differential_evolution
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.fit import transform_mesh, constrained_translation, landmark_transform
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.annular_pipe_envelopes import close_with_shared_disks
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.conforming_assembly_ownership import subtract_assembly
from geometry_tools.algorithms.conforming_fields import replay_corner_fields
from geometry_tools.algorithms.surface_role_replay import split_by_source_role

ap = argparse.ArgumentParser()
ap.add_argument('--aorta', required=True)
ap.add_argument('--pulmonary', required=True)
ap.add_argument('--stage', required=True)
ap.add_argument('--native-intersections',action='store_true')
ap.add_argument('--placement', help='Replay a prior bounded root-group translation, then recheck all envelopes')
ap.add_argument('--actual-envelope-search', action='store_true',
                help='Fit against literal closed envelopes instead of conservative convex proxies')
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
path = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
atlas = json.loads(path.read_text())
parts = {p['id']: p for p in atlas['parts']}
reader = AtlasReader()
origin = np.array([0., 1.375, .03125])
roots, envelopes, groups, inlet_frames = [], [], [], []
root_face_counts={}
for stage, valve in [(args.aorta, 'aortic-valve'), (args.pulmonary, 'pulmonary-valve')]:
    meta = json.loads((h / stage / 'construction.json').read_text())
    assert meta['passed'] and meta['source_catalogue_sha256'] == digest_file(path)
    root = Mesh.load(h / stage / 'finite-root.npz')
    roots.append((meta['source_part_id'], root))
    z = np.load(h / stage / 'source-walls.npz')
    root_face_counts[meta['source_part_id']]=len(z['source_faces'])
    native=reader.read_part(atlas,parts[meta['source_part_id']],frame=root.frame,catalogue_path=path)
    original,lookup=np.unique(native.vertices,axis=0,return_inverse=True)
    if meta.get('root_continuation'):
        inlet=np.eye(4);registration=dict(native_cusp_hinges_preserved=True,maximum_landmark_error_m=0.)
    else:
        assert np.array_equal(lookup[native.faces],z['source_faces'])
        proximal=min(boundary_loops(z['source_faces']),key=lambda ring:z['outer_positions'][ring].mean(0)[1])
        inlet,registration=landmark_transform(original[proximal],z['outer_positions'][proximal])
        assert registration['maximum_landmark_error_m']<1e-12
    inlet_frames.append((inlet,registration))
    v, f, _ = close_with_shared_disks(z['outer_positions'], z['source_faces'])
    envelopes.append(Mesh(v, f, root.frame))
    cusps = next(g['partIds'] for g in atlas['finiteCardiacIntegration']['termGroups'] if g['key'] == valve)
    assert len(cusps) == 3
    groups.append([meta['source_part_id']] + cusps)
# Four bounded transverse translation variables; each root's three cusps
# receive its same transform. The SVC and all other anatomy remain fixed.
protected = np.load(h / 'svc-rounded-rims-v267/source-walls.npz')
pv, pf, _ = close_with_shared_disks(protected['outer_positions'], protected['source_faces'])
svc_envelope = Mesh(pv, pf, envelopes[0].frame)
def conservative_clearance_envelope(mesh):
    # Construction-only convex hull plus a 26-direction dilation. This
    # deliberately overbounds each actual envelope; the final finite-triangle
    # checks still use the uninflated source surfaces.
    directions = np.array([q for q in np.ndindex(3, 3, 3) if q != (1, 1, 1)], float) - 1
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    points = (mesh.vertices[:, None] - origin + .0004 * directions[None]).reshape(-1, 3)
    hull = ConvexHull(points)
    ff = hull.simplices.copy()
    tri = points[ff]
    normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    reverse = np.einsum('ti,ti->t', normal, hull.equations[:, :3]) < 0
    ff[reverse] = ff[reverse][:, [0, 2, 1]]
    return make_solid(points, ff)[0]
proxies = ([make_solid(m.vertices - origin, m.faces)[0] for m in envelopes + [svc_envelope]]
           if args.actual_envelope_search else
           [conservative_clearance_envelope(m) for m in envelopes + [svc_envelope]])
trials, best = [], [float('inf'), None]
def merit(values):
    if time.monotonic() - start > 18:
        atomic_json(out / 'search.json', dict(trials=trials, best_variables_mm=best[1], passed=False), overwrite=True)
        raise TimeoutError('Four-variable protected arterial search deadline')
    translations = [np.array([values[i], 0., values[i + 1]]) / 1000 for i in [0, 2]]
    sizes = np.linalg.norm(translations, axis=1)
    if np.any(sizes > .008): return 1e12 + float(np.sum(sizes ** 2))
    a, p = [solid.translate(t) for solid, t in zip(proxies, translations)]
    amounts = [max(0., float((a ^ p).volume())), max(0., float((a ^ proxies[2]).volume())), max(0., float((p ^ proxies[2]).volume()))]
    score = (1e13 * sum(amounts) + 1e5 * sum(v > 0 for v in amounts)
             + float(np.sum(values ** 2)))
    trials.append(dict(translations_m=[t.tolist() for t in translations], clearance_proxy_overlaps_m3=amounts, score=score))
    if score < best[0]: best[:] = [score, np.asarray(values).tolist()]
    return score
if args.placement:
    prior = json.loads((h / args.placement / 'placement.json').read_text())
    assert prior['passed'] and prior['dependency_groups'] == groups
    matrices = [np.asarray(m) for m in prior['transforms']]
    assert all(np.array_equal(m[:3, :3], np.eye(3)) and m[1, 3] == 0 for m in matrices)
    chosen = np.array([value for m in matrices for value in m[[0, 2], 3] * 1000])
else:
    rng = np.random.default_rng(7286)
    population = rng.uniform(-5., 5., (24, 4))
    population[0] = [2., 5., 7., 3.]
    result = differential_evolution(merit, [(-8., 8.)] * 4, init=population, maxiter=14, seed=7286, polish=False, tol=1e-5)
    chosen = np.array(best[1])
translations = [np.array([chosen[i], 0., chosen[i + 1]]) / 1000 for i in [0, 2]]
distance = max(np.linalg.norm(t) for t in translations)
atomic_json(out / 'search.json', dict(trials=trials, best_variables_mm=chosen.tolist(), elapsed_seconds=time.monotonic() - start))
actual_solids = [make_solid(m.vertices - origin + t, m.faces)[0] for m, t in zip(envelopes, translations)]
svc_solid = make_solid(pv - origin, pf)[0]
assert all(float((x ^ y).volume()) <= 0. for x, y in [(actual_solids[0], actual_solids[1]), (actual_solids[0], svc_solid), (actual_solids[1], svc_solid)]), 'Bounded protected artery search has no disjoint candidate'
transforms = []
for translation in translations:
    matrix, metrics = constrained_translation(translation, [], bounds=[[-.008, .008], [0., 0.], [-.008, .008]])
    transforms.append(matrix)
placed_envelopes = [transform_mesh(m, t) for m, t in zip(envelopes, transforms)]
separations = []
for first, second in [(placed_envelopes[0], placed_envelopes[1]), (placed_envelopes[0], svc_envelope), (placed_envelopes[1], svc_envelope)]:
    _, separation = pair_contacts(first, second, .0003, max_pairs=500000, max_tests=4000000)
    assert separation['complete_within_pair_budget']
    if separation['returned_triangle_pairs']:
        _,exact=pair_contacts(first,second,0.,max_pairs=500000,max_tests=4000000)
        assert exact['complete_within_pair_budget'] and not exact['crossing_face_pairs'] and not exact['coplanar_overlap_face_pairs'], 'Actual vessel envelopes cross'
        separation['zero_distance_check']=exact
    separation['three_tenths_mm_clearance_established']=separation['returned_triangle_pairs']==0
    separations.append(separation)
atomic_json(out / 'placement.json', dict(passed=True, dependency_groups=groups,
    transforms=[m.tolist() for m in transforms], bounded_search=trials,
    cusp_inlet_frames=[dict(transform=frame.tolist(),registration=proof) for frame,proof in inlet_frames],
    fixed_parts=['ZA:Superior vena cava', 'all cardiac wall cells', 'skeleton', 'lungs', 'diaphragm'],
    complete_root_separation=separations, elapsed_seconds=time.monotonic() - start))
print('Maximum root-group translation', round(distance * 1000, 3), 'mm; all three envelopes have zero occupied overlap and no crossing faces', flush=True)
meshes = {}
for ids, (_, root), matrix, (inlet,_) in zip(groups, roots, transforms,inlet_frames):
    meshes[ids[0]] = transform_mesh(root, matrix)
    for pid in ids[1:]:
        source = reader.read_part(atlas, parts[pid], frame=root.frame, catalogue_path=path)
        meshes[pid] = transform_mesh(source, matrix@inlet)
points, faces, fields, textures, cells, roles = [], {}, {}, {}, {}, {}
normal_repairs = []
count = 0
for pid, m in meshes.items():
    sid = pid + '/source-outflow-chart'
    points.extend(m.vertices - origin)
    local_faces=m.faces+count
    count += len(m.vertices)
    tri = m.triangles
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None]
    if 'normals' in m.attributes:
        words = m.attributes['normals'][m.faces].copy()
        if 'normalSigns' in m.attributes:
            words *= m.attributes['normalSigns'][m.faces, None]
        bad = np.einsum('tki,ti->tk', words.astype(float), n) <= 0
        words[bad] = np.broadcast_to(np.rint(n * 32767).astype('<i2')[:, None], words.shape)[bad]
        normal_repairs.append(dict(part=pid, source_nonoutward_corners=int(bad.sum())))
    else:
        words = np.repeat(np.rint(n * 32767).astype('<i2')[:, None], 3, axis=1)
    uv=m.attributes['uvs'][m.faces].copy() if 'uvs' in m.attributes else np.zeros((len(m.faces), 3, 2), '<f4')
    cells[pid]=[]
    spans=[('',0,len(m.faces),'')]
    if pid in root_face_counts:
        n=root_face_counts[pid];spans=[('/outer',0,n,''),('/inner',n,2*n,'arterial-inner:'+pid),('/rims',2*n,len(m.faces),'')]
    for suffix,first,last,role in spans:
        name=sid+suffix;faces[name]=local_faces[first:last];fields[name]=words[first:last];textures[name]=uv[first:last];roles[name]=role;cells[pid].append(dict(surface=name,side=1))
positions, exact = np.unique(np.array(points), axis=0, return_inverse=True)
faces = {sid: exact[fs] for sid, fs in faces.items()}
atomic_npz(out / 'source-arrangement.npz', positions=positions)
steps = []
def own(hosts, cutters):
    global positions, faces, fields, textures, cells, roles
    if time.monotonic() - start > 55:
        raise TimeoutError('Outflow ownership stage deadline')
    result = subtract_assembly(positions, faces, cells, hosts, cutters,exact_source_replay=not args.native_intersections)
    result,roles=split_by_source_role(result,roles)
    fields, textures = replay_corner_fields(positions, faces, fields, textures, result)
    positions, faces, cells = result['positions'], result['surfaces'], result['cells']
    steps.append(dict(hosts=hosts, cutters=cutters, evidence=result['evidence']))
    atomic_json(out / 'progress.json', dict(steps=steps, elapsed_seconds=time.monotonic() - start), overwrite=True)
# Preserve the source union, with explicit computational ownership at native
# coaptation. A shared cusp chart is not a fused anatomical attachment.
for ids in groups:
    own([ids[2]], [ids[1]])
    own(ids[1:3], [ids[3]])
own(groups[0][1:] + groups[1][1:], [ids[0] for ids in groups])
arrays, definitions = dict(positions=positions), []
for i, (sid, ff) in enumerate(faces.items()):
    arrays[f'surface_{i}_faces'] = ff
    arrays[f'surface_{i}_normals'] = fields[sid]
    arrays[f'surface_{i}_uvs'] = textures[sid]
    definitions.append(dict(id=sid, faces=f'surface_{i}_faces', normals=f'surface_{i}_normals', uvs=f'surface_{i}_uvs',role=roles[sid]))
atomic_npz(out / 'arrangement.npz', **arrays)
atomic_json(out / 'ownership.json', dict(
    constructed=True, catalogue_sha256=digest_file(path), origin=origin.tolist(), source_stage=args.stage,
    cells=cells, surfaces=definitions, source_parts={pid: parts[pid] for pid in meshes},
    dependency_groups=groups, translations=[m[:3, 3].tolist() for m in transforms],
    root_stages=[args.aorta, args.pulmonary], separation_search=trials,
    cusp_inlet_frames=[dict(transform=frame.tolist(),registration=proof) for frame,proof in inlet_frames],
    complete_root_separation=separations, source_normal_repairs=normal_repairs, ownership_steps=steps,
    numerical_acceptance_pending=True, host_attachments_pending=True, exact_source_replay=not args.native_intersections,
    elapsed_seconds=time.monotonic() - start, installed=False))
print(json.dumps(dict(constructed=True, cells=len(cells), individual_translation_mm=distance * 1000,
                      elapsed_seconds=time.monotonic() - start)), flush=True)
