"""Source-rim urinary joins with one coherent distal ureter transport.

The existing renal roots and skeleton remain fixed. Only explicit outer and
inner mouth cycles are bridged; no triangle spans a urinary lumen.
"""
from pathlib import Path
import json

import numpy as np

from .algorithms.annular_pipe_envelopes import boundary_loops, face_components
from .algorithms.constrained_boundary_disk import cancel_opposed_faces
from .algorithms.tube_axis_map import axis_transport, endpoint_course
from .algorithms.tube_rim_bridge import bridge_boundaries, end_annulus
from .algorithms.tube_sections import arclength, path_tangents, unit
from .algorithms.tube_station_topology import recover_station_rings
from .audit import audit
from .bladder_carriers import recover
from .common import atomic_json, atomic_npz, digest_file
from .db import Database
from .model import Mesh
from .registry import load_registered
from .review import prepare


def source_tube(db, side, output, cache=None):
    loaded = [load_registered(db, f'RENALDETAIL:male:{side}:ureter-wall:{name}', 'local-ta98')
              for name in ['adventitia', 'muscularis', 'mucosa']]
    meshes = [q[1]for q in loaded]
    key=lambda points:tuple(sorted(tuple(p)for p in points))
    root_keys=set();mucosal_free_keys=set()
    for role,(record,mesh)in zip(['adventitia','muscularis','mucosa'],loaded):
        for g in record['metadata']['record']['tissueDrawGroups']:
            fs=mesh.faces[g['start']//3:(g['start']+g['count'])//3]
            if g.get('neighbor','').startswith(f'RENALDETAIL:male:{side}:upj:'):
                root_keys.update(key(mesh.vertices[t])for t in fs)
            elif role=='mucosa'and not g.get('neighbor'):
                mucosal_free_keys.update(key(mesh.vertices[t])for t in fs)
    vertices = np.concatenate([m.vertices for m in meshes])
    faces = np.concatenate([m.faces + sum(len(q.vertices) for q in meshes[:i]) for i, m in enumerate(meshes)])
    v, inv = np.unique(vertices, axis=0, return_inverse=True)
    f, _ = cancel_opposed_faces(inv[faces])
    used, inv = np.unique(f, return_inverse=True); v = v[used]; f = inv.reshape(-1, 3)
    wall = Mesh(v, f, meshes[0].frame)
    _, native = load_registered(db, 'ZA:Ureter.' + side[0], 'local-ta98')
    nv, inverse = np.unique(native.vertices, axis=0, return_inverse=True)
    rings = recover_station_rings(nv, inverse[native.faces])[0]['rings']
    axis = np.array([nv[q].mean(0) for q in rings])
    if axis[0, 1] > axis[-1, 1]:
        axis = axis[::-1]
    if cache is not None:
        cached=Path(cache);saved=Mesh.load(cached/(side+'-source-wall.npz'))
        method=json.loads((cached/(side+'-sheet-lineage.json')).read_text())
        if method.get('explicit_root_and_mucosal_lineage')is not True:raise ValueError('Cached radial-only source sheet selection is not sufficient')
        with np.load(cached/(side+'-axis.npz'),allow_pickle=False)as z:
            if not np.array_equal(z['axis'],axis):raise ValueError('Cached native axis changed')
        with np.load(cached/(side+'-carriers/carriers.npz'),allow_pickle=False)as z:
            if saved.geometry_sha!=wall.geometry_sha or not np.array_equal(saved.faces,wall.faces):raise ValueError('Cached actual source wall changed')
            inside=z['inner_source_face_ids'].copy()
        result=json.loads((cached/(side+'-carriers/carriers.json')).read_text())
        if result['source_geometry_sha256']!=wall.geometry_sha or not all(q['passed']for q in result['checks'].values()):raise ValueError('Cached source carrier checks do not pass')
        wall.save(output/(side+'-source-wall.npz'));atomic_npz(output/(side+'-axis.npz'),axis=axis,frame=np.array(wall.frame))
        return wall,axis,inside,result
    center = wall.triangles.mean(1); a = axis[:-1]; d = np.diff(axis, axis=0)
    u = np.clip(np.einsum('fki,ki->fk', center[:, None] - a, d) / np.sum(d * d, axis=1), 0, 1)
    q = a + u[:, :, None] * d
    nearest = q[np.arange(len(center)), np.argmin(np.sum((center[:, None] - q) ** 2, axis=2), axis=1)]
    n = np.cross(wall.triangles[:, 1] - wall.triangles[:, 0], wall.triangles[:, 2] - wall.triangles[:, 0])
    radial = center - nearest
    cosine = np.einsum('ij,ij->i', n, radial) / np.linalg.norm(n, axis=1) / np.linalg.norm(radial, axis=1)
    proposed = np.flatnonzero(cosine < -.25)
    _, labels = face_components(f[proposed]); inside = proposed[labels == np.argmax(np.bincount(labels))]
    # The native guide only locates the distal opening. At the transferred
    # renal root, use the ACTUAL reciprocal UPJ annuli and mucosal source chart
    # to separate the complete two physical sheets, without an axis guess.
    distal=min(boundary_loops(f[inside]),key=lambda r:v[r].mean(0)[1])
    distal_cap,_=end_annulus(f,f[inside],distal)
    root_ids=np.array([i for i,t in enumerate(f)if key(v[t])in root_keys],int)
    if len(root_ids)!=len(root_keys):raise ValueError('The complete stored renal-root annulus is missing')
    if len(boundary_loops(f[root_ids]))!=2:raise ValueError('The stored renal-root tissue interface is not an annulus')
    sides=np.setdiff1d(np.arange(len(f)),np.r_[distal_cap,root_ids])
    count,component=face_components(f[sides])
    if count!=2:raise ValueError('The two explicit end annuli do not separate exactly the inner and outer sheets')
    mucosal_labels={int(component[i])for i,face in enumerate(f[sides])if key(v[face])in mucosal_free_keys}
    if len(mucosal_labels)!=1:raise ValueError('The stored mucosa does not establish exactly one inner sheet')
    inside=sides[component==next(iter(mucosal_labels))]
    if len(boundary_loops(f[inside]))!=2:raise ValueError('The complete source lining has unexpected mouths')
    wall.save(output / (side + '-source-wall.npz'))
    result = recover(wall, output / (side + '-carriers'), inner_face_ids=inside,
                     mouth_labels=['distal-ureter', 'renal-root'], cap_forbid_source_diagonals=True)
    atomic_npz(output / (side + '-axis.npz'), axis=axis, frame=np.array(wall.frame))
    atomic_json(output/(side+'-sheet-lineage.json'),dict(explicit_root_and_mucosal_lineage=True,
        source_geometry_sha256=wall.geometry_sha,root_interface_faces=len(root_ids),distal_interface_faces=len(distal_cap),
        inner_faces=len(inside),radial_guide_used_only_for_distal_location=True))
    return wall, axis, inside, result


def curved_prefix(old_axis, tip, normal, length=.045):
    arc = arclength(old_axis); k = int(np.searchsorted(arc, length))
    if k < 3 or k >= len(old_axis) - 1:
        raise ValueError('A source-connected fixed ureter station is required')
    tangent = path_tangents(old_axis); end = old_axis[k]
    distance = float(np.linalg.norm(end - tip)); u = np.linspace(0, 1, 121)[:, None]
    first, last = distance * normal, distance * tangent[k]
    prefix = (2*u**3-3*u**2+1)*tip + (u**3-2*u**2+u)*first + (-2*u**3+3*u**2)*end + (u**3-u**2)*last
    derivative = (6*u**2-6*u)*tip + (3*u**2-4*u+1)*first + (-6*u**2+6*u)*end + (3*u**2-2*u)*last
    prefix[0] = tip; prefix[-1] = end
    new, rotations = endpoint_course(old_axis, prefix, unit(derivative), k)
    return new, rotations, k


def build(bladder_wall, bladder_carriers, output, *, source_cache=None):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    bladder = Mesh.load(bladder_wall)
    with np.load(bladder_carriers, allow_pickle=False) as z:
        if not np.array_equal(bladder.faces, z['wall_faces']):
            raise ValueError('The reviewed whole bladder source faces changed')
        lining = z['inner_side_faces'].copy()
        port_loops = {side: z[side + '_ureteric_ring'].copy() for side in ['left', 'right']}
    body_remove = []; body_rims = {}
    for side, inner in port_loops.items():
        removed, outer = end_annulus(bladder.faces, lining, inner)
        body_remove.extend(removed); body_rims[side] = (outer, inner)
    keep = np.setdiff1d(np.arange(len(bladder.faces)), body_remove)
    vertices = [bladder.vertices]; faces = [bladder.faces[keep]]
    inner_keys = {tuple(sorted(t)) for t in lining}
    inner_ids = [i for i, f in enumerate(faces[0]) if tuple(sorted(f)) in inner_keys]
    rows = []; db = Database(); db.assert_current_catalogues('local-ta98')
    for side in ['left', 'right']:
        wall, axis, inside, source = source_tube(db, side, output, source_cache)
        loops = boundary_loops(wall.faces[inside])
        distal = min(loops, key=lambda q: wall.vertices[q].mean(0)[1])
        removed, tube_outer = end_annulus(wall.faces, wall.faces[inside], distal)
        body_outer, body_inner = body_rims[side]
        xyz = bladder.vertices[body_outer]
        normal = -np.cross(xyz - xyz.mean(0), np.roll(xyz, -1, axis=0) - xyz.mean(0)).sum(0)
        normal /= np.linalg.norm(normal)
        tip = bladder.vertices[body_outer].mean(0) + .008 * normal
        course, rotations, fixed_station = curved_prefix(axis, tip, normal)
        moved, local_rotations = axis_transport(wall.vertices, axis, course, rotations)
        transported = Mesh(moved, wall.faces, wall.frame)
        check = audit(transported, require_closed=True, expected_components=1,
                      expected_genus=[1], max_pairs=1000000)
        transported.save(output / (side + '-transported-wall.npz'))
        if not check['passed']:
            atomic_json(output / (side + '-transport-failure.json'), check)
            raise ValueError('The shared distal ureter transport fails the actual wall checks')
        offset = sum(len(v) for v in vertices); face_start = sum(len(f) for f in faces)
        vertices.append(moved)
        retained = np.setdiff1d(np.arange(len(wall.faces)), removed)
        faces.append(wall.faces[retained] + offset)
        inner_set = set(map(int, inside))
        inner_ids.extend(face_start + i for i, original in enumerate(retained) if int(original) in inner_set)
        positions = np.concatenate(vertices)
        outer_bridge, outer_proof = bridge_boundaries(positions, tube_outer + offset, body_outer)
        faces.append(outer_bridge)
        inner_start = sum(len(f) for f in faces)
        inner_bridge, inner_proof = bridge_boundaries(positions, distal + offset, body_inner)
        faces.append(inner_bridge)
        inner_ids.extend(range(inner_start, inner_start + len(inner_bridge)))
        atomic_npz(output / (side + '-transport.npz'), source_axis=axis, target_axis=course,
                   rotations=rotations, per_control_rotations=local_rotations,
                   bladder_inner_ring=body_inner, bladder_outer_ring=body_outer,
                   ureter_inner_ring=distal+offset, ureter_outer_ring=tube_outer+offset,
                   outer_bridge=outer_bridge, inner_bridge=inner_bridge)
        rows.append(dict(side=side, source_geometry_sha256=wall.geometry_sha,
                         source_axis_sha256=digest_file(output / (side + '-axis.npz')),
                         fixed_from_station=fixed_station, fixed_tail_literal=np.array_equal(axis[fixed_station:], course[fixed_station:]),
                         tip_m=tip.tolist(), outward_normal=normal.tolist(),
                         maximum_displacement_m=float(np.linalg.norm(moved - wall.vertices, axis=1).max()),
                         transported_wall_check=check, outer_bridge=outer_proof, inner_bridge=inner_proof))
    whole = Mesh(np.concatenate(vertices), np.concatenate(faces), bladder.frame)
    whole.save(output / 'urinary-wall.npz')
    check = audit(whole, require_closed=True, expected_components=1, expected_genus=[2], max_pairs=1000000)
    result = dict(schema='urinary-explicit-rim-joins/v1', bladder_source_sha256=bladder.geometry_sha,
                  ureters=rows, whole_check=check, passed=check['passed'],
                  anatomical_acceptance=False, skeleton_modified=False, installed=False,
                  source_faces_removed_only_at_explicit_join_annuli=True,
                  lumen_disks_created=False, remaining=['Fluid carrier replay, full sections and pelvic contacts.',
                    'Coat continuity through the collars; source renal roots and bladder neck.',
                    'Bladder features, ordinary catalogue installation and female transfer.'])
    atomic_npz(output / 'inner-face-ids.npz', inner_face_ids=np.asarray(inner_ids, int))
    atomic_json(output / 'junctions.json', result)
    if check['passed']:
        result['carriers'] = recover(whole, output / 'joined-carriers', inner_face_ids=np.asarray(inner_ids, int),
            mouth_labels=['urethral-neck', 'right-renal-root', 'left-renal-root'], cap_forbid_source_diagonals=True)
        center = bladder.bounds.mean(0)
        result['review'] = prepare(db, dict(frame=whole.frame,
            meshes=[dict(mesh=str(output/'urinary-wall.npz'), name='Joined bladder and ureters', color=[.78,.45,.4,1])],
            views=[dict(name='front-urinary-joins',direction=[0,.1,1]),
                   dict(name='posterior-urinary-joins',direction=[0,.1,-1]),
                   dict(name='bladder-coronal-section',direction=[.1,.1,1],
                        section=dict(normal=[0,0,1],point=center.tolist(),clear_outer=True))], resolution=1000))
        atomic_json(output / 'junctions-reviewed-request.json', result)
    return result


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('bladder_wall')
    p.add_argument('--carriers',required=True);p.add_argument('--output',required=True)
    p.add_argument('--source-cache')
    args=p.parse_args();r=build(args.bladder_wall,args.carriers,args.output,source_cache=args.source_cache)
    print(json.dumps(dict(passed=r['passed'],failures=r['whole_check']['failures'],review=r.get('review')),indent=2))
