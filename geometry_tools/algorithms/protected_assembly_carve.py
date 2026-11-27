"""Carve declared host cells while preserving other components exactly.

A protected component can become exposed where its former neighbour is cut
away. Its old shared chart is then partitioned into a remaining contact and a
free surface. This is not a repair by plugging a hole or moving anatomy.
"""
from collections import defaultdict
import numpy as np
from .finite_tissue_boolean import subtract_occupied, barycentric, make_solid
from .source_face_complement import complement_source_faces, conform_surface_boundary_edges
from .shared_triangle_edges import reconcile_boundaries, refine_surfaces
from .finite_tissue_boolean import CommonPositionBank
from .conforming_fields import replay_corner_fields


def carve_with_protected_components(positions, surfaces, cells, normals, uvs,
        host_names, cutter_positions, cutter_faces, *, cutter_owners=None,
        classifier_cap_faces=(), tolerance=2e-12, prepared_carve=None,
        save_carve=None, boundary_failure=None, exact_source_replay=True,
        require_protected_charts_literal=False,
        prove_cap_exclusion_by_output_ancestry=False):
    p = np.asarray(positions, float)
    hosts = set(host_names)
    if not hosts or not hosts < set(cells):
        raise ValueError('Declare nonempty host and protected component sets')
    protected = set(cells) - hosts
    owners = defaultdict(list)
    for pid, uses in cells.items():
        for use in uses:
            owners[use['surface']].append((pid, use))
    selected = {sid: fs for sid, fs in surfaces.items() if any(pid in hosts for pid, _ in owners[sid])}
    detached = {}
    for pid in hosts:
        detached[pid] = [dict(u, **({'neighbor':None} if u.get('neighbor') in protected else {})) for u in cells[pid]]
    cell_faces = {pid: np.concatenate([selected[u['surface']] if u['side'] == 1 else selected[u['surface']][:, [0, 2, 1]] for u in uses]) for pid, uses in detached.items()}
    parents = {pid: [(u['surface'], i, u['side']) for u in uses for i in range(len(selected[u['surface']]))] for pid, uses in detached.items()}
    rows = {sid: dict(id=sid, kind='shared' if sum(pid in hosts for pid, _ in owners[sid]) == 2 else 'outer', array=sid) for sid in selected}
    result = prepared_carve
    if result is None:
        result = subtract_occupied(p, selected, rows, cell_faces, parents, detached,
            cutter_positions, cutter_faces,
            classifier_cap_faces=() if prove_cap_exclusion_by_output_ancestry else classifier_cap_faces,
            occupied_face_owners=cutter_owners, exact_source_replay=exact_source_replay,
            reconcile_shared_subdivisions=True, intersection_tolerance=tolerance)
        if save_carve: save_carve(result)
    # A cap may touch an existing reciprocal wall boundary. In that explicit
    # mode require the stronger, result-wide exclusion certificate: no actual
    # cut surface of any host may descend from any virtual cap triangle.
    # Reject even a tiny surviving cap fragment; never hide it after the fact.
    cap_ids = np.asarray(classifier_cap_faces, dtype=np.int64)
    if prove_cap_exclusion_by_output_ancestry:
        for row in result['surfaceRows']:
            if row.get('parentSurface') == 'occupied-outer':
                if np.any(np.isin(result['charts'][row['id']]['parentFaces'], cap_ids)):
                    raise ValueError('Virtual cap contributed to an anatomical cut: ' + row['id'])
    # The computational cutter may itself contain almost coincident Boolean
    # controls. Only its/intersection controls may be identified within the
    # declared construction tolerance; the native position prefix is literal.
    rawq=result['positions']
    canonical_bank=CommonPositionBank(p,tolerance)
    remap=np.r_[np.arange(len(p)),np.array([canonical_bank.insert(point) for point in rawq[len(p):]],int)]
    q=canonical_bank.positions
    result=dict(result)
    result['surfaces']=dict(result['surfaces']); result['charts']=dict(result['charts'])
    collapsed=0
    for sid,faces in result['surfaces'].items():
        faces=remap[faces]; tri=q[faces]
        retain=np.any(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])!=0,axis=1)
        collapsed+=int((~retain).sum())
        result['surfaces'][sid]=faces[retain]
        result['charts'][sid]={key:np.asarray(value)[retain] for key,value in result['charts'][sid].items()}
    output = dict(result['surfaces'])
    uses = {pid: [dict(u) for u in value] for pid, value in result['cells'].items()}
    source = dict(selected)
    source['COMPUTE:carve-boundary'] = np.asarray(cutter_faces) + len(p)
    fields = dict(normals)
    textures = dict(uvs)
    tri = np.asarray(cutter_positions)[cutter_faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None]
    fields['COMPUTE:carve-boundary'] = np.repeat(np.rint(n * 32767).astype('<i2')[:, None], 3, axis=1)
    textures['COMPUTE:carve-boundary'] = np.zeros((len(tri), 3, 2), '<f4')
    charts = {}
    for row in result['surfaceRows']:
        sid = row['id']
        family = sid if sid in selected else 'COMPUTE:carve-boundary'
        chart = result['charts'][sid]
        charts[sid] = dict(source_surfaces=[family] * len(output[sid]), source_faces=chart['parentFaces'],
            normal_signs=np.full(len(output[sid]), row.get('parentNormalSign', 1)), barycentrics=chart['barycentrics'])
    words, tex = replay_corner_fields(np.r_[p, cutter_positions], source, fields, textures,
                                     dict(positions=q, surfaces=output, charts=charts))
    # Computational voids are never installed as solids or anatomical neighbours.
    for value in uses.values():
        for u in value:
            if u.get('neighbor') not in hosts:
                u.pop('neighbor', None)
    exposed = []
    requests = []
    for pid in protected:
        uses[pid] = []
        for old in cells[pid]:
            sid = old['surface']
            wall = next((name for name, _ in owners[sid] if name in hosts), None)
            if wall is None:
                output[sid], words[sid], tex[sid] = surfaces[sid], normals[sid], uvs[sid]
                uses[pid].append(dict(old))
                continue
            kept = output[sid]
            par = result['charts'][sid]['parentFaces']
            counts = np.bincount(par, minlength=len(surfaces[sid]))
            changed = counts != 1
            for parent in np.flatnonzero(counts == 1):
                face = kept[np.flatnonzero(par == parent)[0]]
                original = p[surfaces[sid][parent]]
                if not any(np.array_equal(q[face], np.roll(original, k, axis=0)) for k in range(3)):
                    changed[parent] = True
            affected = np.flatnonzero(changed)
            if require_protected_charts_literal and len(affected):
                raise ValueError('Attachment guard failed to retain a whole source chart: '+sid)
            if len(kept):
                uses[pid].append(dict(old))
                for u in uses[wall]:
                    if u['surface'] == sid: u['neighbor'] = pid
            if not len(affected): continue
            lookup = np.full(len(surfaces[sid]), -1, int)
            lookup[affected] = np.arange(len(affected))
            take = np.flatnonzero(changed[par])
            q, free, free_parent, proof = complement_source_faces(q, p, surfaces[sid][affected], kept[take], lookup[par[take]], tolerance=tolerance, planar_scale=1., backend='triangle')
            if not len(free): continue
            free_parent = affected[free_parent]
            name = f'{sid}/exposed:{pid}'
            output[name] = free
            pseudo = dict(positions=q, surfaces={name:free}, charts={name:dict(
                source_surfaces=[sid] * len(free), source_faces=free_parent, normal_signs=np.ones(len(free), int),
                barycentrics=np.array([barycentric(q[face], p[surfaces[sid][parent]]) for face, parent in zip(free, free_parent)]))})
            ww, tt = replay_corner_fields(p, surfaces, normals, uvs, pseudo)
            words.update(ww); tex.update(tt)
            uses[pid].append(dict(surface=name, side=old['side']))
            exposed.append(dict(protected_part=pid, former_neighbor=wall, surface=name, proof=proof))
            for parent in affected:
                requests.append((sid, par.copy(), name, free_parent.copy(), int(parent)))
    # Canonicalise only literally equal coordinates, retaining the source bank
    # prefix and all original field words for subsequent bounded storage checks.
    _, first, inverse = np.unique(q, axis=0, return_index=True, return_inverse=True)
    canonical = first[inverse]
    output = {sid:canonical[fs] for sid, fs in output.items() if len(fs)}
    for pid in uses: uses[pid] = [u for u in uses[pid] if u['surface'] in output]
    # Boolean contact and planar complement triangulations may subdivide the
    # same edge differently. Reconcile actual opposite boundary chains first;
    # the parent-patch test below then checks the complete original perimeter.
    conformed, conformity = conform_surface_boundary_edges(q, output, uses,
        tolerance=tolerance, record_face_parents=True)
    ancestry0 = conformity.pop('sourceFaceParents')
    for sid, fs in conformed.items():
        if np.array_equal(fs, output[sid]): continue
        parent = ancestry0[sid]
        pseudo = dict(positions=q, surfaces={sid:fs}, charts={sid:dict(
            source_surfaces=[sid] * len(fs), source_faces=parent, normal_signs=np.ones(len(fs), int),
            barycentrics=np.array([barycentric(q[face], q[output[sid][i]]) for face, i in zip(fs, parent)]))})
        ww, tt = replay_corner_fields(q, {sid:output[sid]}, {sid:words[sid]}, {sid:tex[sid]}, pseudo)
        words[sid], tex[sid] = ww[sid], tt[sid]
    output = conformed
    plans = defaultdict(set)
    maximum = 0.
    for sid, kept_parent, free_sid, free_parent, parent in requests:
        kept = output[sid][kept_parent[ancestry0[sid]] == parent] if sid in output else np.empty((0,3),int)
        free = output[free_sid][free_parent[ancestry0[free_sid]] == parent]
        actual = np.r_[kept, free]
        expected=canonical[surfaces[sid][parent:parent+1]]
        try:
            maximum = max(maximum, reconcile_boundaries(q, actual, expected, plans, tolerance=tolerance))
        except ValueError:
            if boundary_failure: boundary_failure(q, actual, expected, dict(surface=sid,parent=parent))
            raise
    # Parent-patch validation also proposes subdivisions on virtual original
    # diagonals. The earlier whole-cell conformity may already have propagated
    # every real edge. Do not retessellate an already closed real arrangement
    # merely to refine a non-rendered source comparison triangle.
    closed=True
    for pid,value in uses.items():
        fs=np.concatenate([output[u['surface']] if u['side']==1 else output[u['surface']][:,[0,2,1]] for u in value])
        edges=np.r_[fs[:,[0,1]],fs[:,[1,2]],fs[:,[2,0]]]
        _,inv,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
        balance=np.bincount(inv,weights=np.where(edges[:,0]<edges[:,1],1,-1))
        if np.any(count!=2) or np.any(balance):
            closed=False
            if boundary_failure:
                bad=np.flatnonzero((count!=2)|(balance!=0))
                ids=np.unique(np.flatnonzero(np.isin(inv,bad))%len(fs))
                boundary_failure(q,fs[ids],np.empty((0,3),int),dict(part=pid,
                    nonmanifold_or_unbalanced_edges=len(bad),minimum_count=int(count.min()),
                    maximum_count=int(count.max()),phase='whole-cell-conformity'))
            break
    if closed:
        refined=output;ancestry={sid:np.arange(len(fs)) for sid,fs in output.items()};changed=0
    else:
        bank = CommonPositionBank(q, 0.)
        refined, ancestry, changed = refine_surfaces(bank, output, plans)
        q = bank.positions
    for sid, fs in refined.items():
        if np.array_equal(fs, output[sid]): continue
        parent = ancestry[sid]
        pseudo = dict(positions=q, surfaces={sid:fs}, charts={sid:dict(
            source_surfaces=[sid] * len(fs), source_faces=parent, normal_signs=np.ones(len(fs), int),
            barycentrics=np.array([barycentric(q[face], q[output[sid][i]]) for face, i in zip(fs, parent)]))})
        ww, tt = replay_corner_fields(q, {sid:output[sid]}, {sid:words[sid]}, {sid:tex[sid]}, pseudo)
        words[sid], tex[sid] = ww[sid], tt[sid]
    # The protected union is preserved facewise by complement and refinement;
    # this independent volume check catches a reversed or omitted chart.
    volume_proof = []
    for pid in sorted(protected):
        before = np.concatenate([surfaces[u['surface']] if u['side'] == 1 else surfaces[u['surface']][:, [0,2,1]] for u in cells[pid]])
        after = np.concatenate([refined[u['surface']] if u['side'] == 1 else refined[u['surface']][:, [0,2,1]] for u in uses[pid]])
        a, _ = make_solid(p, before); b, _ = make_solid(q, after)
        error = abs(float(a.volume()) - float(b.volume()))
        if error > max(1e-19, abs(float(a.volume())) * 1e-10):
            raise ValueError('Protected component volume changed: ' + pid)
        volume_proof.append(dict(part=pid, volume_error_m3=error))
    return dict(positions=q, surfaces=refined, cells=uses, normals=words, uvs=tex,
        evidence=dict(carve=result['evidence'], protected_exposures=exposed,
            computational_control_canonicalization_m=canonical_bank.maximum_snap,
            exact_source_intersection_replay=exact_source_replay,
            protected_charts_required_literal=require_protected_charts_literal,
            exact_zero_computational_faces_removed=collapsed,
            source_prefix_literal=np.array_equal(q[:len(p)], p), source_geometry_not_moved=True,
            protected_volume_checks=volume_proof, boundary_edges_refined=len(plans), refined_faces=changed,
            maximum_boundary_refinement_error_m=maximum, computational_fluid_cells_emitted=0,
            protected_boundary_conformity=conformity,
            cap_exclusion_by_complete_output_ancestry=prove_cap_exclusion_by_output_ancestry,
            classifier_cap_contributions=0,
            already_conforming_real_cells=closed,
            geometry_and_anatomical_audits_pending=True))
