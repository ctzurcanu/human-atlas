"""Encode source-bound muscle partitions and their reciprocal old wall cells.

Uses an already saved Boolean arrangement. No new anatomical inference,
deformation, installation or reference-geometry hiding occurs here.
"""
from collections import defaultdict
from pathlib import Path
import json
import numpy as np

from .common import atomic_json, atomic_npz, digest_file
from .model import AtlasReader, Mesh
from .spatial import BVH
from .tissue import pack
from .urinary_delivery import shared_encoding
from .algorithms.finite_tissue_boolean import barycentric
from .algorithms.quantized_surface_topology import conform_stored_surfaces
from .algorithms.float32_diagonal_repair import repair_diagonals


def construct(arrangement_directory, muscle_directory, output):
    folder, muscle_directory, output = map(lambda p: Path(p).resolve(), (arrangement_directory, muscle_directory, output))
    if output.exists():
        raise FileExistsError(output)
    spec = json.loads((folder/'arrangement.json').read_text())
    if digest_file(spec['source_catalogue']) != spec['source_catalogue_sha256']:
        raise ValueError('The installed source catalogue changed')
    atlas = json.loads(Path(spec['source_catalogue']).read_text())
    byid = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    originals = {pid: reader.read_part(atlas, byid[pid], frame=spec['frame']) for pid in spec['source_part_ids']}
    records = {pid: byid[pid] for pid in originals}
    parent_for = {pid: pid for pid in originals}
    host = 'BLADDER:male:neck:adventitia'
    new = 'BLADDER:male:neck:vesicoprostatic'
    parent_for[new] = host
    chart_for = {pid: ('prostate' if pid == 'ZA:Prostate' else 'bladder') for pid in spec['cells']}
    origins = {chart: np.asarray(byid[pid]['geometryOrigin'], float) for chart, pid in [('bladder', host), ('prostate', 'ZA:Prostate')]}
    for name in ['pubovesical-left', 'pubovesical-right']:
        pid = 'BLADDER:male:neck:'+name
        path = muscle_directory/name/('BLADDER_male_'+name+'.npz')
        originals[pid] = Mesh.load(path)
        meta = json.loads((path.parent/'candidate.json').read_text())
        rec = meta['packed']['parts'][meta['part_id']]
        rec = json.loads(json.dumps(rec))
        rec['tissueDrawGroups'][0]['neighbor'] = meta['first_part_id']
        rec['tissueDrawGroups'][-1]['neighbor'] = meta['second_part_id']
        records[pid] = rec
        parent_for[pid] = pid
        chart_for[pid] = name
        origins[name] = np.asarray(meta['geometry_origin'], float)
    with np.load(folder/'arrangement.npz', allow_pickle=False) as z:
        p = z['positions'].copy()
        original_surfaces = {sid: z[key].copy() for sid, key in spec['surfaces'].items()}
    owners = defaultdict(list)
    for pid, uses in spec['cells'].items():
        for use in uses:
            owners[use['surface']].append((pid, use['side']))
    source_lookups = {}
    face_neighbors = {}
    trees = {}
    for pid, mesh in originals.items():
        source_lookups[pid] = {tuple(sorted(map(tuple, t))): i for i, t in enumerate(mesh.triangles)}
        neighbors = np.full(len(mesh.faces), None, dtype=object)
        for group in records[pid]['tissueDrawGroups']:
            neighbors[group['start']//3:(group['start']+group['count'])//3] = group.get('neighbor')
        face_neighbors[pid] = neighbors
    def parents_for(pid, triangles):
        mesh = originals[pid]
        result = np.array([source_lookups[pid].get(tuple(sorted(map(tuple, t))), -1) for t in triangles], int)
        missing = np.flatnonzero(result < 0)
        if len(missing):
            if pid not in trees:
                trees[pid] = BVH(mesh)
            query, _ = trees[pid].closest_points(triangles[missing].mean(1))
            if query['distance_m'].max() > max(4e-12, spec.get('construction_resolution_m', 0)):
                raise ValueError('A proposed surface has no source parent: '+pid)
            result[missing] = query['face_ids']
        return result
    surfaces, cells, lineage = {}, {pid: [] for pid in spec['cells']}, {}
    for sid, fs in original_surfaces.items():
        uses = owners[sid]
        first, sign = uses[0]
        assert sign == 1
        if len(uses) == 2 and {who for who, _ in uses} == {host, new}:
            surfaces[sid] = fs
            lineage[sid] = None
            for who, side in uses:
                cells[who].append(dict(surface=sid, side=side, neighbor=next(a for a, _ in uses if a != who)))
            continue
        pubic_owner = next((who for who, _ in uses if ':pubovesical-' in who), None)
        parent = pubic_owner or parent_for[first]
        parents = parents_for(parent, p[fs])
        if len(uses) == 1:
            neighbors = face_neighbors[parent][parents]
            for peer in sorted(set(neighbors), key=str):
                if peer in cells:
                    raise ValueError('A source peer became an unpaired free surface: '+str((first, peer)))
                take = np.flatnonzero(neighbors == peer)
                name = sid+':'+str(peer)
                surfaces[name] = fs[take]
                lineage[name] = dict(part=parent, faces=parents[take])
                cells[first].append(dict(surface=name, side=1, neighbor=peer))
        else:
            surfaces[sid] = fs
            lineage[sid] = dict(part=parent, faces=parents)
            for who, side in uses:
                cells[who].append(dict(surface=sid, side=side, neighbor=next(a for a, _ in uses if a != who)))
    encoded, encoding = shared_encoding(p, surfaces, cells, origins, chart_for=chart_for)
    surfaces, kept, conformity = conform_stored_surfaces(encoded, surfaces, source_positions=p)
    for sid, source in lineage.items():
        if source is not None:
            source['faces'] = source['faces'][kept[sid]]
    # Preserve continuous source planes. Only a coplanar internal diagonal may
    # change, with the whole directed boundary and all controls retained.
    diagonal_repairs = {}
    for sid, fs in surfaces.items():
        source = lineage[sid]
        tri, stored = p[fs], encoded[fs]
        normal = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        stored_normal = np.cross(stored[:, 1]-stored[:, 0], stored[:, 2]-stored[:, 0])
        if np.all(np.einsum('ij,ij->i', normal, stored_normal) > 0):
            continue
        first = next(pid for pid, uses in cells.items() if any(u['surface'] == sid for u in uses))
        chart = chart_for[first]
        surfaces[sid], diagonal_repairs[sid] = repair_diagonals(p, fs, origins[chart], stored_positions=encoded)
        if source is not None:
            source['faces'] = parents_for(source['part'], p[surfaces[sid]])
    words, uvs, maximum_parent_plane_error = {}, {}, 0.
    prostatic_surfaces = {use['surface'] for use in cells['ZA:Prostate']}
    for sid, fs in surfaces.items():
        source = lineage[sid]
        if source is None:
            tri = encoded[fs]
            normal = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
            normal /= np.linalg.norm(normal, axis=1)[:, None]
            words[sid] = np.repeat(np.rint(normal*32767).astype('<i2')[:, None, :], 3, axis=1)
        else:
            mesh, parents = originals[source['part']], source['faces']
            src = mesh.triangles[parents]
            original_normal = np.cross(src[:, 1]-src[:, 0], src[:, 2]-src[:, 0])
            actual = p[fs]
            actual_normal = np.cross(actual[:, 1]-actual[:, 0], actual[:, 2]-actual[:, 0])
            orientation = np.where(np.einsum('ij,ij->i', original_normal, actual_normal) > 0, 1, -1)
            weights = np.array([barycentric(t, parent) for t, parent in zip(p[fs], src)])
            replay = np.einsum('tck,tkj->tcj', weights, src)
            maximum_parent_plane_error = max(maximum_parent_plane_error, float(np.linalg.norm(replay-p[fs], axis=2).max()))
            field = mesh.attributes['normals'].astype(int)*mesh.attributes.get('normalSigns', np.ones(len(mesh.vertices), int))[:, None]
            values = np.einsum('tck,tkj->tcj', weights, field[mesh.faces[parents]].astype(float))
            values /= np.linalg.norm(values, axis=2)[:, :, None]
            values = np.rint(values*32767).astype('<i2')
            for i, (tri, parent) in enumerate(zip(p[fs], src)):
                for k, point in enumerate(tri):
                    same = np.flatnonzero(np.all(parent == point, axis=1))
                    if len(same):
                        values[i, k] = field[mesh.faces[parents[i], same[0]]]
            values = (values.astype(int)*orientation[:, None, None]).astype('<i2')
            words[sid] = values
        if sid in prostatic_surfaces:
            mesh = originals['ZA:Prostate']
            parents = parents_for('ZA:Prostate', p[fs])
            weights = np.array([barycentric(t, parent) for t, parent in zip(p[fs], mesh.triangles[parents])])
            uvs[sid] = np.einsum('tck,tkj->tcj', weights, mesh.attributes['uvs'][mesh.faces[parents]])
    output.mkdir(parents=True)
    arrays, surface_keys, word_keys = {'positions': p, 'encoded_positions': encoded}, {}, {}
    for i, sid in enumerate(surfaces):
        surface_keys[sid] = 'surface'+str(i)
        word_keys[sid] = 'words'+str(i)
        arrays[surface_keys[sid]] = surfaces[sid]
        arrays[word_keys[sid]] = words[sid]
    atomic_npz(output/'encoded-arrangement.npz', **arrays)
    groups = {}
    for chart, origin in origins.items():
        local = {pid: [dict(u, neighbor=u.get('neighbor') if u.get('neighbor') in chart_for and chart_for[u['neighbor']] == chart else None) for u in uses]
                 for pid, uses in cells.items() if chart_for[pid] == chart}
        used = {u['surface'] for uses in local.values() for u in uses}
        result = pack(encoded, {sid: surfaces[sid] for sid in used}, local, {pid: pid for pid in local}, spec['frame'], output/chart,
                      origin=origin, corner_words={sid: words[sid] for sid in used},
                      corner_uvs={sid: uvs[sid] for sid in used} if chart == 'prostate' else None, max_pairs=1000000)
        groups[chart] = result
        for pid, part in result['parts'].items():
            for draw, use in zip(part['tissueDrawGroups'], cells[pid]):
                if use.get('neighbor'):
                    draw['neighbor'] = use['neighbor']
        atomic_json(output/chart/'delivery.json', result)
        print('chart', chart, 'passed', result['passed'], 'cells', len(local), flush=True)
    report = dict(schema='bladder-muscle-delivery/v1', source_catalogue=spec['source_catalogue'], source_catalogue_sha256=spec['source_catalogue_sha256'],
                  arrangement=str(folder/'arrangement.npz'), arrangement_sha256=digest_file(folder/'arrangement.npz'), frame=spec['frame'],
                  source=str(output/'encoded-arrangement.npz'), cells=cells, surfaces=surface_keys, corner_words=word_keys,
                  chart_for=chart_for, origins={k: v.tolist() for k, v in origins.items()}, encoding=encoding, conformity=conformity,
                  diagonal_repairs=diagonal_repairs, maximum_parent_plane_error_m=maximum_parent_plane_error,
                  groups=groups, passed=all(g['passed'] for g in groups.values()), anatomical_acceptance=False, installed=False)
    atomic_json(output/'delivery.json', report)
    return report
