"""Source-bound posterior bladder-muscle partition, without installation.

The classifier is an explicit illustrative construction input. Anatomical
identity, dimensions and rendered acceptance remain separate from CSG checks.
"""
from pathlib import Path
import json
import numpy as np
import manifold3d as md

from .common import atomic_json, atomic_npz, digest_file
from .model import Mesh
from .tissue import pack
from .algorithms.finite_cell_partition import partition_closed_cell
from .algorithms.finite_tissue_boolean import barycentric


def construct(source_manifest, output, bounds):
    source_manifest, output = Path(source_manifest).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    source = json.loads(source_manifest.read_text())
    if digest_file(source['source_catalogue']) != source['source_catalogue_sha256']:
        raise ValueError('The source urinary assembly changed')
    host_id = 'BLADDER:male:neck:adventitia'
    new_id = 'BLADDER:male:neck:vesicoprostatic'
    row = next(r for r in source['sources'] if r['part_id'] == host_id)
    if digest_file(row['path']) != row['sha256']:
        raise ValueError('The finite source compartment changed')
    host = Mesh.load(row['path'])
    low, high = np.asarray(bounds, float)
    if low.shape != (3,) or high.shape != (3,) or not np.isfinite(bounds).all() or np.any(high <= low):
        raise ValueError('Explicit finite classifier bounds required')
    classifier = md.Manifold.cube((high-low).tolist()).translate(low.tolist()).to_mesh64()
    result = partition_closed_cell(host.vertices, host.faces,
                                  np.asarray(classifier.vert_properties)[:, :3],
                                  np.asarray(classifier.tri_verts))
    positions = result['positions']
    groups = row['record']['tissueDrawGroups']
    normal_words = host.attributes['normals'].astype(int)*host.attributes.get('normalSigns', np.ones(len(host.vertices), int))[:, None]
    source_has_uvs = 'uvs' in host.attributes
    surfaces, words, uvs, neighbors = {}, {}, {}, {}
    cells = {'inside': [], 'outside': []}
    contact_areas = {}
    for owner, (faces, parents) in result['source_fragments'].items():
        for ordinal, group in enumerate(groups):
            take = (parents >= group['start']//3) & (parents < (group['start']+group['count'])//3)
            if not take.any():
                continue
            sid = f'{owner}:source:{ordinal}'
            triangles, ancestry = faces[take], parents[take]
            surfaces[sid] = triangles
            neighbors[sid] = group.get('neighbor')
            # External peers are deliberately absent from this two-cell pack;
            # the source-bound neighbour ledger carries their required joins.
            cells[owner].append(dict(surface=sid, side=1))
            weights = np.array([barycentric(positions[t], host.triangles[k]) for t, k in zip(triangles, ancestry)])
            field = np.einsum('tck,tkj->tcj', weights, normal_words[host.faces[ancestry]].astype(float))
            field /= np.linalg.norm(field, axis=2)[:, :, None]
            field = np.rint(field*32767).astype('<i2')
            for i, (actual, parent) in enumerate(zip(positions[triangles], host.triangles[ancestry])):
                for j, point in enumerate(actual):
                    literal = np.flatnonzero(np.all(parent == point, axis=1))
                    if len(literal):
                        field[i, j] = normal_words[host.faces[ancestry[i], literal[0]]]
            words[sid] = field
            if source_has_uvs:
                uvs[sid] = np.einsum('tck,tkj->tcj', weights, host.attributes['uvs'][host.faces[ancestry]])
            if owner == 'inside':
                t = positions[triangles]
                contact_areas[sid] = float(np.linalg.norm(np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), axis=1).sum()/2)
    actual_peers = {neighbors[sid] for sid in contact_areas if neighbors[sid] is not None}
    required_peers = {'BLADDER:male:neck:neck-external-longitudinal', 'ZA:Prostate'}
    if actual_peers != required_peers:
        raise ValueError(f'Posterior muscle must touch only the declared tissue peers, found {actual_peers}')
    sid = 'inside:outside'
    cut = result['cut_faces']
    surfaces[sid] = cut
    t = positions[cut]
    field = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
    field /= np.linalg.norm(field, axis=1)[:, None]
    words[sid] = np.repeat(np.rint(field*32767).astype('<i2')[:, None, :], 3, axis=1)
    if source_has_uvs:
        uvs[sid] = np.zeros((len(cut), 3, 2))
    cells['inside'].append(dict(surface=sid, side=1, neighbor='outside'))
    cells['outside'].append(dict(surface=sid, side=-1, neighbor='inside'))
    output.mkdir(parents=True)
    packed = pack(positions, surfaces, cells, {'inside': new_id, 'outside': host_id}, source['frame'], output/'packed',
                  origin=row['record']['geometryOrigin'], corner_words=words,
                  corner_uvs=uvs if source_has_uvs else None, max_pairs=1000000)
    connected = all(len(check['topology']['components']) == 1 for check in packed['checks']['cells'].values())
    atomic_npz(output/'partition.npz', positions=positions, cut_faces=cut,
               inside_source_faces=result['source_fragments']['inside'][0],
               inside_source_parents=result['source_fragments']['inside'][1],
               outside_source_faces=result['source_fragments']['outside'][0],
               outside_source_parents=result['source_fragments']['outside'][1])
    report = dict(schema='bladder-posterior-muscle-partition/v1', source_manifest=str(source_manifest),
                  source_manifest_sha256=digest_file(source_manifest), source_catalogue_sha256=source['source_catalogue_sha256'],
                  source_host_id=host_id, new_muscle_id=new_id, ta98_term='A08.3.01.021M', classifier_bounds=[low.tolist(), high.tolist()],
                  source_neighbors=neighbors, contact_areas_m2=contact_areas, evidence=result['evidence'], packed=packed,
                  passed=bool(packed['passed'] and connected), both_cells_connected=connected,
                  adjoining_source_subdivisions_installed=False,
                  illustrative_dimensions=True, anatomical_acceptance=False, installed=False)
    atomic_json(output/'candidate.json', report)
    return report
