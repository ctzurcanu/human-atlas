"""Read literal native pulmonary curves; do not change an atlas or build tissue."""
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader


def main():
    began = time.monotonic()
    source = Path('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')
    census_path = Path('reports/ta98-blood-circulation-inventory-installed-v12.json')
    census = json.loads(census_path.read_text())
    source_sha = digest_file(source)
    assert source_sha == census['source_catalogue_sha256']
    atlas = json.loads(source.read_text())
    parts = {p['id']: p for p in atlas['parts']}
    hidden = set(atlas['defaultHidden'])
    terms = [t for t in census['terms'] if t['id'].startswith(('A12.2.01.', 'A12.3.02.'))]
    reader = AtlasReader()
    candidate_ids = sorted({p['id'] for p in atlas['parts'] if p['id'].startswith('ZA:') and
                            ('lung' in p['name'].lower() and p['system'] in ('arterial', 'venous') or
                             'pulmonary vein' in p['name'].lower() or
                             p['id'] in ('ZA:Right pulmonary artery', 'ZA:Left pulmonary artery'))})
    out = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
    out.mkdir(parents=True, exist_ok=True)
    assert not (out / 'evidence.json').exists(), 'Preserve previous evidence'
    records = []
    arrays = {}
    for pi, ident in enumerate(candidate_ids):
        p = parts[ident]
        mesh = reader.read_part(atlas, p, frame='model:male:assembled-metres')
        v, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
        f = inverse[mesh.faces]
        edges = np.sort(np.vstack((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]])), axis=1)
        unique_edges, incidence = np.unique(edges, axis=0, return_counts=True)
        graph = coo_matrix((np.ones(len(unique_edges)), (unique_edges[:, 0], unique_edges[:, 1])), shape=(len(v), len(v)))
        component_count, labels = connected_components(graph, directed=False)
        boundary = unique_edges[incidence == 1]
        components = []
        for ci in range(component_count):
            face_ids = np.flatnonzero(labels[f[:, 0]] == ci)
            used = np.unique(f[face_ids])
            b = boundary[labels[boundary[:, 0]] == ci]
            rings = []
            if len(b):
                bv, bi = np.unique(b, return_inverse=True)
                be = bi.reshape(-1, 2)
                bg = coo_matrix((np.ones(len(be)), (be[:, 0], be[:, 1])), shape=(len(bv), len(bv)))
                bc, bl = connected_components(bg, directed=False)
                for ri in range(bc):
                    rv = bv[bl == ri]
                    points = v[rv]
                    center = points.mean(axis=0)
                    singular = np.linalg.svd(points-center, full_matrices=False)
                    radius = np.linalg.norm(points-center, axis=1)
                    edges_in_ring = b[np.isin(b[:, 0], rv)]
                    degree = np.bincount(edges_in_ring.ravel(), minlength=len(v))[rv]
                    array_key = f'p{pi}_c{ci}_r{ri}_literal_points'
                    arrays[array_key] = points
                    rings.append(dict(array=array_key, point_count=len(points),
                                      all_boundary_degrees_two=bool(np.all(degree == 2)),
                                      center_m=center.tolist(), radius_min_max_m=[float(radius.min()), float(radius.max())],
                                      plane_residual_max_m=float(np.abs((points-center)@singular[2][-1]).max())))
            arrays[f'p{pi}_c{ci}_source_face_ids'] = face_ids
            components.append(dict(index=ci, faces=len(face_ids), bounds_m=[v[used].min(0).tolist(), v[used].max(0).tolist()],
                                   boundary_rings=rings, two_ring_tube_outer_surface_only=len(rings) == 2))
        records.append(dict(id=ident, name=p['name'], active=ident not in hidden and not p.get('suppressed', False),
                            license=p.get('license'), provenance=p.get('provenance'), geometry_sha256=mesh.geometry_sha,
                            bank_sha256=reader.file_sha(reader.asset(atlas['chunks'][p['chunk']]['url'])),
                            faces=len(f), component_count=component_count,
                            nonmanifold_edges=int(np.sum(incidence > 2)),
                            components=components))
    np.savez_compressed(out / 'literal-components.npz', **arrays)
    inventory = [{k: t[k] for k in ('id', 'name_en', 'parent_id', 'representation_status', 'active_exact_source_part_ids',
                                  'resolved_descendant_only_part_ids')} for t in terms]
    evidence = dict(schema='pulmonary-circulation-literal-source-frontier/v1', source_catalogue=str(source),
                    source_catalogue_sha256=source_sha, census_sha256=digest_file(census_path),
                    literal_arrays_sha256=digest_file(out / 'literal-components.npz'),
                    family_inventory=inventory, sources=records,
                    geometry_changed=False, catalogue_changed=False, anatomy_accepted=False,
                    limitations=['Tube components are outer surfaces, without finite walls or proven connected lumina.',
                                 'No component is assigned an intra/intersegmental identity by proximity alone.',
                                 'New semantic controls need source course and bronchovascular segment reference review.',
                                 'No native component is silently moved to a different recipient or registered by bounds.'],
                    elapsed_seconds=time.monotonic()-began)
    (out / 'evidence.json').write_text(json.dumps(evidence, indent=2)+'\n')
    print(json.dumps(dict(output=str(out), source_count=len(records), component_count=sum(r['component_count'] for r in records),
                          source_catalogue_sha256=source_sha, evidence_sha256=digest_file(out / 'evidence.json'),
                          elapsed_seconds=evidence['elapsed_seconds'])))


if __name__ == '__main__':
    main()
