"""Recover native tube-ring centres and an explicitly non-topological proximity graph."""
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.algorithms.tube_station_topology import recover_station_rings, StationTopologyError
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader

BASE = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began = time.monotonic()
    audit = json.loads((BASE/'evidence.json').read_text())
    contract = json.loads((BASE/'family-contract-v2.json').read_text())
    source = Path(audit['source_catalogue'])
    assert digest_file(source) == contract['source_catalogue_sha256']
    atlas = json.loads(source.read_text())
    parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    arrays = {}
    curves = []
    unsupported = []
    for pi, record in enumerate(audit['sources']):
        part = parts[record['id']]
        mesh = reader.read_part(atlas, part, frame='model:male:assembled-metres')
        assert mesh.geometry_sha == record['geometry_sha256']
        v, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
        f = inverse[mesh.faces]
        try:
            components = recover_station_rings(v, f)
        except StationTopologyError as error:
            unsupported.append(dict(id=part['id'], geometry_sha256=mesh.geometry_sha, reason=str(error),
                                    invented_centerline=False))
            continue
        for ci, component in enumerate(components):
            rings = np.array(component['rings'], dtype=int)
            section_points = v[rings]
            center = section_points.mean(1)
            radius = np.linalg.norm(section_points-center[:,None], axis=2)
            key = f'p{pi}_c{ci}'
            arrays[key+'_centers'] = center
            arrays[key+'_literal_sections'] = section_points
            arrays[key+'_source_face_ids'] = component['face_ids']
            arrays[key+'_radius_min_max'] = np.c_[radius.min(1), radius.max(1)]
            arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(center, axis=0), axis=1))]
            index = len(curves)
            source_flow_hint = 0 if radius[0].mean() >= radius[-1].mean() else 1
            curves.append(dict(index=index, key=key, part_id=part['id'], component=ci,
                               geometry_sha256=mesh.geometry_sha, system=part['system'], stations=len(center),
                               faces=len(component['face_ids']), source_ring_size=rings.shape[1],
                               outer_cast_centerline_not_blood_compartment=True,
                               endpoints_m=[center[0].tolist(), center[-1].tolist()],
                               endpoint_mean_radius_m=[float(radius[0].mean()), float(radius[-1].mean())],
                               larger_radius_endpoint=source_flow_hint, physiological_flow_direction_established=False,
                               larger_to_smaller_radius_net_y_mm=float((center[(-1,0)[source_flow_hint]]-center[(0,-1)[source_flow_hint]])[1]*1000),
                               length_m=float(arc[-1])))
    output = BASE/'native-curves-v3'
    output.mkdir()
    np.savez_compressed(output/'literal-section-centres.npz', **arrays)
    # Endpoint -> other curve is a geometric candidate. It neither joins faces nor
    # establishes branch parenthood, lumen communication, drainage or flow direction.
    graph_rows = []
    for curve in curves:
        center = arrays[curve['key']+'_centers']
        rr = arrays[curve['key']+'_radius_min_max'].mean(1)
        for ei, endpoint in enumerate((center[0], center[-1])):
            nearby = []
            for other in curves:
                if other['index'] == curve['index'] or other['system'] != curve['system']:
                    continue
                path = arrays[other['key']+'_centers']
                radii = arrays[other['key']+'_radius_min_max'].mean(1)
                seg = np.diff(path, axis=0)
                squared = np.einsum('ij,ij->i', seg, seg)
                fraction = np.clip(np.einsum('ij,ij->i', endpoint-path[:-1], seg)/np.maximum(squared, 1e-30), 0., 1.)
                point = path[:-1]+fraction[:,None]*seg
                distance = np.linalg.norm(point-endpoint, axis=1)
                at = int(np.argmin(distance))
                other_radius = float(radii[at]*(1-fraction[at])+radii[at+1]*fraction[at])
                endpoint_radius = float(rr[0 if ei == 0 else -1])
                if distance[at] <= endpoint_radius+other_radius+.001:
                    nearby.append(dict(other_curve=other['index'], other_part_id=other['part_id'],
                                       other_source_component=other['component'],
                                       other_segment=at, other_segment_fraction=float(fraction[at]),
                                       center_distance_m=float(distance[at]),
                                       source_radii_sum_m=endpoint_radius+other_radius,
                                       actual_shared_coordinates_verified=False,
                                       communicating_lumen_verified=False))
            nearby.sort(key=lambda x: (x['center_distance_m'], x['other_curve']))
            graph_rows.append(dict(curve=curve['index'], endpoint=ei, candidates=nearby))
    evidence = dict(schema='native-pulmonary-section-centres-and-proximity/v3',
                    source_catalogue_sha256=contract['source_catalogue_sha256'],
                    family_contract_sha256=digest_file(BASE/'family-contract-v2.json'),
                    source_audit_sha256=digest_file(BASE/'evidence.json'),
                    arrays_sha256=digest_file(output/'literal-section-centres.npz'),
                    curves=curves, unsupported_sources=unsupported,
                    endpoint_proximity_graph=graph_rows,
                    anatomical_assignments=[], geometry_changed=False, catalogue_changed=False,
                    native_sources_count=len(audit['sources']), section_centerline_count=len(curves),
                    graph_is_not_anatomical_or_mesh_topology=True,
                    elapsed_seconds=time.monotonic()-began)
    (output/'evidence.json').write_text(json.dumps(evidence, indent=2)+'\n')
    print(json.dumps(dict(output=str(output), curves=len(curves), unsupported=len(unsupported),
                          endpoint_candidates=sum(len(r['candidates'])for r in graph_rows),
                          evidence_sha256=digest_file(output/'evidence.json'), elapsed_seconds=evidence['elapsed_seconds'])))


if __name__ == '__main__':
    main()
