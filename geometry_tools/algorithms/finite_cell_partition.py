"""Partition a finite cell without inventing another overlapping coat.

The supplied closed classifier only selects a region. Its internal cut is a
single canonical surface used with opposite signs; original surface fragments
retain per-triangle ancestry for adjoining-cell and attribute reconstruction.
No anatomical identity, classifier dimensions or tissue properties are inferred.
"""
from collections import defaultdict
import numpy as np

from .finite_tissue_boolean import (
    CommonPositionBank, FiniteBooleanError, barycentric, boundary_signature,
    make_solid, provenance,
)
from geometry_tools.distance import triangle_pairs


def partition_closed_cell(positions, faces, classifier_positions, classifier_faces,
                          *, intersection_tolerance_m=2e-12):
    positions, faces = np.asarray(positions, float), np.asarray(faces, int)
    host, host_id = make_solid(positions, faces)
    classifier, _ = make_solid(classifier_positions, classifier_faces)
    # The two complementary operations retain the same input face identities.
    solids = {'inside': host ^ classifier, 'outside': host - classifier}
    bank = CommonPositionBank(positions, intersection_tolerance_m)
    source, cut, volumes = {}, {}, {}
    for name, solid in solids.items():
        if str(solid.status()) != 'Error.NoError' or solid.is_empty():
            raise FiniteBooleanError(f'{name}: classifier must form two finite nonempty cells')
        mesh = solid.to_mesh64()
        origins, parents = provenance(mesh)
        ids = np.array([bank.insert(v) for v in np.asarray(mesh.vert_properties)[:, :3]], int)
        actual = ids[np.asarray(mesh.tri_verts)]
        take = origins == host_id
        source[name] = (actual[take], parents[take])
        cut[name] = actual[~take]
        volumes[name] = float(solid.volume())
    # The implementation may choose different diagonals, but both sides must
    # have the same directed cut boundary before one triangulation is reused.
    canonical = cut['inside']
    if not len(canonical) or boundary_signature(canonical) != boundary_signature(cut['outside'][:, [0, 2, 1]]):
        raise FiniteBooleanError('Complementary classifier cuts have incompatible boundaries')
    p = bank.positions
    by_parent = defaultdict(lambda: defaultdict(list))
    for name, (triangles, parents) in source.items():
        for triangle, parent in zip(triangles, parents):
            if not 0 <= parent < len(faces):
                raise FiniteBooleanError('Unknown original cell face')
            by_parent[int(parent)][name].append(triangle)
    maximum_area_residual = 0.
    minimum_parent_weight = 1.
    maximum_plane_distance = 0.
    literal_source_faces = 0
    for parent, original in enumerate(positions[faces]):
        fragments = by_parent[parent]
        all_faces = np.concatenate([np.asarray(fs, int) for fs in fragments.values()]) if fragments else np.empty((0, 3), int)
        if not len(all_faces):
            raise FiniteBooleanError(f'Original cell face {parent} disappeared')
        actual = p[all_faces]
        # Literal cyclic equality proves plane, area, winding and complete
        # ownership without solving a nearly singular barycentric system.
        # The check deliberately rejects opposite winding and extra fragments.
        if len(actual)==1 and any(np.array_equal(actual[0],np.roll(original,k,axis=0))for k in range(3)):
            literal_source_faces+=1
            continue
        weights = barycentric(actual.reshape(-1, 3), original)
        replay = weights @ original
        plane_distance = np.linalg.norm(replay-actual.reshape(-1, 3), axis=1).max()
        minimum_parent_weight = min(minimum_parent_weight, float(weights.min()))
        maximum_plane_distance = max(maximum_plane_distance, float(plane_distance))
        altitude_scale = max(np.linalg.norm(original-np.roll(original, -1, axis=0), axis=1).max(), np.finfo(float).tiny)
        # Coordinate error, not an anatomical clearance or free area tolerance.
        if plane_distance > intersection_tolerance_m*4 or weights.min() < -intersection_tolerance_m*4/altitude_scale:
            raise FiniteBooleanError(f'Fragment leaves original triangle {parent}')
        cross = np.cross(original[1]-original[0], original[2]-original[0])
        actual_cross = np.cross(actual[:, 1]-actual[:, 0], actual[:, 2]-actual[:, 0])
        if np.any(actual_cross @ cross <= 0):
            raise FiniteBooleanError(f'Reversed or degenerate fragment of source face {parent}')
        original_area = np.linalg.norm(cross)/2
        actual_area = np.linalg.norm(actual_cross, axis=1).sum()/2
        perimeter = np.linalg.norm(original-np.roll(original, -1, axis=0), axis=1).sum()
        residual = abs(actual_area-original_area)
        maximum_area_residual = max(maximum_area_residual, float(residual))
        if residual > intersection_tolerance_m*perimeter*4 + 1e-20:
            raise FiniteBooleanError(f'Source face {parent} is not completely partitioned')
        if len(fragments) == 2:
            a, b = [p[np.asarray(fragments[k], int)] for k in ['inside', 'outside']]
            first = np.repeat(a, len(b), axis=0)
            second = np.tile(b, (len(a), 1, 1))
            relations = triangle_pairs(first, second)['relation']
            if np.any(np.isin(relations, ['crossing', 'coplanar-overlap'])):
                raise FiniteBooleanError(f'Complementary fragments overlap on source face {parent}')
    cell_faces = {
        'inside': np.r_[source['inside'][0], canonical],
        'outside': np.r_[source['outside'][0], canonical[:, [0, 2, 1]]],
    }
    roundoff = 2048*np.finfo(float).eps*max(np.ptp(positions, axis=0).max(), np.finfo(float).tiny)**3*len(faces)
    canonical_volumes = {}
    for name, triangles in cell_faces.items():
        solid, _ = make_solid(p, triangles)
        canonical_volumes[name] = float(solid.volume())
        if abs(canonical_volumes[name]-volumes[name]) > roundoff:
            raise FiniteBooleanError(f'Canonical interface changed {name} volume')
    if abs(sum(canonical_volumes.values())-host.volume()) > roundoff:
        raise FiniteBooleanError('Partition changed the whole occupied volume')
    return dict(
        positions=p, source_fragments=source, cut_faces=canonical, cell_faces=cell_faces,
        evidence=dict(literal_source_controls=np.array_equal(p[:len(positions)], positions),
                      source_faces_completely_partitioned=len(faces),
                      literal_source_faces_replayed=literal_source_faces,
                      maximum_source_face_area_residual_m2=maximum_area_residual,
                      minimum_source_barycentric_weight=minimum_parent_weight,
                      maximum_source_plane_distance_m=maximum_plane_distance,
                      source_volume_m3=float(host.volume()), cell_volumes_m3=canonical_volumes,
                      numerical_volume_bound_m3=float(roundoff),
                      intersection_tolerance_m=intersection_tolerance_m,
                      maximum_intersection_snap_m=bank.maximum_snap,
                      one_reciprocal_interface=True, anatomy_inferred=False))
