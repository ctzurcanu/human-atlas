"""Open3D acceleration and rigid registration proposals for file-backed meshes.

Raycasting uses normalized Float32 copies. Returned nearest-face witnesses are
recomputed on the literal Float64 source triangles; ``certify=True`` additionally
searches the existing Float64 BVH. Fast results alone are not acceptance gates.
No smoothing, welding by tolerance, reconstruction or atlas installation occurs.
"""
from functools import lru_cache
from importlib import metadata, util
from pathlib import Path
import time

import numpy as np

from .common import atomic_npz, digest_file
from .distance import point_triangle, triangle_pairs


def available():
    return util.find_spec('open3d') is not None


@lru_cache(maxsize=1)
def library():
    try:
        import open3d
    except ImportError as exc:
        raise RuntimeError('Install geometry_tools/requirements.txt in the geometry Python environment to use Open3D') from exc
    return open3d


def status():
    out = dict(available=available(), device='CPU:0', lazy_import=True)
    if out['available']:
        out['version'] = metadata.version('open3d')
    return out


def _points(value, name='points'):
    q = np.asarray(value, dtype=np.float64)
    if q.ndim != 2 or q.shape[1] != 3 or not np.isfinite(q).all():
        raise ValueError(f'{name} must have finite shape [N,3] in metres')
    return q


class SurfaceIndex:
    """Reusable native triangle index; primitive IDs retain source face order."""
    def __init__(self, mesh, *, nthreads=2):
        if not isinstance(nthreads, int) or not 1 <= nthreads <= 64:
            raise ValueError('Open3D thread count must be between 1 and 64')
        if not len(mesh.faces):
            raise ValueError('Open3D surface queries require nonempty triangles')
        self.mesh, self.nthreads = mesh, nthreads
        used, inv = np.unique(mesh.faces, return_inverse=True)
        v = mesh.vertices[used]
        self.origin = (v.min(0) + v.max(0)) * .5
        self.scale = float(np.ptp(v, axis=0).max())
        if self.scale <= 0:
            raise ValueError('A point-only mesh cannot form a raycasting scene')
        local = ((v - self.origin) / self.scale).astype(np.float32)
        f = inv.reshape(-1, 3).astype(np.uint32)
        t = local[f].astype(np.float64)
        if np.any(np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) == 0):
            raise ValueError('Open3D Float32 normalization collapses a source triangle; use the Float64 BVH')
        self.vertex_rounding_m = float(np.linalg.norm(local.astype(float) * self.scale + self.origin - v, axis=1).max())
        o3d = library()
        self.scene = o3d.t.geometry.RaycastingScene(nthreads=nthreads)
        self.scene.add_triangles(o3d.core.Tensor(local), o3d.core.Tensor(f))

    def evidence(self):
        return dict(backend='open3d', version=library().__version__, device='CPU:0',
                    geometry_sha256=self.mesh.geometry_sha, frame=self.mesh.frame,
                    local_origin_m=self.origin.tolist(), normalization_scale_m=self.scale,
                    maximum_vertex_rounding_m=self.vertex_rounding_m,
                    source_triangle_order_preserved=True, anatomical_acceptance=False)

    def closest_points(self, points, *, certify=False, bvh=None):
        q = _points(points)
        if len(q):
            local = ((q - self.origin) / self.scale).astype(np.float32)
            if not np.isfinite(local).all():
                raise ValueError('Query coordinates exceed the normalized Float32 range')
            raw = self.scene.compute_closest_points(library().core.Tensor(local), nthreads=self.nthreads)
            ids = raw['primitive_ids'].numpy().astype(np.int64)
            if np.any(ids < 0) or np.any(ids >= len(self.mesh.faces)):
                raise RuntimeError('Open3D did not return a valid source triangle')
            d2, p, weights = point_triangle(q, self.mesh.triangles[ids])
        else:
            ids = np.empty(0, dtype=np.int64)
            d2, p, weights = np.empty(0), np.empty((0, 3)), np.empty((0, 3))
        arrays = dict(distance_m=np.sqrt(d2), points=p, face_ids=ids, barycentric=weights)
        if certify:
            from .spatial import BVH
            bvh = bvh or BVH(self.mesh)
            if bvh.mesh.geometry_sha != self.mesh.geometry_sha:
                raise ValueError('Certification BVH belongs to different source triangles')
            rows = [bvh.closest_point(point, seed_face=int(face)) for point, face in zip(q, ids)]
            arrays = dict(distance_m=np.array([r['distance_m'] for r in rows]),
                          points=np.array([r['point'] for r in rows]).reshape(-1, 3),
                          face_ids=np.array([r['face_id'] for r in rows], dtype=np.int64),
                          barycentric=np.array([r['barycentric'] for r in rows]).reshape(-1, 3),
                          tested_triangles=np.array([r['tested_triangles'] for r in rows], dtype=np.int64))
        return arrays, {**self.evidence(), 'queries': len(q),
                        'nearest_triangle_verified': bool(certify),
                        'float64_source_witnesses': True,
                        'distance_is_upper_bound_only': not certify,
                        'precision': 'Float32 proposals; Float64 original-triangle witnesses' + (' and BVH search' if certify else '')}

    def cast_rays(self, origins, directions):
        origins, directions = _points(origins, 'origins'), _points(directions, 'directions')
        if origins.shape != directions.shape:
            raise ValueError('Each ray origin needs one direction')
        lengths = np.linalg.norm(directions, axis=1)
        if np.any(lengths == 0):
            raise ValueError('Ray directions must be nonzero')
        directions = directions / lengths[:, None]
        rays = np.c_[(origins - self.origin) / self.scale, directions].astype(np.float32)
        if not np.isfinite(rays).all():
            raise ValueError('Rays exceed the normalized Float32 range')
        if len(rays):
            raw = self.scene.cast_rays(library().core.Tensor(rays), nthreads=self.nthreads)
            hit = np.isfinite(raw['t_hit'].numpy())
            ids = raw['primitive_ids'].numpy().astype(np.int64)
            uv = raw['primitive_uvs'].numpy().astype(float)
        else:
            hit, ids, uv = np.zeros(0, bool), np.empty(0, dtype=np.int64), np.empty((0, 2))
        ids[~hit] = -1
        weights = np.c_[1 - uv.sum(1), uv]
        points = np.full(origins.shape, np.nan)
        points[hit] = np.einsum('ni,nij->nj', weights[hit], self.mesh.triangles[ids[hit]])
        travel = np.full(len(origins), np.nan)
        travel[hit] = np.einsum('ij,ij->i', points[hit] - origins[hit], directions[hit])
        residual = np.full(len(origins), np.nan)
        residual[hit] = np.linalg.norm(points[hit] - origins[hit] - travel[hit, None] * directions[hit], axis=1)
        return dict(hit=hit, face_ids=ids, points=points, barycentric=weights,
                    distance_m=travel, ray_residual_m=residual), {
                        **self.evidence(), 'queries': len(origins), 'hits': int(hit.sum()),
                        'maximum_ray_residual_m': float(residual[hit].max()) if hit.any() else None,
                        'nearest_hit_verified': False, 'approximate_query': True,
                        'miss_representation': 'hit=False, face_ids=-1, points/distance NaN in NPZ'}


def diagnostics(mesh, *, self_intersections=True, include_pairs=False, maximum_export_pairs=100000):
    """Independent C++ diagnostics on exact-coordinate-welded copies only."""
    if include_pairs and (not self_intersections or not isinstance(maximum_export_pairs, int) or maximum_export_pairs < 1):
        raise ValueError('Full source pairs require self-intersection testing and a positive export budget')
    used = np.unique(mesh.faces)
    v, inverse = np.unique(mesh.vertices[used], axis=0, return_inverse=True)
    mapping = np.full(len(mesh.vertices), -1, dtype=np.int64)
    mapping[used] = inverse
    o3d = library()
    native = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(v),
                                     o3d.utility.Vector3iVector(mapping[mesh.faces].astype(np.int32)))
    crossings = np.asarray(native.get_self_intersecting_triangles()) if self_intersections else np.empty((0, 2), int)
    confirmed = []; rejected = 0; examples = []
    for start in range(0, len(crossings), 2048):
        pairs = crossings[start:start + 2048]
        actual = triangle_pairs(mesh.triangles[pairs[:, 0]], mesh.triangles[pairs[:, 1]])
        defects = np.isin(actual['relation'], ['crossing', 'coplanar-overlap'])
        confirmed.extend(pairs[defects].tolist()); rejected += int((~defects).sum())
        for pair, relation, distance in zip(pairs, actual['relation'], actual['distance_m']):
            if len(examples) < 20:
                examples.append(dict(source_faces=pair.tolist(), strict_relation=str(relation), distance_m=float(distance)))
    edges = np.asarray(native.get_non_manifold_edges(allow_boundary_edges=True))
    vertices = np.asarray(native.get_non_manifold_vertices())
    if include_pairs and len(confirmed) > maximum_export_pairs:
        raise ValueError('Complete source defect list exceeds the declared export budget')
    result = dict(backend='open3d', version=o3d.__version__, geometry_sha256=mesh.geometry_sha,
                frame=mesh.frame, edge_manifold=bool(native.is_edge_manifold(allow_boundary_edges=True)),
                vertex_manifold=bool(native.is_vertex_manifold()), orientable=bool(native.is_orientable()),
                self_intersections_checked=self_intersections, self_intersection_pairs=len(confirmed),
                source_face_pair_examples=confirmed[:20], native_candidate_pairs=len(crossings),
                native_candidates_rejected_by_source_predicates=rejected, native_pair_replay_examples=examples,
                candidate_confirmation='Original Float64 triangles and existing strict crossing/coplanar predicates',
                nonmanifold_edges=len(edges), nonmanifold_vertices=len(vertices),
                exact_coordinate_welding_only=True, source_modified=False,
                acceptance_gate=False, scope='Independent C++ diagnostic; retain exact source predicates and anatomical review')
    if include_pairs:
        result.update(source_face_pairs=confirmed, source_face_pairs_complete=True)
    return result


def _samples(mesh, count, seed):
    t = mesh.triangles
    normals = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    areas = np.linalg.norm(normals, axis=1)
    if not areas.sum() > 0:
        raise ValueError('ICP requires a positive-area triangle surface')
    rng = np.random.default_rng(seed)
    ids = rng.choice(len(t), size=count, p=areas / areas.sum())
    uv = rng.random((count, 2)); root = np.sqrt(uv[:, 0])
    weights = np.c_[1 - root, root * (1 - uv[:, 1]), root * uv[:, 1]]
    return np.einsum('ni,nij->nj', weights, t[ids]), normals[ids] / areas[ids, None]


def icp_transform(source, target, initial, *, target_frame, max_correspondence_m,
                  voxel_sizes_m, max_correction_m, max_rotation_degrees,
                  sample_count=12000, seed=41, method='point-to-point',
                  iterations=40, minimum_fitness=.7):
    """Bounded ICP refinement of an explicit initialization, never installation."""
    from .fit import transform_mesh
    if target.frame != target_frame:
        raise ValueError('ICP target frame must be explicitly declared and match its mesh')
    initial = np.asarray(initial, float)
    initialized = transform_mesh(source, initial, target_frame)
    levels = np.asarray(voxel_sizes_m, float)
    if levels.ndim != 1 or not len(levels) or not np.isfinite(levels).all() or np.any(levels <= 0) or np.any(np.diff(levels) >= 0):
        raise ValueError('ICP voxel sizes must be positive, finite and strictly decreasing in metres')
    if not np.isfinite([max_correspondence_m, max_correction_m, max_rotation_degrees, minimum_fitness]).all() or max_correspondence_m <= 0 or max_correction_m < 0 or not 0 <= max_rotation_degrees <= 180 or not 0 < minimum_fitness <= 1:
        raise ValueError('Explicit finite ICP distance, correction, rotation and fitness limits are required')
    if not isinstance(sample_count, int) or not 100 <= sample_count <= 200000 or not isinstance(iterations, int) or not 1 <= iterations <= 1000:
        raise ValueError('ICP needs 100–200000 samples and 1–1000 iterations per level')
    if method not in ['point-to-point', 'point-to-plane']:
        raise ValueError('ICP method must be point-to-point or point-to-plane')
    o3d = library(); reg = o3d.pipelines.registration
    a, _ = _samples(initialized, sample_count, seed)
    b, normals = _samples(target, sample_count, seed)
    cloud = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(a))
    host = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(b))
    host.normals = o3d.utility.Vector3dVector(normals)
    estimate = reg.TransformationEstimationPointToPoint(with_scaling=False) if method == 'point-to-point' else reg.TransformationEstimationPointToPlane(reg.HuberLoss(max_correspondence_m / 2))
    delta = np.eye(4); history = []; began = time.perf_counter()
    before = reg.evaluate_registration(cloud, host, max_correspondence_m, delta)
    for voxel in levels:
        src, dst = cloud.voxel_down_sample(float(voxel)), host.voxel_down_sample(float(voxel))
        if min(len(src.points), len(dst.points)) < 6:
            raise ValueError('ICP voxel size leaves fewer than six surface samples')
        result = reg.registration_icp(src, dst, max_correspondence_m, delta, estimate,
                                      reg.ICPConvergenceCriteria(relative_fitness=1e-8, relative_rmse=1e-8, max_iteration=iterations))
        delta = np.asarray(result.transformation)
        history.append(dict(voxel_m=float(voxel), source_samples=len(src.points), target_samples=len(dst.points),
                            fitness=float(result.fitness), inlier_rmse_m=float(result.inlier_rmse)))
    matrix = delta @ initial
    candidate = transform_mesh(source, matrix, target_frame)
    correction = float(np.linalg.norm(candidate.vertices - initialized.vertices, axis=1).max())
    angle = float(np.degrees(np.arccos(np.clip((np.trace(delta[:3, :3]) - 1) / 2, -1, 1))))
    final = reg.evaluate_registration(cloud, host, max_correspondence_m, delta)
    if correction > max_correction_m or angle > max_rotation_degrees:
        raise ValueError(f'ICP correction exceeds declared bounds: {correction:g} m, {angle:g} degrees')
    if final.fitness < minimum_fitness or final.inlier_rmse > before.inlier_rmse and before.fitness >= minimum_fitness:
        raise ValueError('ICP has insufficient overlap or worsens the initialized fit')
    return matrix, dict(backend='open3d', version=o3d.__version__, method=method,
                        source_geometry_sha256=source.geometry_sha, target_geometry_sha256=target.geometry_sha,
                        target_frame=target_frame, sample_count=sample_count, sample_seed=seed,
                        initial_transform=initial.tolist(), initial_fitness=float(before.fitness),
                        initial_inlier_rmse_m=float(before.inlier_rmse), fitness=float(final.fitness),
                        inlier_rmse_m=float(final.inlier_rmse), maximum_correction_m=correction,
                        rotation_correction_degrees=angle, levels=history,
                        elapsed_seconds=time.perf_counter() - began, scale_refined=False,
                        anatomical_correspondence_verified=False)


def run_query(args, db):
    from .cli import load_ref
    _, mesh = load_ref(db, args.mesh, args.model, args.frame)
    output = Path(args.output).resolve()
    if output.exists():
        raise FileExistsError(output)
    with np.load(args.queries, allow_pickle=False) as z:
        if str(z['frame']) != mesh.frame:
            raise ValueError('Surface queries and mesh must have the same explicit coordinate frame')
        index = SurfaceIndex(mesh, nthreads=args.threads)
        if args.rays:
            if args.certify:
                raise ValueError('Ray queries are approximate; --certify is available for nearest points')
            arrays, out = index.cast_rays(z['origins'], z['directions'])
        else:
            arrays, out = index.closest_points(z['points'], certify=args.certify)
    atomic_npz(output, **arrays, frame=np.array(mesh.frame), geometry_sha256=np.array(mesh.geometry_sha))
    return {**out, 'schema': 'surface-query/v1', 'output': str(output),
            'output_sha256': digest_file(output), 'queries_sha256': digest_file(args.queries)}


def add_parser(sub):
    p = sub.add_parser('open3d-status', help='Inspect the optional native backend without importing its GUI')
    p.set_defaults(handler=lambda args, db: status())
    p = sub.add_parser('surface-query', help='Batched source-triangle nearest points or ray intersections')
    p.add_argument('mesh'); p.add_argument('queries'); p.add_argument('--output', required=True)
    p.add_argument('--model'); p.add_argument('--frame'); p.add_argument('--rays', action='store_true')
    p.add_argument('--certify', action='store_true'); p.add_argument('--threads', type=int, default=2)
    p.set_defaults(handler=run_query)
