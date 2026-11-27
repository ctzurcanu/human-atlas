"""Bounded C++ clearance proposals for rigid profile centres.

Occupancy assumes closed, disjoint recipient envelopes. Native proposals are
never a clearance certificate; emitted indexed faces need independent checks.
"""
import time
import numpy as np
from geometry_tools.open3d_backend import library


def project_profile_centres(index, points, support_m, movable, *, gap_m=.00025,
                            maximum_steps=8, maximum_step_m=.003, wall_seconds=2.):
    start = time.monotonic()
    original = np.asarray(points, float)
    q = original.copy(); support = np.asarray(support_m, float); movable = np.asarray(movable, bool)
    if support.shape != (len(q),) or movable.shape != (len(q),):
        raise ValueError('Explicit profile supports and protected/moving mask required')
    steps = []
    for _ in range(maximum_steps):
        closest, _ = index.closest_points(q)
        local = ((q-index.origin)/index.scale).astype(np.float32)
        inside = index.scene.compute_occupancy(library().core.Tensor(local), nthreads=2, nsamples=3).numpy() > .5
        distance = closest['distance_m']
        signed = distance*np.where(inside, -1., 1.)
        remaining = support+gap_m-signed
        active = movable & (remaining > 1e-8)
        steps.append(dict(activeProfiles=int(active.sum()), occupiedProfiles=int(inside.sum()),
                          maximumSampledDeficit_m=float(np.maximum(remaining, 0).max()),
                          protectedViolations=int((~movable & (remaining > 1e-8)).sum())))
        if not active.any():
            break
        direction = q-closest['points']
        direction[inside] *= -1
        length = np.linalg.norm(direction, axis=1)
        near = length < 1e-12
        if near.any():
            t = index.mesh.triangles[closest['face_ids'][near]]
            direction[near] = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
            length[near] = np.linalg.norm(direction[near], axis=1)
        if np.any(length[active] == 0):
            raise ValueError('Recipient has no nonzero proposal direction')
        direction = np.divide(direction, length[:, None], out=np.zeros_like(direction), where=length[:, None] > 0)
        q[active] += direction[active]*np.minimum(remaining[active], maximum_step_m)[:, None]
        q[~movable] = original[~movable]
        if time.monotonic()-start > wall_seconds:
            raise TimeoutError('Native profile-clearance projection deadline')
    return q, dict(steps=steps, maximumSteps=maximum_steps, maximumStep_m=maximum_step_m,
                   protectedCentresLiteral=bool(np.array_equal(q[~movable], original[~movable])),
                   maximumCentreCorrection_m=float(np.linalg.norm(q-original, axis=1).max()),
                   nativeProposalOnly=True, wholeSurfaceClearanceCertified=False,
                   elapsedSeconds=time.monotonic()-start)
