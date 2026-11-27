"""SAT exit proposals filtered by exact rigid-station motion support."""
import numpy as np


def reachable_exit_axes(first, second, first_stations, second_stations, *, maximum_motion_m, gap_m):
    a, b = np.asarray(first, float), np.asarray(second, float)
    sa, sb = np.asarray(first_stations), np.asarray(second_stations)
    if a.shape != b.shape or a.ndim != 3 or a.shape[1:] != (3, 3) or sa.shape != a.shape[:2] or sb.shape != b.shape[:2]:
        raise ValueError('Literal paired triangles and their three station IDs required')
    ea, eb = np.roll(a, -1, axis=1)-a, np.roll(b, -1, axis=1)-b
    na, nb = np.cross(ea[:, 0], ea[:, 1]), np.cross(eb[:, 0], eb[:, 1])
    axes = np.concatenate([na[:, None], nb[:, None], np.cross(ea[:, :, None], eb[:, None]).reshape(-1, 9, 3),
                           np.cross(na[:, None], ea), np.cross(nb[:, None], eb)], axis=1)
    lengths = np.linalg.norm(axes, axis=2)
    axes = np.divide(axes, lengths[:, :, None], out=np.zeros_like(axes), where=lengths[:, :, None] > 1e-25)
    valid_axes = lengths > 1e-25
    axes = np.concatenate([axes, -axes], axis=1)
    valid_axes = np.concatenate([valid_axes, valid_axes], axis=1)
    delta = a[:, :, None]-b[:, None, :]
    needed = gap_m-np.einsum('pvwc,pac->pavw', delta, axes)
    support = maximum_motion_m*((sa[:, :, None] >= 0).astype(int)+(sb[:, None, :] >= 0).astype(int))
    support = np.where(sa[:, :, None] == sb[:, None, :], 0., support)
    reachable = valid_axes & np.all(needed <= support[:, None]+1e-12, axis=(2, 3))
    relative = np.divide(np.maximum(needed, 0), support[:, None], out=np.zeros_like(needed), where=support[:, None] > 0)
    cost = relative.max(axis=(2, 3))
    cost[~reachable] = np.inf
    selected = cost.argmin(1)
    index = np.arange(len(a))
    available = np.isfinite(cost[index, selected])
    return axes[index, selected], available, dict(pairCount=len(a),
           pairsWithNecessaryMotionSupport=int(available.sum()), pairsUnreachableWithinDeclaredStations=int((~available).sum()),
           perPairSupportFraction=np.where(available, cost[index, selected], 0).tolist(),
           scope='Necessary independent corner-motion support only; coupled feasibility and actual triangles still require verification')
