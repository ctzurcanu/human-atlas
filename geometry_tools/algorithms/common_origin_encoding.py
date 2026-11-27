"""Encode shared controls into compatible, explicitly declared F32 origins.

Different source banks may use different origins. Shared controls must decode
to exactly the same world coordinates in every owning bank; independent
rounding can open an otherwise matched anatomical interface.
"""
import numpy as np


def encode_positions(points, origin_masks, origins):
    p = np.asarray(points, float)
    masks = np.asarray(origin_masks, np.uint64)
    origins = np.asarray(origins, float)
    if (p.ndim != 2 or p.shape[1] != 3 or masks.shape != (len(p),)
            or origins.ndim != 2 or origins.shape[1] != 3
            or not 1 <= len(origins) <= 63
            or not np.isfinite(p).all() or not np.isfinite(origins).all()
            or np.any(masks >> np.uint64(len(origins)))):
        raise ValueError('Finite world controls, valid origin masks and explicit origins required')
    out = p.copy()
    maximum_error = maximum_bound = 0.
    shared = iterations = 0
    for mask in np.unique(masks):
        if not mask:
            continue  # Unused source controls are not emitted.
        ids = np.flatnonzero(masks == mask)
        owns = [i for i in range(len(origins)) if int(mask) & (1 << i)]
        original = p[ids]
        selected_origins = origins[owns]
        scale = np.maximum(abs(original), abs(selected_origins).max(axis=0))
        evaluation_bound = 8 * abs(np.spacing(scale))
        if len(owns) == 1:
            local = original - selected_origins[0]
            if np.any(abs(local) > np.finfo(np.float32).max):
                raise ValueError('Coordinates exceed the declared F32 bank range')
            encoded = local.astype('<f4')
            actual = encoded.astype(float) + selected_origins[0]
            bound = abs(np.spacing(encoded).astype(float)) / 2 + evaluation_bound
        else:
            shared += len(ids)
            anchor = selected_origins[0]
            actual = original.copy()
            step = np.zeros_like(original)
            for iteration in range(4):
                local = actual[None] - selected_origins[:, None]
                if np.any(abs(local) > np.finfo(np.float32).max):
                    raise ValueError('Shared coordinates exceed a declared F32 bank range')
                cast = local.astype('<f4')
                step = np.maximum(step, abs(np.spacing(cast).astype(float)).max(axis=0))
                step = np.maximum(step, abs(np.spacing(original)))
                differences = (selected_origins[:, None] - anchor) / step
                if not np.all(differences == np.rint(differences)):
                    raise ValueError('Declared source origins have incompatible shared F32 lattices')
                actual = np.rint((original - anchor) / step) * step + anchor
                decoded = ((actual[None] - selected_origins[:, None]).astype('<f4').astype(float)
                           + selected_origins[:, None])
                if np.all(decoded == actual[None]):
                    iterations = max(iterations, iteration + 1)
                    break
            else:
                raise ValueError('Shared F32 encoding did not stabilize within four bin checks')
            bound = step / 2 + evaluation_bound
        if not np.isfinite(actual).all() or np.any(abs(actual - original) > bound):
            raise ValueError('Actual coordinate rounding exceeds its derived F32 bound')
        for origin in selected_origins:
            decoded = (actual - origin).astype('<f4').astype(float) + origin
            if not np.array_equal(decoded, actual):
                raise ValueError('An owning bank does not decode the canonical coordinates exactly')
        out[ids] = actual
        maximum_error = max(maximum_error, float(np.linalg.norm(actual - original, axis=1).max(initial=0.)))
        maximum_bound = max(maximum_bound, float(np.linalg.norm(bound, axis=1).max(initial=0.)))
    return out, dict(shared_controls=shared, origin_count=len(origins),
                     maximum_world_rounding_error_m=maximum_error,
                     maximum_derived_rounding_bound_m=maximum_bound,
                     maximum_bin_checks=iterations, all_owning_decoders_exact=True,
                     no_geometric_proximity_welding=True)
