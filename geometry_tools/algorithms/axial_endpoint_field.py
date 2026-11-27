"""Bounded endpoint translation with literal fixed controls and normal transport.

Every plane perpendicular to ``axis`` receives one uniform translation: its
profile is preserved exactly. This does not certify calibre in oblique clinical
sections, mouth correspondence, or anatomy. The analytic positive-Jacobian
certificate is independent of sampled mesh quality checks.
"""
import numpy as np


def endpoint_translation(vertices, normals, axis, start, stop, delta, *, minimum_determinant=.7):
    v = np.asarray(vertices, float)
    n = np.asarray(normals, float)
    axis = np.asarray(axis, float)
    delta = np.asarray(delta, float)
    if v.ndim != 2 or v.shape[1] != 3 or n.shape != v.shape:
        raise ValueError('Corresponding vertex/normal vectors required')
    if axis.shape != (3,) or delta.shape != (3,) or not np.isfinite(np.r_[v.ravel(), n.ravel(), axis, delta, start, stop]).all():
        raise ValueError('Finite positions, normals and field parameters required')
    if not stop > start or np.linalg.norm(axis) == 0 or minimum_determinant <= 0:
        raise ValueError('Ordered endpoint planes and positive determinant bound required')
    axis = axis / np.linalg.norm(axis)
    # Quintic smoothstep has maximum derivative 15/8 at t=1/2.
    derivative_max = 1.875 / (stop - start)
    dot = float(axis @ delta)
    determinant_min = min(1., 1. + derivative_max * dot)
    if determinant_min < minimum_determinant:
        raise ValueError('The complete continuous endpoint field violates its determinant bound')
    axial = v @ axis
    t = np.clip((axial - start) / (stop - start), 0., 1.)
    weight = t**3 * (10. - 15. * t + 6. * t*t)
    derivative = 30. * t*t * (1. - t)**2 / (stop - start)
    moved = v + weight[:, None] * delta
    # Sherman-Morrison inverse transpose of I + delta (axis*w')^T.
    transported = n - axis * (derivative * (n @ delta) / (1. + derivative * dot))[:, None]
    length = np.linalg.norm(transported, axis=1)
    if np.any(length == 0):
        raise ValueError('Cannot transport zero source normals')
    transported /= length[:, None]
    fixed = axial <= start
    moved[fixed] = v[fixed]
    transported[fixed] = n[fixed]
    return moved, transported, dict(
        continuous_determinant_lower_bound=determinant_min,
        continuous_jacobian_deviation_upper_bound=float(derivative_max * np.linalg.norm(delta)),
        maximum_displacement_metres=float(np.linalg.norm(delta)),
        fixed_control_count=int(fixed.sum()), fixed_controls_literal=bool(np.array_equal(moved[fixed], v[fixed])),
        parallel_section_profiles_preserved=True,
        oblique_section_calibre_verified=False,
    )
