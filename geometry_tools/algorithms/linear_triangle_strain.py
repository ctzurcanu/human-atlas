"""Exact triangle stretch along a prescribed linear displacement field.

The original 3-by-2 triangle bases are factored once.  Subsequent scalar
evaluations use the two eigenvalues of their 2-by-2 Gram matrices; they do not
repeat thousands of small SVDs.  This measures a supplied field, not anatomy.
"""
import numpy as np


class LinearTriangleStrain:
    def __init__(self, vertices, faces, displacement):
        p = np.asarray(vertices, dtype=float)
        f = np.asarray(faces, dtype=int)
        d = np.asarray(displacement, dtype=float)
        if (p.ndim != 2 or p.shape[1] != 3 or d.shape != p.shape
                or f.ndim != 2 or f.shape[1] != 3 or not len(f)
                or f.min() < 0 or f.max() >= len(p)
                or not np.isfinite(p).all() or not np.isfinite(d).all()):
            raise ValueError('Finite source triangles and corresponding field required')
        tri, delta = p[f], d[f]
        source = (tri[:, 1:] - tri[:, :1]).transpose(0, 2, 1)
        u, singular, vh = np.linalg.svd(source, full_matrices=False)
        if np.any(singular[:, 1] <= 0):
            raise ValueError('A degenerate source triangle has no stretch basis')
        changed = (delta[:, 1:] - delta[:, :1]).transpose(0, 2, 1)
        normalized = changed @ (vh.transpose(0, 2, 1) / singular[:, None, :])
        cross = u.transpose(0, 2, 1) @ normalized
        self.linear = cross + cross.transpose(0, 2, 1)
        self.quadratic = normalized.transpose(0, 2, 1) @ normalized

    def stretches(self, weight):
        t = float(weight)
        if not np.isfinite(t):
            raise ValueError('Finite field weight required')
        gram = np.eye(2)[None] + t * self.linear + t*t * self.quadratic
        a, b, c = gram[:, 0, 0], gram[:, 0, 1], gram[:, 1, 1]
        half = (a+c) / 2
        radius = np.hypot((a-c) / 2, b)
        eigenvalues = np.stack((half-radius, half+radius), axis=1)
        if np.any(eigenvalues < -1e-12):
            raise ValueError('Numerically invalid source deformation Gram matrix')
        return np.sqrt(np.maximum(eigenvalues, 0))

    def maximum_strain(self, weight):
        return float(np.abs(self.stretches(weight) - 1).max())
