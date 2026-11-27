"""Exact source singular values and local proper-pose derivatives.

Sequential halfspaces propose a correction; they are not global bounds on
nonlinear section poses. All original facets must pass exact replay before a
candidate can be emitted. The original reference and inherited budgets never
change between iterations.
"""
import numpy as np
from scipy import sparse
from scipy.spatial.transform import Rotation
from .rigid_section_twists import replay_section_twists


def skew(v):
    v = np.asarray(v, float)
    out = np.zeros(v.shape[:-1] + (3, 3))
    out[..., 0, 1] = -v[..., 2]; out[..., 0, 2] = v[..., 1]
    out[..., 1, 0] = v[..., 2]; out[..., 1, 2] = -v[..., 0]
    out[..., 2, 0] = -v[..., 1]; out[..., 2, 1] = v[..., 0]
    return out


def right_rotation_jacobian(v):
    v = np.asarray(v, float)
    theta2 = np.einsum('...i,...i->...', v, v)
    small = theta2 < 1e-8
    a = np.empty_like(theta2); b = np.empty_like(theta2)
    a[small] = .5-theta2[small]/24+theta2[small]**2/720
    b[small] = 1/6-theta2[small]/120+theta2[small]**2/5040
    theta = np.sqrt(theta2[~small])
    a[~small] = (1-np.cos(theta))/theta2[~small]
    b[~small] = (theta-np.sin(theta))/(theta2[~small]*theta)
    K = skew(v)
    return np.eye(3)-a[..., None, None]*K+b[..., None, None]*(K@K)


class OriginalSingularPoseRows:
    def __init__(self, reference, current, groups, centers, levers, *, nominal_budget=.075):
        if reference.frame != current.frame or not np.array_equal(reference.faces, current.faces):
            raise ValueError('Literal matching original source faces in one frame required')
        self.current = current
        self.groups = np.asarray(groups, int)
        self.centers = np.asarray(centers, float)
        self.levers = np.asarray(levers, float)
        self.count = len(self.centers)
        if not 0 < nominal_budget < 1 or self.groups.shape != (len(current.vertices),):
            raise ValueError('Explicit source profile ownership and physical budget required')
        E = (reference.triangles[:, 1:]-reference.triangles[:, :1]).transpose(0, 2, 1)
        _, s, v = np.linalg.svd(E, full_matrices=False)
        if np.any(s <= 0):
            raise ValueError('Nondegenerate original source facets required')
        self.chart = v.transpose(0, 2, 1)/s[:, None, :]
        zero = np.zeros((2*self.count, 3))
        self.initial_singular = self.values(zero)
        self.budget = np.maximum(nominal_budget, abs(self.initial_singular-1).max(axis=1))

    def values(self, parameters):
        q, _ = replay_section_twists(self.current.vertices, self.groups,
                                     self.centers, self.levers, parameters)
        t = q[self.current.faces]
        return np.linalg.svd((t[:, 1:]-t[:, :1]).transpose(0, 2, 1)@self.chart,
                             compute_uv=False)

    def violations(self, parameters, *, tolerance=1e-8):
        s = self.values(parameters)
        excess = abs(s-1).max(axis=1)-self.budget
        bad = np.flatnonzero(excess > tolerance)
        return bad, dict(complete=True, allOriginalFacetsChecked=len(s),
            violatingOriginalFacets=len(bad), maximumBudgetExcess=float(excess.max(initial=0)),
            maximumOriginalStrain=float(abs(s-1).max(initial=0)),
            originalReferenceAndBudgetsUnchanged=True, exactProperPoseReplay=True)

    def singular_gradients(self, parameters, face_ids):
        x = np.asarray(parameters, float)
        ids = np.asarray(face_ids, int)
        q, _ = replay_section_twists(self.current.vertices, self.groups,
                                     self.centers, self.levers, x)
        faces = self.current.faces[ids]
        t = q[faces]
        F = (t[:, 1:]-t[:, :1]).transpose(0, 2, 1)@self.chart[ids]
        u, s, vh = np.linalg.svd(F, full_matrices=False)
        phi = x[self.count:]/self.levers[:, None]
        rotations = Rotation.from_rotvec(phi).as_matrix()
        jacobian = right_rotation_jacobian(phi)
        rr, cc, vv = [], [], []
        for singular in range(2):
            G = np.einsum('ni,nj->nij', u[:, :, singular], vh[:, singular, :])
            H = G@self.chart[ids].transpose(0, 2, 1)
            point_gradient = np.stack([-H.sum(axis=2), H[:, :, 0], H[:, :, 1]], axis=1)
            for corner in range(3):
                groups = self.groups[faces[:, corner]]
                live = groups >= 0
                g = groups[live]
                v = self.current.vertices[faces[live, corner]]-self.centers[g]
                J = -rotations[g]@skew(v)@jacobian[g]/self.levers[g, None, None]
                translation = point_gradient[live, corner]
                rotation = np.einsum('ni,nij->nj', translation, J)
                rows = np.repeat(np.flatnonzero(live)*2+singular, 3)
                rr.extend([rows, rows])
                cc.extend([(g[:, None]*3+np.arange(3)).ravel(),
                           (3*self.count+g[:, None]*3+np.arange(3)).ravel()])
                vv.extend([translation.ravel(), rotation.ravel()])
        gradient = sparse.coo_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))),
                                      shape=(2*len(ids), 6*self.count)).tocsr()
        gradient.eliminate_zeros()
        return s, gradient

    def local_rows(self, parameters, face_ids):
        ids = np.asarray(face_ids, int)
        s, gradient = self.singular_gradients(parameters, ids)
        budget = np.repeat(self.budget[ids], 2)
        linear = gradient@np.asarray(parameters).ravel()
        norm = np.sqrt(np.asarray(gradient.power(2).sum(axis=1)).ravel())
        live = norm > 1e-12
        if np.any(~live & (abs(s.ravel()-1) > budget+1e-8)):
            raise ValueError('A violating facet has no resolved pose derivative')
        # Normalize to displacement units, preserving each literal inequality.
        # The halfspace replay permits at most 1e-10 in row units. Keep the
        # implied singular-value residual <=1e-9 even on very small facets.
        # Positive scaling changes conditioning, never the exact inequality.
        scale = np.maximum(1., norm[live]/10.)
        multiplier = scale/norm[live]
        G = sparse.diags(multiplier)@gradient[live]
        upper = (s.ravel()[live]-(1+budget[live])-linear[live])*multiplier
        lower = ((1-budget[live])-s.ravel()[live]+linear[live])*multiplier
        return sparse.vstack([-G, G], format='csr'), np.r_[upper, lower], dict(
            complete=True, sourceFaceIds=ids.tolist(), localLinearizationOnly=True,
            originalReferenceAndBudgetsUnchanged=True, exactReplayRequired=True,
            maximumImpliedSingularResidualAt1e10RowTolerance=1e-9,
            positiveRowScalingDoesNotChangeInequalities=True,
            derivative='analytic singular-value and proper SO(3) right-Jacobian chain rule')
