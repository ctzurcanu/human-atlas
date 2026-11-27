"""Shared, smooth digit-course proposal; apply its poses to rigid tube rings.

Named source/target bone guides determine correspondence, not anatomical
acceptance. The carrier does not shrink a vessel radius. Recipient surface
clearance, discrete element validity and conduit interfaces need later checks.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree
from .tube_sections import arclength, rotation_minimizing_frames, unit


def proper_polar(matrices):
    u, _, vh = np.linalg.svd(matrices)
    sign = np.ones(matrices.shape[:-2]+(3,))
    sign[..., -1] = np.linalg.det(u@vh)
    return (u*sign[..., None, :])@vh


class DigitCourseCarrier:
    def __init__(self, courses, *, blend_width_m=.0025, proximal_fade_fraction=.3,
                 samples=256, initial_reference=(0., 0., 1.), radial_surfaces=None):
        if not courses or blend_width_m <= 0 or not 0 < proximal_fade_fraction < 1:
            raise ValueError('Explicit courses, positive blend width and bounded fade required')
        if not 32 <= samples <= 2048:
            raise ValueError('Bounded course sample count required')
        self.width = blend_width_m
        self.fade = proximal_fade_fraction
        self.radial_surfaces = radial_surfaces
        if radial_surfaces is not None and len(radial_surfaces) != len(courses):
            raise ValueError('One source/recipient surface pair per declared digit required')
        self.courses = []
        for entry in courses:
            a, b = np.asarray(entry['source'], float), np.asarray(entry['target'], float)
            if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3:
                raise ValueError('Corresponding finite 3D guide points required')
            s = arclength(a); s /= s[-1]
            old, new = PchipInterpolator(s, a), PchipInterpolator(s, b)
            t = np.linspace(0., 1., samples)
            x, y = old(t), new(t)
            f = rotation_minimizing_frames(x, unit(old.derivative()(t)), initial_reference)
            g = rotation_minimizing_frames(y, unit(new.derivative()(t)), initial_reference)
            rotations = g@f.transpose(0, 2, 1)
            self.courses.append((entry.get('digit'), old, new, t, x, rotations, cKDTree(x)))

    def poses(self, points, profile_support_m=None):
        p = np.asarray(points, float)
        if p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(p).all():
            raise ValueError('Finite profile centres required')
        support = np.zeros(len(p)) if profile_support_m is None else np.asarray(profile_support_m, float)
        if support.shape != (len(p),) or not np.isfinite(support).all() or np.any(support < 0):
            raise ValueError('One finite nonnegative actual profile support per centre required')
        distance, deltas, rotations, fades, radial = [], [], [], [], []
        for ci, (_, old, new, t, x, R, tree) in enumerate(self.courses):
            _, nearest = tree.query(p)
            starts = np.stack([np.maximum(nearest-1, 0), np.minimum(nearest, len(t)-2)], axis=1)
            a, b = x[starts], x[starts+1]
            e = b-a
            fraction = np.clip(np.einsum('nij,nij->ni', p[:, None]-a, e)/np.sum(e*e, axis=2), 0., 1.)
            proposals = a+fraction[..., None]*e
            choice = np.argmin(np.sum((proposals-p[:, None])**2, axis=2), axis=1)
            k = starts[np.arange(len(p)), choice]
            z = fraction[np.arange(len(p)), choice]
            parameter = t[k]+z*(t[k+1]-t[k])
            centre = old(parameter)
            pose = proper_polar(R[k]*(1-z[:, None, None])+R[k+1]*z[:, None, None])
            target_axis = new(parameter)
            offset = p-centre
            posed_offset = np.einsum('nij,nj->ni', pose, offset)
            posed = target_axis+posed_offset
            radius = np.linalg.norm(offset, axis=1)
            radial_info = dict(applied=0, sourceAndRecipientHits=0, maximumCorrection_m=0.)
            if self.radial_surfaces is not None:
                src, dst = self.radial_surfaces[ci]
                valid = (radius > 1e-7) & (radius < .025)
                ids = np.flatnonzero(valid)
                direction = offset[ids]/radius[ids, None]
                new_direction = posed_offset[ids]/radius[ids, None]
                first, _ = src.cast_rays(centre[ids], direction)
                second, _ = dst.cast_rays(target_axis[ids], new_direction)
                hits = first['hit'] & second['hit']
                old_radius, new_radius = first['distance_m'], second['distance_m']
                # A near-zero hit can be an internal partition cap at a guide
                # point. Do not turn it into an asserted outer cortex radius.
                usable = hits & (old_radius > 5e-5) & (new_radius > 5e-5)
                usable &= (old_radius < radius[ids]+.004) & (new_radius < .015)
                desired_gap = np.maximum(radius[ids]-old_radius, support[ids]+.00025)
                correction = new_radius+desired_gap-radius[ids]
                usable &= np.abs(correction) <= .008
                active = ids[usable]
                posed[active] += new_direction[usable]*correction[usable, None]
                radial_info = dict(applied=int(usable.sum()), sourceAndRecipientHits=int(hits.sum()),
                                   maximumCorrection_m=float(np.max(abs(correction[usable]))) if usable.any() else 0.,
                                   rayQueriesAreProposalsOnly=True, minimumWorkingGapBeyondProfile_m=.00025)
            radial.append(radial_info)
            u = np.clip(parameter/self.fade, 0., 1.)
            fade = u*u*u*(10.+u*(-15.+6.*u))
            distance.append(np.sum((centre-p)**2, axis=1))
            deltas.append(posed-p); rotations.append(pose); fades.append(fade)
        squared = np.array(distance).T
        weights = np.exp(-(squared-squared.min(axis=1)[:, None])/(2*self.width**2))
        weights /= weights.sum(axis=1)[:, None]
        effective = weights*np.array(fades).T
        target = p+np.einsum('ni,nij->nj', effective, np.stack(deltas, axis=1))
        matrix = np.eye(3)+np.einsum('ni,nijk->njk', effective,
                                     np.stack(rotations, axis=1)-np.eye(3))
        pose = proper_polar(matrix)
        return target, pose, dict(course_weights=weights, effective_weights=effective,
                                  source_distance_m=np.sqrt(squared), radial_surface_proposals=radial)
