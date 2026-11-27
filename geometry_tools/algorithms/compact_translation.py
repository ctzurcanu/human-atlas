"""Compact C2 translations with a global distortion certificate.

The certificate concerns the continuous spatial map. Transformed triangle meshes
still require their own topology, intersection, tissue and anatomical checks.
"""
from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class CompactTranslation:
    center: np.ndarray
    radii: np.ndarray
    displacement: np.ndarray
    frame: str

    def __post_init__(self):
        for name in ('center', 'radii', 'displacement'):
            value = np.array(getattr(self, name), dtype=float, copy=True)
            if value.shape != (3,) or not np.isfinite(value).all():
                raise ValueError(f'{name} requires three finite metre coordinates')
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if not isinstance(self.frame, str) or not self.frame.strip():
            raise ValueError('An explicit coordinate frame is required')
        if np.any(self.radii <= 0):
            raise ValueError('Support radii must be positive')
        with np.errstate(over='ignore'):
            low, high = self.center-self.radii, self.center+self.radii
        if (not np.isfinite(low).all() or not np.isfinite(high).all()
                or np.any(low >= self.center) or np.any(high <= self.center)):
            raise ValueError('Support bounds cannot be represented at this precision')
        if not self.lipschitz_bound < 1:
            raise ValueError('Displacement exceeds the globally injective bound')

    @property
    def lipschitz_bound(self):
        # max |phi'(r)| = 20*(1/4)*(3/4)^3 = 135/64.
        return (135/64)*math.hypot(*self.displacement)/float(self.radii.min())

    @property
    def minimum_stretch_bound(self):
        return 1-self.lipschitz_bound

    def weights_and_gradients(self, points):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError('Points require a finite [N,3] array')
        weights = np.zeros(len(points))
        gradients = np.zeros(points.shape)
        box = np.all((points >= self.center-self.radii)
                     & (points <= self.center+self.radii), axis=1)
        ids = np.flatnonzero(box)
        offset = points[ids]-self.center
        local = offset/self.radii
        r = np.linalg.norm(local, axis=1)
        live = r < 1
        ids, offset, r = ids[live], offset[live], r[live]
        weights[ids] = (1-r)**4*(4*r+1)
        gradients[ids] = -20*(1-r)[:,None]**3*offset/self.radii**2
        return weights, gradients

    def map(self, points):
        weights, _ = self.weights_and_gradients(points)
        return np.asarray(points, dtype=float)+weights[:,None]*self.displacement

    def jacobians(self, points):
        _, gradient = self.weights_and_gradients(points)
        return np.eye(3)[None,:,:]+self.displacement[None,:,None]*gradient[:,None,:]

    def inverse_transposes(self, points):
        _, gradient = self.weights_and_gradients(points)
        determinant = 1+gradient@self.displacement
        return np.eye(3)[None,:,:]-gradient[:,:,None]*self.displacement[None,None,:]/determinant[:,None,None]

    def transform_normals(self, points, normals):
        normals = np.asarray(normals, dtype=float)
        if normals.shape != np.asarray(points).shape or not np.isfinite(normals).all():
            raise ValueError('Each point requires a finite normal')
        maximum = np.max(np.abs(normals), axis=1)
        if np.any(maximum == 0):
            raise ValueError('Zero normals cannot be transported')
        unit = normals/maximum[:,None]
        unit /= np.linalg.norm(unit, axis=1)[:,None]
        mapped = np.einsum('nij,nj->ni', self.inverse_transposes(points), unit)
        return mapped/np.linalg.norm(mapped, axis=1)[:,None]

    def bound_boxes(self, lower, upper):
        """Enclose mapped boxes and bound the derivative throughout each box.

        The one-dimensional radial derivative peaks at r=1/4. Bounds cover
        complete boxes, rather than just their corners or sampled mesh points.
        Float64 arithmetic includes scale-dependent outward padding.
        """
        lower, upper = np.asarray(lower, float), np.asarray(upper, float)
        if (lower.ndim != 2 or lower.shape[1] != 3 or lower.shape != upper.shape
                or not len(lower) or not np.isfinite(lower).all()
                or not np.isfinite(upper).all() or np.any(lower > upper)):
            raise ValueError('A nonempty set of finite ordered [N,3] boxes is required')
        scale = np.maximum.reduce((np.abs(lower), np.abs(upper),
                                   np.broadcast_to(np.abs(self.center),lower.shape),
                                   np.broadcast_to(self.radii,lower.shape)))
        padding = 128*np.finfo(float).eps*scale
        low, high = lower-padding, upper+padding
        near = np.maximum(np.maximum(low-self.center,self.center-high),0.)/self.radii
        far = np.maximum(abs(low-self.center),abs(high-self.center))/self.radii
        rmin = np.linalg.norm(near,axis=1)
        rmax = np.linalg.norm(far,axis=1)
        a, b = np.clip(rmin,0.,1.), np.clip(rmax,0.,1.)
        maximum = (1-a)**4*(4*a+1)
        minimum = (1-b)**4*(4*b+1)
        peak = np.minimum(np.maximum(.25,a),b)
        gradient_bound = 20*peak*(1-peak)**3/float(self.radii.min())
        gradient_bound *= 1+256*np.finfo(float).eps
        norm_bound = math.hypot(*self.displacement)*gradient_bound
        # J=I+d*g^T has smallest singular value
        # (sqrt((|d||g|)^2+4*(1+d.g))-|d||g|)/2, capped by 1.
        # Bound d.g as well as |g| so shear does not consume the same
        # distortion allowance as an equally large pure contraction.
        products_low = self.displacement*(low-self.center)/self.radii**2
        products_high = self.displacement*(high-self.center)/self.radii**2
        terms = np.maximum(products_low,products_high)
        dot_offset_max = terms.sum(1)+32*np.finfo(float).eps*np.abs(terms).sum(1)
        kmax = 20*(1-a)**3*(1+256*np.finfo(float).eps)
        kmin = 20*(1-b)**3*(1-256*np.finfo(float).eps)
        dot_lower = np.where(dot_offset_max>=0,-kmax*dot_offset_max,-kmin*dot_offset_max)
        dot_lower = np.maximum(dot_lower,-norm_bound)
        dot_lower -= 64*np.finfo(float).eps*norm_bound
        determinant_lower = 1+dot_lower
        lower_stretch = np.zeros(len(lower))
        safe = determinant_lower>0
        lower_stretch[safe] = np.minimum(1.,(np.sqrt(norm_bound[safe]**2+4*determinant_lower[safe])-norm_bound[safe])/2)
        lower_stretch[norm_bound>0] = np.nextafter(lower_stretch[norm_bound>0],-np.inf)
        positive = self.displacement >= 0
        delta_low = self.displacement*np.where(positive,minimum[:,None],maximum[:,None])
        delta_high = self.displacement*np.where(positive,maximum[:,None],minimum[:,None])
        mapped_low = np.nextafter(low+delta_low,-np.inf)
        mapped_high = np.nextafter(high+delta_high,np.inf)
        return mapped_low,mapped_high,lower_stretch

    def evidence(self):
        return dict(center_m=self.center.tolist(), radii_m=self.radii.tolist(),
                    displacement_m=self.displacement.tolist(), frame=self.frame,
                    displacement_lipschitz_bound=self.lipschitz_bound,
                    continuous_minimum_stretch_bound=self.minimum_stretch_bound,
                    continuous_map_globally_injective=True,
                    discrete_mesh_quality_verified=False, anatomical_acceptance=False)


def bound_composition(fields, lower, upper):
    """Carry a complete box cover and its cumulative material stretch bound."""
    lower, upper = np.array(lower, float, copy=True), np.array(upper, float, copy=True)
    stretch = np.ones(len(lower))
    frame = None
    for field in fields:
        if frame is not None and frame != field.frame:
            raise ValueError('Compose fields only in the same explicit frame')
        frame = field.frame
        lower,upper,step_stretch = field.bound_boxes(lower,upper)
        stretch *= step_stretch
    return lower,upper,stretch
