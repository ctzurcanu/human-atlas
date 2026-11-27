"""Explicit compact maps with analytic positive-Jacobian bounds.

Every point in a dependent assembly must receive the same ordered maps.
Continuous positivity does not establish discrete triangle quality or anatomy.
"""
import numpy as np
from numpy.polynomial import Polynomial
from functools import lru_cache


def compact_radial_map(points, *, center, core_m, support_m, kind,
                       translation_m=None, axis=None, angle_radians=None, growth=None):
    p=np.asarray(points,float); c=np.asarray(center,float)
    core=float(core_m); width=float(support_m)-core
    if p.ndim!=2 or p.shape[1]!=3 or c.shape!=(3,) or core<0 or width<=0:
        raise ValueError('Finite points, center and ordered positive support radii are required')
    if not np.isfinite(np.r_[p.ravel(),c,core,width]).all():
        raise ValueError('Finite compact-map inputs are required')
    x=p-c; r=np.linalg.norm(x,axis=1); t=np.clip((r-core)/width,0,1)
    w=(1-t)**4*(4*t+1)
    derivative=np.where((r>core)&(r<core+width),-20*t*(1-t)**3/width,0.)
    grad=np.divide(x,r[:,None],out=np.zeros_like(x),where=r[:,None]>0)*derivative[:,None]
    eye=np.broadcast_to(np.eye(3),(len(p),3,3)).copy()
    if kind=='translation':
        d=np.asarray(translation_m,float)
        if d.shape!=(3,) or not np.isfinite(d).all(): raise ValueError('Finite translation required')
        lower=1-2.109375*np.linalg.norm(d)/width
        if lower<=0:raise ValueError('Translation support cannot certify a positive Jacobian')
        out=p+w[:,None]*d
        jac=eye+d[None,:,None]*grad[:,None,:]
    elif kind=='twist':
        a=np.asarray(axis,float); angle=float(angle_radians)
        if a.shape!=(3,) or not np.isfinite(np.r_[a,angle]).all() or np.linalg.norm(a)==0:
            raise ValueError('Finite nonzero twist axis and angle required')
        a/=np.linalg.norm(a)
        cross=np.array([[0.,-a[2],a[1]],[a[2],0.,-a[0]],[-a[1],a[0],0.]])
        theta=angle*w; cosine=np.cos(theta); sine=np.sin(theta)
        rot=eye*cosine[:,None,None]+np.outer(a,a)[None]*(1-cosine)[:,None,None]+cross[None]*sine[:,None,None]
        y=np.einsum('nij,nj->ni',rot,x)
        jac=rot+np.cross(a,y)[:,:,None]*(angle*grad)[:,None,:]
        out=c+y
        # Rotation preserves radius; its radial derivative is orthogonal to
        # the radius. The determinant lemma therefore gives det(J) = 1.
        lower=1.
    elif kind=='dilation':
        amp=float(growth)-1
        if not np.isfinite(amp) or 1+amp<=0:raise ValueError('Positive finite radial growth required')
        u=Polynomial([0.,1.]); weight=(1-u)**4*(4*u+1)
        radial=Polynomial([1.])+amp*(weight+(core/width+u)*weight.deriv())
        critical=[0.,1.]+[float(q.real)for q in radial.deriv().roots()if abs(q.imag)<1e-10 and 0<q.real<1]
        radial_lower=min(float(radial(q))for q in critical)
        tangential_lower=min(1.,1+amp)
        lower=radial_lower*tangential_lower**2
        if lower<=0:raise ValueError('Radial growth cannot certify a positive Jacobian')
        scale=1+amp*w
        out=c+x*scale[:,None]
        jac=eye*scale[:,None,None]+amp*x[:,:,None]*grad[:,None,:]
    else:raise ValueError('Unknown compact-map kind')
    # Preserve literal fixed controls outside compact support. Computing
    # c + (p - c) can otherwise change last bits even for the identity map.
    outside=r>=core+width
    out[outside]=p[outside]
    jac[outside]=eye[outside]
    return out,jac,dict(kind=kind,analytic_determinant_lower_bound=float(lower),
                        core_m=core,support_m=core+width,
                        fixed_outside_support=True,discrete_mesh_checks_pending=True)


def compact_ellipsoidal_translation(points, *, center, support_axes_m,
                                     core_fraction, translation_m):
    """Translate a deep corridor while fixing points outside an ellipsoid.

    The diagonal coordinate change preserves the determinant. The determinant
    lemma gives det(J) = 1 + gradient(weight).translation, bounded below by
    1 - 2.109375 * norm(translation / axes) / (1 - core_fraction).
    """
    p = np.asarray(points, float)
    c, axes, d = (np.asarray(value, float) for value in (center, support_axes_m, translation_m))
    core = float(core_fraction)
    if p.ndim != 2 or p.shape[1] != 3 or any(v.shape != (3,) for v in (c, axes, d)):
        raise ValueError('Finite N by 3 coordinates and three-component ellipsoid controls required')
    if not np.isfinite(np.r_[p.ravel(), c, axes, d, core]).all() or np.any(axes <= 0) or not 0 <= core < 1:
        raise ValueError('Positive finite support axes and a core fraction in [0, 1) required')
    width = 1 - core
    lower = 1 - 2.109375 * np.linalg.norm(d / axes) / width
    if lower <= 0:
        raise ValueError('Ellipsoid cannot certify a positive Jacobian for this translation')
    q = (p - c) / axes
    radius = np.linalg.norm(q, axis=1)
    t = np.clip((radius - core) / width, 0, 1)
    weight = (1 - t)**4 * (4*t + 1)
    derivative = np.where((radius > core) & (radius < 1), -20*t*(1-t)**3 / width, 0.)
    gradient = np.divide(q, radius[:, None], out=np.zeros_like(q), where=radius[:, None] > 0)
    gradient *= derivative[:, None] / axes
    out = p + weight[:, None] * d
    jac = np.broadcast_to(np.eye(3), (len(p), 3, 3)).copy() + d[None, :, None] * gradient[:, None, :]
    outside = radius >= 1
    out[outside] = p[outside]
    jac[outside] = np.eye(3)
    return out, jac, dict(kind='ellipsoidal-translation',
        analytic_determinant_lower_bound=float(lower), support_axes_m=axes.tolist(),
        core_fraction=core, fixed_outside_support=True, discrete_mesh_checks_pending=True)


@lru_cache(maxsize=32)
def _corridor_charts(stations):
    from scipy.interpolate import PchipInterpolator
    values = np.asarray(stations, float)
    if values.ndim != 2 or values.shape[1] != 3 or len(values) < 2 or not np.isfinite(values).all() or np.any(np.diff(values[:, 0]) <= 0):
        raise ValueError('Strictly ordered finite Y/X/Z corridor stations required')
    x = PchipInterpolator(values[:, 0], values[:, 1])
    z = PchipInterpolator(values[:, 0], values[:, 2])
    return x, x.derivative(), z, z.derivative()


def separable_corridor_translation(points, *, centerline_stations, center_y_m,
                                   core_axes_m, support_axes_m, translation_m=None,
                                   y_core_interval_m=None, y_support_interval_m=None,
                                   translation_stations=None):
    """Transport complete transverse tube sections in a finite Y corridor.

    Inside its X/Z core each Y section receives one rigid translation, even
    while the translation ramps along Y. The direction must be transverse to
    Y; then centreline slopes do not enter the determinant bound.
    """
    p = np.asarray(points, float)
    core, support = (np.asarray(value, float) for value in (core_axes_m, support_axes_m))
    cy = float(center_y_m)
    if p.ndim != 2 or p.shape[1] != 3 or any(v.shape != (3,) for v in (core, support)):
        raise ValueError('Finite N by 3 points and three-axis corridor controls required')
    if not np.isfinite(np.r_[p.ravel(), core, support, cy]).all() or np.any(core < 0) or np.any(support <= core):
        raise ValueError('Positive transition widths and a strictly transverse translation required')
    if (translation_m is None) == (translation_stations is None):
        raise ValueError('Exactly one constant translation or longitudinal translation course required')
    if translation_stations is None:
        constant = np.asarray(translation_m, float)
        if constant.shape != (3,) or not np.isfinite(constant).all() or constant[1] != 0:
            raise ValueError('Finite strictly transverse translation required')
        maximum = np.abs(constant[[0, 2]])
    else:
        shifts = tuple((float(s['y_m']), float(s['x_m']), float(s['z_m'])) for s in translation_stations)
        tx, tdx, tz, tdz = _corridor_charts(shifts)
        # Zero endpoint slopes let constant continuation stay C1. Reducing
        # the endpoint slopes of PCHIP preserves its interval bounds, which
        # are also the translation maxima used in the determinant bound.
        from scipy.interpolate import CubicHermiteSpline
        def flat_ends(curve):
            slopes=curve.derivative()(curve.x);slopes[[0,-1]]=0.
            return CubicHermiteSpline(curve.x,curve(curve.x),slopes)
        tx,tz=flat_ends(tx),flat_ends(tz);tdx,tdz=tx.derivative(),tz.derivative()
        maximum = np.max(np.abs(np.asarray(shifts)[:, [1, 2]]), axis=0)
    width = support - core
    lower = 1 - 2.109375 * np.sum(maximum / width[[0, 2]])
    if lower <= 0:
        raise ValueError('Corridor cannot certify a positive Jacobian')
    stations = tuple((float(s['y_m']), float(s['x_m']), float(s['z_m'])) for s in centerline_stations)
    px, dx, pz, dz = _corridor_charts(stations)
    y = p[:, 1]
    if translation_stations is None:
        d = np.broadcast_to(constant, p.shape)
        derivative_d = np.zeros_like(p)
    else:
        clipped = np.clip(y, tx.x[0], tx.x[-1])
        d = np.column_stack((tx(clipped), np.zeros(len(p)), tz(clipped)))
        derivative_d = np.column_stack((tdx(clipped), np.zeros(len(p)), tdz(clipped)))
        derivative_d[(y < tx.x[0]) | (y > tx.x[-1])] = 0
    def chart(curve, derivative):
        clipped = np.clip(y, curve.x[0], curve.x[-1])
        slope = derivative(clipped)
        return curve(clipped) + slope*(y-clipped), slope
    cx, sx = chart(px, dx)
    cz, sz = chart(pz, dz)
    offsets = np.column_stack((p[:, 0]-cx, y-cy, p[:, 2]-cz))
    absolute = np.abs(offsets)
    t = np.clip((absolute-core)/width, 0, 1)
    weights = (1-t)**4*(4*t+1)
    derivatives = np.where((absolute > core) & (absolute < support),
        -20*t*(1-t)**3/width*np.sign(offsets), 0.)
    if y_core_interval_m is not None or y_support_interval_m is not None:
        yc, ys = np.asarray(y_core_interval_m, float), np.asarray(y_support_interval_m, float)
        if yc.shape != (2,) or ys.shape != (2,) or not np.isfinite(np.r_[yc, ys]).all() or not ys[0] < yc[0] <= yc[1] < ys[1]:
            raise ValueError('Ordered finite asymmetric Y support/core intervals required')
        lower_side = y < yc[0]
        gap = np.where(lower_side, yc[0]-ys[0], ys[1]-yc[1])
        distance = np.where(lower_side, yc[0]-y, y-yc[1])
        ty = np.clip(distance/gap, 0, 1)
        weights[:, 1] = (1-ty)**4*(4*ty+1)
        derivatives[:, 1] = np.where((distance > 0) & (distance < gap),
            -20*ty*(1-ty)**3/gap*np.where(lower_side, -1., 1.), 0.)
    weight = np.prod(weights, axis=1)
    gx = derivatives[:, 0]*weights[:, 1]*weights[:, 2]
    gz = weights[:, 0]*weights[:, 1]*derivatives[:, 2]
    gy = weights[:, 0]*derivatives[:, 1]*weights[:, 2] - sx*gx - sz*gz
    gradient = np.column_stack((gx, gy, gz))
    mapped = p + weight[:, None]*d
    jac = np.broadcast_to(np.eye(3), (len(p), 3, 3)).copy() + d[:, :, None]*gradient[:, None, :]
    jac[:, :, 1] += weight[:, None]*derivative_d
    outside = np.any(absolute >= support, axis=1)
    if y_support_interval_m is not None:
        outside = (absolute[:, 0] >= support[0]) | (absolute[:, 2] >= support[2]) | (y <= ys[0]) | (y >= ys[1])
    mapped[outside] = p[outside]
    jac[outside] = np.eye(3)
    return mapped, jac, dict(kind='separable-corridor-translation' if translation_stations is None else 'separable-corridor-course',
        analytic_determinant_lower_bound=float(lower),
        transverse_sections_rigid_inside_XZ_core=True, Y_fixed=True,
        translation_course_C1=True,
        fixed_outside_support=True, discrete_mesh_checks_pending=True)


def transported_corridor_translation(points, *, centerline_stations, core_axes_m,
        support_axes_m, translation_m, y_core_interval_m, y_support_interval_m,
        steps):
    """Compose small positive maps around a moving, rigid transverse core.

    A large one-step translation can compress an adjacent vessel. Transporting
    the compact support with its core permits a much narrower support. Every
    finite step has the same analytic determinant certificate as the static
    corridor; this is a composition, not an approximate ODE integration.
    Y is fixed, so a source transverse section in the core moves rigidly even
    through the longitudinal ramp. No per-object inclusion mask is accepted.
    """
    p = np.asarray(points, float)
    core, support, d = (np.asarray(v, float) for v in (core_axes_m, support_axes_m, translation_m))
    yc, ys = np.asarray(y_core_interval_m, float), np.asarray(y_support_interval_m, float)
    if p.ndim != 2 or p.shape[1] != 3 or any(v.shape != (3,) for v in (core, support, d)):
        raise ValueError('N by 3 points and three-component corridor controls required')
    if (not np.isfinite(np.r_[p.ravel(), core, support, d, yc, ys]).all()
            or np.any(core < 0) or np.any(support <= core) or d[1] != 0
            or yc.shape != (2,) or ys.shape != (2,)
            or not ys[0] < yc[0] <= yc[1] < ys[1]
            or isinstance(steps, bool) or int(steps) != steps or not 1 <= steps <= 256):
        raise ValueError('Finite transverse controls, ordered intervals and 1..256 integer steps required')
    steps = int(steps)
    width = support-core
    step_d = d/steps
    lower_step = 1-2.109375*np.sum(np.abs(step_d[[0, 2]])/width[[0, 2]])
    if lower_step <= 0:
        raise ValueError('Each transported step must certify a positive Jacobian')
    stations = tuple((float(s['y_m']), float(s['x_m']), float(s['z_m'])) for s in centerline_stations)
    px, dx, pz, dz = _corridor_charts(stations)
    y = p[:, 1]
    def chart(curve, derivative):
        clipped = np.clip(y, curve.x[0], curve.x[-1])
        slope = derivative(clipped)
        return curve(clipped)+slope*(y-clipped), slope
    cx, sx = chart(px, dx)
    cz, sz = chart(pz, dz)
    lo = y < yc[0]
    gap = np.where(lo, yc[0]-ys[0], ys[1]-yc[1])
    distance = np.where(lo, yc[0]-y, y-yc[1])
    ty = np.clip(distance/gap, 0, 1)
    wy = (1-ty)**4*(4*ty+1)
    dwy = np.where((distance > 0)&(distance < gap),
        -20*ty*(1-ty)**3/gap*np.where(lo, -1., 1.), 0.)
    # This is the complete swept support, rather than a source-object filter.
    active = ((y > ys[0])&(y < ys[1])
        &(np.abs(p[:, 0]-cx) < support[0]+abs(d[0]))
        &(np.abs(p[:, 2]-cz) < support[2]+abs(d[2])))
    ids = np.flatnonzero(active)
    out = p.copy()
    jac = np.broadcast_to(np.eye(3), (len(p), 3, 3)).copy()
    q = p[ids].copy()
    j = jac[ids].copy()
    for k in range(steps if len(ids) else 0):
        time_fraction = k/steps
        offset = q[:, [0, 2]]-np.c_[cx[ids], cz[ids]]-time_fraction*wy[ids, None]*d[[0, 2]]
        absolute = np.abs(offset)
        t = np.clip((absolute-core[[0, 2]])/width[[0, 2]], 0, 1)
        weights = (1-t)**4*(4*t+1)
        derivatives = np.where((absolute > core[[0, 2]])&(absolute < support[[0, 2]]),
            -20*t*(1-t)**3/width[[0, 2]]*np.sign(offset), 0.)
        weight = weights[:, 0]*weights[:, 1]*wy[ids]
        gx = derivatives[:, 0]*weights[:, 1]*wy[ids]
        gz = weights[:, 0]*derivatives[:, 1]*wy[ids]
        gy = weights[:, 0]*weights[:, 1]*dwy[ids]
        gy -= (sx[ids]+time_fraction*dwy[ids]*d[0])*gx
        gy -= (sz[ids]+time_fraction*dwy[ids]*d[2])*gz
        gradient = np.c_[gx, gy, gz]
        local = np.eye(3)[None]+step_d[None, :, None]*gradient[:, None, :]
        q += weight[:, None]*step_d
        j = local@j
    out[ids], jac[ids] = q, j
    return out, jac, dict(kind='transported-corridor-translation', steps=steps,
        analytic_determinant_lower_bound=float(lower_step**steps),
        analytic_each_step_determinant_lower_bound=float(lower_step),
        finite_composition_not_approximate_flow=True, transverse_sections_rigid_inside_XZ_core=True,
        Y_fixed=True, fixed_outside_swept_support=True, discrete_mesh_checks_pending=True)


def compose_compact_radial_maps(points, fields):
    """Apply ordered common fields and their complete chain-rule Jacobian.

    Each field uses coordinates after the preceding field. No per-object mask
    is accepted: every member of an anatomical closure receives the same map.
    The product of positive determinant bounds certifies the continuous map;
    discrete source triangles and anatomical interfaces still need checking.
    """
    mapped = np.asarray(points, dtype=float)
    if mapped.ndim != 2 or mapped.shape[1] != 3 or not np.isfinite(mapped).all():
        raise ValueError('Finite N by 3 points are required')
    if not isinstance(fields, (list, tuple)) or not fields:
        raise ValueError('A nonempty ordered field sequence is required')
    jacobian = np.broadcast_to(np.eye(3), (len(mapped), 3, 3)).copy()
    proofs = []
    lower = 1.0
    for field in fields:
        if field.get('kind') == 'ellipsoidal-translation':
            mapped, local, proof = compact_ellipsoidal_translation(mapped, **{k:v for k,v in field.items() if k != 'kind'})
        elif field.get('kind') in ['separable-corridor-translation', 'separable-corridor-course']:
            mapped, local, proof = separable_corridor_translation(mapped, **{k:v for k,v in field.items() if k != 'kind'})
        elif field.get('kind') == 'transported-corridor-translation':
            mapped, local, proof = transported_corridor_translation(mapped, **{k:v for k,v in field.items() if k != 'kind'})
        else:
            mapped, local, proof = compact_radial_map(mapped, **field)
        jacobian = local @ jacobian
        lower *= proof['analytic_determinant_lower_bound']
        proofs.append(proof)
    return mapped, jacobian, {
        'kind': 'ordered-composition',
        'field_count': len(fields),
        'analytic_determinant_lower_bound': float(lower),
        'fields': proofs,
        'discrete_mesh_checks_pending': True,
    }
