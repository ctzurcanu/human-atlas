"""Order-preserving, continuous bounds for differently subdivided closed rims.

The two polylines receive a monotone normalized arc-length correspondence.
Every interval is split at both curves' corners. The difference is affine on
each resulting interval, so its norm attains its maximum at an endpoint.
This proves a whole-curve distance bound, rather than sampling a vertex set.
It neither joins surfaces nor establishes an anatomical attachment.
"""
import numpy as np


def _curve(points):
    p = np.asarray(points, float)
    if p.ndim != 2 or p.shape[1] != 3 or not 3 <= len(p) <= 256 or not np.isfinite(p).all():
        raise ValueError('Three to 256 finite ordered rim points required')
    lengths = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    if np.any(lengths <= 0) or len(np.unique(p, axis=0)) != len(p):
        raise ValueError('A rim must have distinct points and positive edges')
    return p, np.r_[0., np.cumsum(lengths)] / lengths.sum()


def _at(p, knots, t):
    t = np.mod(np.asarray(t), 1.)
    i = np.minimum(np.searchsorted(knots, t, side='right') - 1, len(p) - 1)
    f = (t - knots[i]) / (knots[i + 1] - knots[i])
    return p[i] + f[..., None] * (p[(i + 1) % len(p)] - p[i])


def closed_curve_correspondence(first, second, *, opposite=True):
    """Return a conservative Float64 bound for a cyclic monotone pairing.

    Reversal is explicit: opposite source face winding is checked separately
    by the caller. Candidate phases project one fixed corner onto each edge;
    finding a passing phase is sufficient, and failure is not proof that no
    other monotone correspondence could pass.
    """
    p, a = _curve(first)
    q, b = _curve(np.asarray(second)[::-1] if opposite else second)
    edge = np.roll(q, -1, axis=0) - q
    f = np.clip(np.einsum('ij,ij->i', p[0] - q, edge) /
                np.einsum('ij,ij->i', edge, edge), 0., 1.)
    phases = np.unique(np.mod(b[:-1] + f * np.diff(b), 1.))
    margin = float(512 * np.finfo(float).eps * max(1., np.abs(p).max(), np.abs(q).max()))
    best = None
    for phase in phases:
        knots = np.unique(np.r_[0., 1., a[:-1], np.mod(b[:-1] - phase, 1.)])
        first_points, second_points = _at(p, a, knots), _at(q, b, knots + phase)
        maximum = float(np.linalg.norm(first_points - second_points, axis=1).max())
        row = dict(maximumGap_m=maximum + margin, maximumEndpointGap_m=maximum,
                   floatingPointMargin_m=margin, phase=float(phase),
                   commonIntervals=len(knots) - 1, firstControls=len(p), secondControls=len(q),
                   orderPreserved=True, wholeCurveBound=True,
                   method='monotone cyclic arc-length; convex norm on every common linear interval')
        if best is None or row['maximumGap_m'] < best['maximumGap_m']:
            best = row
    return best


def align_existing_subdivision(fine, coarse, *, maximum_motion_m):
    """Align an explicit coarse rim with a subdivided copy of the same rim.

    Keep every fine corner corresponding to a coarse corner fixed. Move only
    the extra fine edge points onto those unchanged polygon edges and the
    coarse corners onto their uniquely corresponding fine corners. No control
    is inserted/deleted and no triangle or cap is constructed.
    """
    fine,_=_curve(fine);coarse,_=_curve(coarse)
    if not len(fine)>len(coarse) or not 0<maximum_motion_m<=1e-5:
        raise ValueError('Explicit finer/coarser rims and a bound at most 10 micrometres required')
    distance=np.linalg.norm(coarse[:,None]-fine[None],axis=2)
    corner=np.argmin(distance,axis=1)
    if len(np.unique(corner))!=len(coarse) or distance[np.arange(len(coarse)),corner].max()>maximum_motion_m:
        raise ValueError('Every coarse corner needs a unique existing fine corner within the bound')
    # A matching source interface has opposite winding. Nearest assignments
    # may select the phase, but cannot reorder or skip a coarse corner.
    order=sorted(corner.tolist())
    source=[int(np.flatnonzero(corner==i)[0])for i in order]
    if not all((source[(i+1)%len(source)]-source[i])%len(source)==len(source)-1 for i in range(len(source))):
        raise ValueError('Declared corners do not form opposite cyclic source orders')
    target_fine=fine.copy();target_coarse=fine[corner].copy();extra=[]
    for i in range(len(fine)):
        if i in corner:continue
        left=max((j for j in corner if j<i),default=max(corner))
        right=min((j for j in corner if j>i),default=min(corner))
        edge=fine[right]-fine[left]
        fraction=float((fine[i]-fine[left])@edge/(edge@edge))
        if not 0<fraction<1:raise ValueError('A fine subdivision point is not inside its declared edge')
        target_fine[i]=fine[left]+fraction*edge
        extra.append(dict(fineControl=i,firstCorner=int(left),secondCorner=int(right),fraction=fraction))
    movement=max(float(np.linalg.norm(target_fine-fine,axis=1).max()),float(np.linalg.norm(target_coarse-coarse,axis=1).max()))
    if movement>maximum_motion_m:raise ValueError('Subdivision alignment exceeds the explicit motion bound')
    proof=closed_curve_correspondence(target_fine,target_coarse)
    if proof['maximumGap_m']>1e-10:raise ValueError('Aligned complete polylines do not coincide')
    return target_fine,target_coarse,dict(fineCornerForCoarse=corner.tolist(),extraFineControls=extra,
        maximumMotion_m=movement,fineCornersLiteral=True,facesAndControlCountsUnchanged=True,
        completeCurveProof=proof,anatomicalAcceptance=False)
