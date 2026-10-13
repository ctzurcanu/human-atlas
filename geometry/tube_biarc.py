"""Tangent-continuous circular 3D tube courses with explicit radius bounds.

Independent NumPy implementation of the biarc equations described by
Ryan Juckett (2014), https://www.ryanjuckett.com/biarc-interpolation/ , and
Song et al. (2014), doi:10.2312/pgs.20141248, section 2, equation 4.
No source code copied. Tangents are never reversed to conceal a bad port.
A bounded course is not a wall, self-intersection or body-clearance certificate.
"""
import numpy as np
from .tube_sections import TubeGeometryError, unit


def _arc(start, tangent, end, max_step, max_turn, chord_error):
    delta=end-start;length=np.linalg.norm(delta)
    if length<=1e-14:raise TubeGeometryError('zero-length arc')
    along=float(delta@tangent);side=delta-along*tangent;squared=float(side@side)
    if squared<=length*length*1e-24:
        if along<=0:raise TubeGeometryError('backward collinear arc requires an explicit wider course')
        n=max(1,int(np.ceil(length/max_step)))
        points=start+np.linspace(0,1,n+1)[:,None]*delta
        tangents=np.broadcast_to(tangent,points.shape).copy()
        return points,tangents,np.zeros_like(points),dict(radius_m=None,angle_rad=0.,length_m=float(length))
    center=start+side*(length*length/(2*squared));radius=float(np.linalg.norm(start-center))
    radial=(start-center)/radius;last=(end-center)/radius
    angle=float(np.arctan2(last@tangent,last@radial)%(2*np.pi))
    sagitta_turn=2*np.arccos(np.clip(1-chord_error/radius,-1,1))
    if sagitta_turn<=0:raise TubeGeometryError('chord error too small for floating-point arc sampling')
    n=max(1,int(np.ceil(radius*angle/max_step)),int(np.ceil(angle/min(max_turn,sagitta_turn))))
    theta=np.linspace(0,angle,n+1)[:,None]
    direction=np.cos(theta)*radial+np.sin(theta)*tangent
    points=center+radius*direction;tangents=-np.sin(theta)*radial+np.cos(theta)*tangent
    points[[0,-1]]=[start,end];tangents[0]=tangent
    return points,tangents,-direction/radius,dict(radius_m=radius,angle_rad=angle,length_m=radius*angle)


def biarc_path(start, start_tangent, end, end_tangent, *,
               first_tangent_distance_m=None, minimum_bend_radius_m=0.,
               max_step_m=.002, max_turn_rad=np.deg2rad(5), chord_error_m=1e-5):
    """Return points, analytic tangents/curvature and construction diagnostics.

    A positive custom first tangent distance controls the arc balance. The
    default uses equal distances. Equal tangents perpendicular to the chord
    use two semicircles. Infeasible/degenerate data raise; no tube is thinned.
    Radius bounds apply to both whole circular arcs, not just sampled curvature.
    """
    p0,p1=np.asarray(start,float),np.asarray(end,float)
    t0,t1=unit(start_tangent),unit(end_tangent)
    if any(q.shape!=(3,) or not np.isfinite(q).all() for q in (p0,p1,t0,t1)):
        raise TubeGeometryError('finite 3D endpoints and tangents required')
    if not np.isfinite([minimum_bend_radius_m,max_step_m,max_turn_rad,chord_error_m]).all() or minimum_bend_radius_m<0 or max_step_m<=0 or not 0<max_turn_rad<np.pi or chord_error_m<=0:
        raise TubeGeometryError('finite nonnegative radius and positive sampling tolerances required')
    v=p1-p0;vv=float(v@v)
    if vv<=1e-28:raise TubeGeometryError('distinct endpoints required')
    scale=np.sqrt(vv);dot=float(np.clip(t0@t1,-1,1));a=2*(1-dot);b=float(v@(t0+t1))
    if first_tangent_distance_m is not None:
        d0=float(first_tangent_distance_m)
        if not np.isfinite(d0) or d0<=0:raise TubeGeometryError('positive first tangent distance required')
        denominator=float(v@t1+d0*(1-dot))
        if abs(denominator)<=scale*1e-14:raise TubeGeometryError('singular biarc shape parameter; choose another distance')
        d1=(vv/2-d0*float(v@t0))/denominator
        if d1<=0 or not np.isfinite(d1):raise TubeGeometryError('shape parameter creates a reversed or degenerate join')
        middle=((p0+d0*t0)*d1+(p1-d1*t1)*d0)/(d0+d1)
    elif a<1e-12:
        projection=float(v@t0)
        if abs(projection)<=scale*1e-12:
            # Both equal tangent lines are perpendicular to the chord.
            d0=d1=None;middle=(p0+p1)/2
        elif projection<0:
            raise TubeGeometryError('equal tangents point behind the end port; supply a wider explicit course')
        else:
            d0=d1=vv/(4*projection);middle=(p0+p1)/2
    else:
        root=np.sqrt(b*b+a*vv)
        d0=d1=vv/(root+b) if b>=0 else (root-b)/a
        middle=(p0+p1+d0*(t0-t1))/2
    first=_arc(p0,t0,middle,max_step_m,max_turn_rad,chord_error_m)
    reverse=_arc(p1,-t1,middle,max_step_m,max_turn_rad,chord_error_m)
    pp,tt,kk,one=first;rp,rt,rk,two=reverse
    join_error=float(np.linalg.norm(tt[-1]+rt[-1]))
    if join_error>1e-9:raise TubeGeometryError('circular arcs fail their shared tangent')
    for part in (one,two):
        if part['radius_m'] is not None and part['radius_m']<minimum_bend_radius_m:
            raise TubeGeometryError(f"arc radius {part['radius_m']:.9g} m is below requested {minimum_bend_radius_m:.9g} m")
    points=np.r_[pp,rp[::-1][1:]];tangents=np.r_[tt,-rt[::-1][1:]];curvature=np.r_[kk,rk[::-1][1:]]
    points[[0,-1]]=[p0,p1];tangents[[0,-1]]=[t0,t1]
    return points,tangents,curvature,dict(arcs=[one,two],join_station=len(pp)-1,
        tangent_join_error=join_error,first_tangent_distance_m=d0,second_tangent_distance_m=d1,
        exact_length_m=one['length_m']+two['length_m'],minimum_bend_radius_m=minimum_bend_radius_m,
        max_step_m=max_step_m,max_turn_rad=max_turn_rad,chord_error_m=chord_error_m,
        tangents_flipped=False,radii_scaled=False,whole_mesh_intersections_checked=False,
        anatomical_acceptance=False)
