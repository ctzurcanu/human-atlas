"""Section-preserving tube bends, with explicit construction tolerances.

NumPy only. No anatomical dimensions are inferred and no radii are silently
reduced. All interfaces at a station receive the SAME proper rigid transform.

Algorithms independently implemented from:
Wang et al., 2008, Table I (double-reflection rotation-minimizing frames):
https://doi.org/10.1145/1330511.1330513
Yuksel et al., 2009, Fig. 3 (centripetal Catmull-Rom):
https://cemyuksel.com/research/catmullrom_param/catmullrom.pdf
Local spline regularity does not prove that an entire swept tube is disjoint.
"""
from dataclasses import dataclass, asdict
import numpy as np


class TubeGeometryError(ValueError):
    """An invalid request, rather than permission to collapse a tube."""


@dataclass(frozen=True)
class TubeTolerances:
    area_relative: float = 0.01
    diameter_relative: float = 0.01
    thickness_relative: float = 0.01
    absolute_m: float = 1e-7
    planarity_m: float = 1e-4
    curvature_support_limit: float = 0.85
    max_station_turn_rad: float = np.deg2rad(15)
    minimum_station_spacing_m: float = 1e-6
    minimum_wall_thickness_m: float = 1e-6
    maximum_length_change_relative: float | None = None

    def __post_init__(self):
        for name, value in asdict(self).items():
            if value is not None and (not np.isfinite(value) or value < 0):
                raise TubeGeometryError(f'{name} must be finite and nonnegative')
        if not 0 < self.curvature_support_limit < 1:
            raise TubeGeometryError('curvature_support_limit must be between 0 and 1')


def _array(a, name, ndim=None):
    q = np.asarray(a, dtype=float)
    if (ndim is not None and q.ndim != ndim) or not np.isfinite(q).all():
        raise TubeGeometryError(f'{name}: invalid shape or nonfinite coordinates')
    return q


def unit(v):
    v = _array(v, 'direction')
    lengths = np.linalg.norm(v, axis=-1, keepdims=True)
    if np.any(lengths <= 1e-14):
        raise TubeGeometryError('zero direction or degenerate tangent')
    return v / lengths


def arclength(points):
    p = _array(points, 'path', 2)
    if p.shape[1:] != (3,) or len(p) < 2:
        raise TubeGeometryError('path needs at least two 3D points')
    lengths = np.linalg.norm(np.diff(p, axis=0), axis=1)
    if np.any(lengths <= 1e-14):
        raise TubeGeometryError('duplicate consecutive path points')
    return np.r_[0., np.cumsum(lengths)]


def path_tangents(points):
    p = _array(points, 'path', 2)
    s = arclength(p)
    # Physical spacing, rather than station indices, controls derivatives.
    return unit(np.gradient(p, s, axis=0, edge_order=2 if len(p) > 2 else 1))


def resample_path(points, count, tangents=None):
    """Arc-length samples of a polyline; endpoints are retained literally.

    For curved input, first use adaptive centripetal_curve/filleted_path. The
    polyline approximation error is separate from sampling density.
    """
    p = _array(points, 'path', 2); s = arclength(p)
    if not isinstance(count, (int, np.integer)) or count < 2:
        raise TubeGeometryError('sample count must be an integer >= 2')
    target = np.linspace(0, s[-1], count)
    q = np.column_stack([np.interp(target, s, p[:, d]) for d in range(3)])
    q[[0, -1]] = p[[0, -1]]
    t = path_tangents(q) if tangents is None else unit(np.column_stack([
        np.interp(target, s, _array(tangents, 'tangents', 2)[:, d]) for d in range(3)]))
    return q, t


def rotation_minimizing_frames(points, tangents=None, initial_reference=None,
                               final_reference=None):
    """Wang double reflections; frame columns are r, s, tangent.

    An optional final reference imposes the least absolute endpoint roll,
    distributed linearly over physical arclength. Straight portions and
    inflections are regular; zero steps and tangent reversals are rejected.
    """
    p = _array(points, 'path', 2); arc = arclength(p)
    t = path_tangents(p) if tangents is None else unit(tangents)
    if t.shape != p.shape:
        raise TubeGeometryError('one tangent per path station is required')
    if np.any(np.sum(t[:-1]*t[1:], axis=1) < -1+1e-10):
        raise TubeGeometryError('antiparallel tangents: round the turn first')
    ref = np.eye(3)[np.argmin(np.abs(t[0]))] if initial_reference is None else _array(initial_reference, 'initial reference')
    r = unit(ref-t[0]*np.dot(ref,t[0]))
    frames = np.empty((len(p),3,3))
    frames[0] = np.column_stack([r,np.cross(t[0],r),t[0]])
    for i in range(len(p)-1):
        v1 = p[i+1]-p[i]; c1 = np.dot(v1,v1)
        rl = r-2*np.dot(v1,r)/c1*v1
        tl = t[i]-2*np.dot(v1,t[i])/c1*v1
        v2 = t[i+1]-tl; c2 = np.dot(v2,v2)
        if c2 <= 1e-24:
            raise TubeGeometryError('singular frame step: refine or repair the course')
        r = rl-2*np.dot(v2,rl)/c2*v2
        r = unit(r-t[i+1]*np.dot(r,t[i+1]))
        frames[i+1] = np.column_stack([r,np.cross(t[i+1],r),t[i+1]])
    if final_reference is not None:
        ref = _array(final_reference, 'final reference')
        ref = unit(ref-t[-1]*np.dot(ref,t[-1]))
        angle = np.arctan2(frames[-1,:,1]@ref,frames[-1,:,0]@ref)
        roll = angle*arc/arc[-1]
        r,s = frames[:,:,0].copy(),frames[:,:,1].copy()
        frames[:,:,0] = np.cos(roll)[:,None]*r+np.sin(roll)[:,None]*s
        frames[:,:,1] = -np.sin(roll)[:,None]*r+np.cos(roll)[:,None]*s
    return frames


def centripetal_curve(waypoints, chord_tolerance_m=1e-4, max_step_m=0.002,
                     max_subdivisions=20):
    """Adaptively flatten centripetal cubic segments with a control-hull bound.

    Returns points and analytic tangent directions. No local-segment cusp
    guarantee is extrapolated to global tube or whole-curve nonintersection.
    """
    p = _array(waypoints, 'waypoints', 2); arclength(p)
    if chord_tolerance_m <= 0 or max_step_m <= 0:
        raise TubeGeometryError('positive chord tolerance and step required')
    ext = np.r_[2*p[:1]-p[1:2],p,2*p[-1:]-p[-2:-1]]
    points, tangents = [], []
    def flatten(b, depth):
        chord = b[3]-b[0]; ll = chord@chord
        fraction = np.clip((b-b[0])@chord/max(ll,1e-30),0,1)
        error = np.linalg.norm(b-(b[0]+fraction[:,None]*chord),axis=1).max()
        if error <= chord_tolerance_m and np.sqrt(ll) <= max_step_m:
            points.append(b[3]);tangents.append(unit(b[3]-b[2]));return
        if depth >= max_subdivisions:
            raise TubeGeometryError('curve approximation did not meet its error bound')
        a=(b[:-1]+b[1:])/2;c=(a[:-1]+a[1:])/2;m=(c[0]+c[1])/2
        flatten(np.array([b[0],a[0],c[0],m]),depth+1)
        flatten(np.array([m,c[1],a[2],b[3]]),depth+1)
    for i in range(len(p)-1):
        a,b,c,d = ext[i:i+4]
        h0,h1,h2 = np.sqrt(np.linalg.norm(np.diff(ext[i:i+4],axis=0),axis=1))
        m1=h1*((b-a)/h0-(c-a)/(h0+h1)+(c-b)/h1)
        m2=h1*((c-b)/h1-(d-b)/(h1+h2)+(d-c)/h2)
        bezier=np.array([b,b+m1/3,c-m2/3,c])
        if not points:points.append(b.copy());tangents.append(unit(m1))
        flatten(bezier,0)
    return np.array(points),np.array(tangents)


def filleted_path(waypoints, bend_radius_m, max_step_m=0.002,
                  max_turn_rad=np.deg2rad(5)):
    """Circular, tangent-continuous corners; scalar or per-corner bend radius.

    Radius is the centreline radius. Infeasible overlapping trims are rejected
    rather than changing the requested radius or tube section covertly.
    """
    p=_array(waypoints,'waypoints',2);lengths=np.diff(arclength(p))
    if max_step_m <= 0 or not 0 < max_turn_rad < np.pi:
        raise TubeGeometryError('positive spacing and a turn step below pi required')
    radii=np.broadcast_to(np.asarray(bend_radius_m,float),(max(len(p)-2,0),))
    if not np.isfinite(radii).all() or np.any(radii<=0):
        raise TubeGeometryError('positive finite corner radii required')
    directions=unit(np.diff(p,axis=0));corners=[];trims=np.zeros(len(p))
    for i in range(1,len(p)-1):
        u,v=directions[i-1:i+1];angle=float(np.arccos(np.clip(u@v,-1,1)))
        if angle<1e-10:corners.append(None);continue
        if np.pi-angle<1e-8:raise TubeGeometryError('180 degree corner needs a wider explicit course')
        r=radii[i-1];trims[i]=r*np.tan(angle/2)
        axis=unit(np.cross(u,v));start=p[i]-trims[i]*u
        center=start+r*np.cross(axis,u)
        corners.append((start,p[i]+trims[i]*v,center,axis,angle,u,r))
    if np.any(trims[:-1]+trims[1:] >= lengths-1e-12):
        raise TubeGeometryError('requested bend radius does not fit between the waypoints')
    points=[p[0].copy()];tangents=[directions[0].copy()]
    def line(end, tangent):
        start=points[-1];n=int(np.ceil(np.linalg.norm(end-start)/max_step_m))
        for j in range(1,n+1):points.append(start+(end-start)*j/n);tangents.append(tangent)
    for i,corner in enumerate(corners,1):
        if corner is None:line(p[i],directions[i]);continue
        start,end,center,axis,angle,u,r=corner;line(start,u)
        n=max(1,int(np.ceil(max(r*angle/max_step_m,angle/max_turn_rad))))
        radial=start-center
        for a in np.linspace(0,angle,n+1)[1:]:
            def rotate(q):return q*np.cos(a)+np.cross(axis,q)*np.sin(a)+axis*(axis@q)*(1-np.cos(a))
            points.append(center+rotate(radial));tangents.append(rotate(u))
        points[-1]=end
    line(p[-1],directions[-1]);points[-1]=p[-1].copy()
    return np.array(points),unit(np.array(tangents))


def profile_frames(fields):
    """Ordered field shape: [interfaces, stations, section-controls, xyz]."""
    p=_array(fields,'section fields',4)
    if p.shape[-1] != 3 or p.shape[2]<3 or p.shape[1]<2:
        raise TubeGeometryError('at least two stations and three section controls required')
    centers=p[0].mean(axis=1);guide=path_tangents(centers);frames=[]
    for ring,tangent in zip(p[0]-centers[:,None],guide):
        vector=np.cross(ring,np.roll(ring,-1,axis=0)).sum(axis=0)
        normal=unit(vector)
        if normal@tangent<0:normal=-normal
        first=unit(ring[0]-normal*(ring[0]@normal))
        frames.append(np.column_stack([first,np.cross(normal,first),normal]))
    return centers,np.array(frames)


def section_metrics(fields):
    p=_array(fields,'section fields',4);c=p.mean(axis=2,keepdims=True);q=p-c
    area=np.linalg.norm(np.cross(q,np.roll(q,-1,axis=2)).sum(axis=2),axis=2)/2
    edge=np.linalg.norm(np.roll(p,-1,axis=2)-p,axis=3)
    pair=np.linalg.norm(p[:,:,:,None]-p[:,:,None,:],axis=4)
    # Maximum diameter, perimeter and full pairwise shape catch an area-preserving
    # ellipse collapse/shear. Thickness is correspondence distance, not histology.
    return dict(area=area,diameter=pair.max(axis=(2,3)),perimeter=edge.sum(axis=2),
                shape_distances=pair,thickness=np.linalg.norm(np.diff(p,axis=0),axis=3))


def bend_diagnostics(points, fields, tangents=None, limit=0.85):
    c=_array(points,'centers',2);s=arclength(c)
    t=path_tangents(c) if tangents is None else unit(tangents)
    derivative=np.gradient(t,s,axis=0,edge_order=2 if len(c)>2 else 1)
    curvature=derivative-t*np.sum(derivative*t,axis=1)[:,None]
    radial=_array(fields,'fields',4)[0]-c[:,None]
    supports=np.abs(np.einsum('si,svi->sv',curvature,radial)).max(axis=1)
    kappa=np.linalg.norm(curvature,axis=1);radius=np.full(len(c),np.inf)
    np.divide(1,kappa,out=radius,where=kappa>1e-14)
    required_radius=np.zeros(len(c))
    np.divide(supports,kappa*limit,out=required_radius,where=kappa>1e-14)
    angles=np.arccos(np.clip(np.sum(t[:-1]*t[1:],axis=1),-1,1))
    return dict(curvature_per_m=kappa,curvature_support=supports,bend_radius_m=radius,
        estimated_minimum_bend_radius_m=required_radius,
        failed_stations=np.flatnonzero(supports>limit),station_turn_rad=angles,
        minimum_sampled_sweep_jacobian=1-float(supports.max()))


def redirect_local_sections(fields, target_centers, dtype=np.float32):
    """Change a local course while retaining each source section's own roll.

    Minimum proper rotations map the old course tangent onto the new one.
    Sections on unchanged course intervals stay byte-for-byte unchanged.
    This is a construction primitive: it does not establish wall/packing validity.
    Use sweep_sections for a complete rotation-minimizing reconstruction.
    """
    p=_array(fields,'section fields',4)
    if p.shape[-1]!=3 or p.shape[1]<2 or p.shape[2]<3:
        raise TubeGeometryError('at least two stations and three section controls required')
    old=p[0].mean(1);new=_array(target_centers,'target centers',2)
    if new.shape!=old.shape:raise TubeGeometryError('one target per source section required')
    before,after=path_tangents(old),path_tangents(new)
    rotations=np.broadcast_to(np.eye(3),(len(old),3,3)).copy()
    cross=np.cross(before,after);dot=np.sum(before*after,axis=1)
    if np.any(dot<-1+1e-8):
        raise TubeGeometryError('local course reverses: supply a wider explicit path')
    for i,(v,c) in enumerate(zip(cross,dot)):
        x,y,z=v;skew=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
        rotations[i]+=skew+skew@skew/(1+c)
    result=(np.einsum('sij,ksvj->ksvi',rotations,p-old[None,:,None])+new[None,:,None]).astype(dtype)
    same=np.all(new==old,axis=1)&np.all(after==before,axis=1)
    result[:,same]=p[:,same].astype(dtype)
    return result,rotations


def sweep_sections(fields, target_centers, tangents=None, tolerances=None,
                   final_reference=None, dtype=np.float32, preserve_source_roll=False):
    """Transport every wall through one rigid section transform; report faults.

    Does not freeze endpoint rings after transport (that would conceal a bad
    transition). Explicit shared collars need geometrically compatible frames.
    A false report prohibits acceptance; it does not suppress diagnostic output.
    """
    tol=tolerances or TubeTolerances();p=_array(fields,'section fields',4)
    oldc,basis=profile_frames(p);c=_array(target_centers,'target centers',2)
    if c.shape!=oldc.shape:raise TubeGeometryError('one target per source section required')
    t=path_tangents(c) if tangents is None else unit(tangents)
    local=np.einsum('sji,ksvj->ksvi',basis,p-oldc[None,:,None])
    if preserve_source_roll:
        if tangents is not None or final_reference is not None:
            raise TubeGeometryError('local redirection derives course tangents; explicit frame guides require a complete sweep')
        q,rotation=redirect_local_sections(p,c,dtype=dtype)
        frame=rotation@basis
    else:
        frame=rotation_minimizing_frames(c,t,basis[0,:,0],final_reference)
        q=(np.einsum('sij,ksvj->ksvi',frame,local)+c[None,:,None]).astype(dtype)
    before,after=section_metrics(p),section_metrics(q)
    def relative(key):
        if not before[key].size:return 0.
        return float(np.max(np.abs(after[key]-before[key])/np.maximum(np.abs(before[key]),tol.absolute_m**(2 if key=='area' else 1))))
    areaerr,diamerr,thickerr=relative('area'),relative('diameter'),relative('thickness')
    shapeerr=float(np.max(np.abs(after['shape_distances']-before['shape_distances'])))
    planarity=float(np.max(np.abs(local[:,:,:,2])))
    bends=bend_diagnostics(c,q,t,tol.curvature_support_limit)
    spacing=np.diff(arclength(c));lengthchange=abs(arclength(c)[-1]/arclength(oldc)[-1]-1)
    failures=[]
    for bad,reason in [(areaerr>tol.area_relative,'section area changed'),
        (np.any(before['area']<=tol.absolute_m**2),'source section has zero area'),
        (diamerr>tol.diameter_relative,'section diameter changed'),
        (thickerr>tol.thickness_relative,'wall correspondence thickness changed'),
        (len(p)>1 and before['thickness'].min()<tol.minimum_wall_thickness_m,'source wall has zero or insufficient correspondence thickness'),
        (shapeerr>tol.absolute_m,'section shape changed beyond absolute tolerance'),
        (planarity>tol.planarity_m,'source profile is insufficiently planar'),
        (len(bends['failed_stations'])>0,'turn too tight for the actual section'),
        (bends['station_turn_rad'].max()>tol.max_station_turn_rad,'course needs more samples or rounder turns'),
        (spacing.min()<tol.minimum_station_spacing_m,'station spacing too small'),
        (tol.maximum_length_change_relative is not None and lengthchange>tol.maximum_length_change_relative,'course length changed too much')]:
        if bad:failures.append(reason)
    report=dict(section_transport_pass=not failures,failures=failures,tolerances=asdict(tol),
        maximum_relative_area_error=areaerr,maximum_relative_diameter_error=diamerr,
        maximum_relative_thickness_error=thickerr,maximum_shape_distance_error_m=shapeerr,
        maximum_source_planarity_error_m=planarity,maximum_curvature_support=float(bends['curvature_support'].max()),
        tight_bend_stations=bends['failed_stations'].tolist(),minimum_sampled_sweep_jacobian=bends['minimum_sampled_sweep_jacobian'],
        bend_corrections=[dict(station=int(i),estimated_radius_m=float(bends['bend_radius_m'][i]),
            required_radius_m=float(bends['estimated_minimum_bend_radius_m'][i])) for i in bends['failed_stations']],
        minimum_station_spacing_m=float(spacing.min()),maximum_station_turn_rad=float(bends['station_turn_rad'].max()),
        source_length_m=float(arclength(oldc)[-1]),target_length_m=float(arclength(c)[-1]),length_change_relative=float(lengthchange),
        all_wall_interfaces_share_transforms=True,radii_scaled=False,
        frame_method='local minimum-angle redirection' if preserve_source_roll else 'double-reflection rotation-minimizing sweep',
        whole_mesh_intersections_checked=False,anatomical_acceptance=False)
    return q,frame,report
