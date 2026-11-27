"""Exact tube mouths, annular tissue rims and finite-wall checks.

Matching is limited to explicitly supplied mouth controls. It never welds an
entire body by distance, caps a lumen, or mistakes overlapping tubes for joins.
"""
from collections import Counter, defaultdict
import numpy as np
from .tube_sections import TubeGeometryError, section_metrics, unit


def _faces(faces, count):
    raw=np.asarray(faces)
    if raw.ndim!=2 or raw.shape[1]!=3 or not np.issubdtype(raw.dtype,np.integer):
        raise TubeGeometryError('faces must be an integer [n,3] array')
    if not len(raw) or raw.min()<0 or raw.max()>=count:
        raise TubeGeometryError('face index outside the vertex bank')
    return raw.astype(np.int64,copy=False)


def boundary_loops(faces):
    """Directed face-boundary cycles, retaining supplied control identities."""
    f=np.asarray(faces,dtype=np.int64)
    counts=Counter(tuple(sorted((int(a),int(b)))) for tri in f for a,b in zip(tri,np.roll(tri,-1)))
    if any(n>2 for n in counts.values()):raise TubeGeometryError('nonmanifold carrier edge')
    directed=[(int(a),int(b)) for tri in f for a,b in zip(tri,np.roll(tri,-1)) if counts[tuple(sorted((int(a),int(b))))]==1]
    nexts=defaultdict(list);previous=Counter()
    for a,b in directed:nexts[a].append(b);previous[b]+=1
    if any(len(v)!=1 or previous[k]!=1 for k,v in nexts.items()):
        raise TubeGeometryError('a mouth boundary is branched or inconsistently wound')
    remaining=set(nexts);rings=[]
    while remaining:
        start=min(remaining);ring=[];current=start
        while current in remaining:
            ring.append(current);remaining.remove(current);current=nexts[current][0]
        if current!=start or len(ring)<3:raise TubeGeometryError('a boundary is not a simple closed mouth')
        rings.append(np.array(ring,dtype=np.int64))
    return rings


def continuity_audit(vertices, faces, expected_mouths=2, expected_components=1, expected_euler=None,
                     expected_collars=None, collar_tolerance_m=0):
    """Topology of the actual lumen sheet, using ONLY exact-coordinate welds.

    Disjoint overlapping sheets remain disjoint. Geometry intersection and
    anatomical patency require additional surface/section and lumen checks.
    """
    v=np.asarray(vertices,float)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all():raise TubeGeometryError('invalid vertices')
    f=_faces(faces,len(v));unique,inverse=np.unique(v,axis=0,return_inverse=True);f=inverse[f]
    edges=Counter();directed=Counter();adj=defaultdict(set);links=defaultdict(lambda:defaultdict(set))
    for tri in f:
        for a,b,c in zip(tri,np.roll(tri,-1),np.roll(tri,-2)):
            a,b,c=int(a),int(b),int(c);edges[tuple(sorted((a,b)))]+=1;directed[a,b]+=1
            adj[a].add(b);adj[b].add(a);links[a][b].add(c);links[a][c].add(b)
    groups=[];unseen=set(adj)
    while unseen:
        stack=[min(unseen)];group=set()
        while stack:
            a=stack.pop()
            if a not in unseen:continue
            unseen.remove(a);group.add(a);stack.extend(adj[a])
        groups.append(group)
    badlinks=[]
    for a,graph in links.items():
        todo=[next(iter(graph))];seen=set()
        while todo:
            b=todo.pop()
            if b in seen:continue
            seen.add(b);todo.extend(graph[b]-seen)
        degrees=[len(ns) for ns in graph.values()]
        if len(seen)!=len(graph) or not (all(d==2 for d in degrees) or degrees.count(1)==2 and all(d in (1,2) for d in degrees)):
            badlinks.append(a)
    nonmanifold=sum(n>2 for n in edges.values())
    winding=sum(n==2 and (directed[a,b]!=1 or directed[b,a]!=1) for (a,b),n in edges.items())
    duplicate=len(f)-len({tuple(sorted(t)) for t in f})
    tri=unique[f];degenerate=int(np.sum(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)==0))
    try:rings=boundary_loops(f);boundary_error=None
    except TubeGeometryError as err:rings=[];boundary_error=str(err)
    eulers=[len(g)-sum(a in g for a,b in edges)+sum(int(t[0]) in g for t in f) for g in groups]
    expected_euler=2*expected_components-expected_mouths if expected_euler is None else expected_euler
    collar_errors=[]
    if collar_tolerance_m<0:raise TubeGeometryError('nonnegative collar tolerance required')
    if expected_collars is not None:
        if len(expected_collars)!=expected_mouths:raise TubeGeometryError('provide every expected mouth collar')
        remaining=set(range(len(rings)))
        for expected in expected_collars:
            expected=np.asarray(expected,float)
            if expected.ndim!=2 or expected.shape[1]!=3 or not np.isfinite(expected).all():raise TubeGeometryError('invalid expected collar')
            candidates=[]
            for i in remaining:
                actual=unique[rings[i]]
                if len(actual)!=len(expected):continue
                distances=np.linalg.norm(actual[:,None]-expected[None],axis=2)
                error=float(max(distances.min(axis=0).max(),distances.min(axis=1).max()))
                candidates.append((error,i))
            if not candidates:collar_errors.append(None);continue
            error,index=min(candidates);remaining.remove(index);collar_errors.append(error)
    collars_pass=all(e is not None and e<=collar_tolerance_m for e in collar_errors)
    passed=not any([nonmanifold,winding,duplicate,degenerate,badlinks,boundary_error]) and len(rings)==expected_mouths and len(groups)==expected_components and sum(eulers)==expected_euler and collars_pass
    return dict(topological_continuity_pass=bool(passed),components=len(groups),mouths=len(rings),
        expected_components=expected_components,expected_mouths=expected_mouths,nonmanifold_edges=nonmanifold,
        inconsistent_winding_edges=winding,duplicate_faces=duplicate,degenerate_faces=degenerate,
        invalid_vertex_links=badlinks,boundary_error=boundary_error,
        mouth_details=[dict(controls=len(r),center_m=unique[r].mean(0).tolist(),bounds_m=[unique[r].min(0).tolist(),unique[r].max(0).tolist()]) for r in rings],
        component_euler=eulers,expected_euler=expected_euler,
        expected_collars_checked=expected_collars is not None,collar_maximum_errors_m=collar_errors,
        exact_coordinate_welding_only=True,whole_mesh_intersections_checked=False,anatomical_acceptance=False)


def match_mouths(first, second, position_tolerance_m=1e-6, area_tolerance_relative=.01):
    """Return phase/reversal correspondence minimizing worst ALL-interface gap.

    first/second shape [interfaces, controls,3]. No interpolation, scaling,
    reflection of anatomy, or arbitrary whole-surface weld is performed.
    """
    a,b=np.asarray(first,float),np.asarray(second,float)
    if a.shape!=b.shape or a.ndim!=3 or a.shape[2]!=3 or a.shape[1]<3:
        raise TubeGeometryError('mouths need matching interface/control counts; resample explicitly first')
    if not np.isfinite(a).all() or not np.isfinite(b).all():raise TubeGeometryError('nonfinite mouth')
    if not np.isfinite([position_tolerance_m,area_tolerance_relative]).all() or position_tolerance_m<0 or area_tolerance_relative<0:raise TubeGeometryError('finite nonnegative tolerances required')
    best=None
    for reverse in (False,True):
        base=np.arange(a.shape[1])[::-1] if reverse else np.arange(a.shape[1])
        for phase in range(a.shape[1]):
            order=np.roll(base,phase);error=float(np.linalg.norm(a-b[:,order],axis=2).max())
            if best is None or error<best[0]:best=(error,order,reverse,phase)
    error,order,reverse,phase=best
    ma=section_metrics(a[:,None])['area'];mb=section_metrics(b[:,None])['area']
    if np.any(ma<=0) or np.any(mb<=0):raise TubeGeometryError('mouth contains a zero-area section')
    areaerror=float(np.max(np.abs(ma-mb)/ma))
    if error>position_tolerance_m or areaerror>area_tolerance_relative:
        raise TubeGeometryError(f'mouth mismatch: gap {error*1000:.6g} mm; relative area error {areaerror:.6g}')
    return order,dict(maximum_gap_m=error,maximum_relative_area_error=areaerror,phase=phase,reverse=reverse,
        position_tolerance_m=position_tolerance_m,area_tolerance_relative=area_tolerance_relative)


def join_mouths(first_fields, first_faces, first_mouth, second_fields,
                second_faces, second_mouth, position_tolerance_m=1e-6,
                area_tolerance_relative=.01,tangent_tolerance_rad=np.deg2rad(15)):
    """Shared bank at explicitly matched open ports; first collar is canonical.

    Input fields [interfaces,vertices,3], in the SAME documented frame. Moves
    second collar by at most the supplied tolerance. Unrelated controls stay
    literal. Output has no annulus or disk across the physiological passage.
    """
    a,b=np.asarray(first_fields),np.asarray(second_fields)
    if a.ndim!=3 or b.ndim!=3 or a.shape[0]!=b.shape[0] or a.shape[2]!=3 or b.shape[2]!=3:
        raise TubeGeometryError('both banks need matching interface counts')
    fa,fb=_faces(first_faces,a.shape[1]),_faces(second_faces,b.shape[1])
    am,bm=np.asarray(first_mouth,int),np.asarray(second_mouth,int)
    for mouth,faces,bank in [(am,fa,a),(bm,fb,b)]:
        if mouth.ndim!=1 or len(set(mouth))!=len(mouth) or mouth.min()<0 or mouth.max()>=bank.shape[1]:
            raise TubeGeometryError('invalid mouth controls')
        if not any(set(mouth)==set(r) for r in boundary_loops(faces)):
            raise TubeGeometryError('specified controls are not an entire open boundary')
    order,proof=match_mouths(a[:,am],b[:,bm],position_tolerance_m,area_tolerance_relative)
    def interior_center(vertices,faces,mouth):
        nearby=faces[np.isin(faces,mouth).any(axis=1)]
        neighbors=np.setdiff1d(np.unique(nearby),mouth)
        if not len(neighbors):raise TubeGeometryError('mouth has no adjoining carrier')
        return vertices[neighbors].mean(axis=0)
    ta=unit(a[0,am].mean(axis=0)-interior_center(a[0],fa,am))
    tb=unit(interior_center(b[0],fb,bm)-b[0,bm].mean(axis=0))
    angle=float(np.arccos(np.clip(ta@tb,-1,1)))
    if not 0<=tangent_tolerance_rad<np.pi:raise TubeGeometryError('tangent tolerance must be in [0, pi)')
    if angle>tangent_tolerance_rad:raise TubeGeometryError(f'adjoining course tangents differ by {np.rad2deg(angle):.6g} degrees')
    proof.update(estimated_course_tangent_angle_rad=angle,tangent_tolerance_rad=float(tangent_tolerance_rad))
    mapping=np.full(b.shape[1],-1,dtype=np.int64);mapping[bm[order]]=am
    unused=np.flatnonzero(mapping<0);mapping[unused]=a.shape[1]+np.arange(len(unused))
    bank=np.concatenate([a,b[:,unused]],axis=1);faces=np.r_[fa,mapping[fb]]
    topology=continuity_audit(bank[0],faces)
    if not topology['topological_continuity_pass']:
        raise TubeGeometryError('joined mouths fail carrier topology: '+str(topology))
    proof.update(shared_controls=len(am),interfaces=a.shape[0],central_lumen_capped=False,
        first_collar_literal=True,unrelated_controls_literal=True,topology=topology)
    return bank,faces,mapping,proof


def loft_faces(stations, section_controls):
    """Outward side triangles for positive ordered section winding."""
    if stations<2 or section_controls<3:raise TubeGeometryError('invalid loft dimensions')
    f=[]
    for i in range(stations-1):
        for j in range(section_controls):
            a=i*section_controls+j;b=i*section_controls+(j+1)%section_controls
            c=a+section_controls;d=b+section_controls
            f.extend([[a,b,d],[a,d,c]])
    return np.array(f,dtype=np.int64)


def annular_rim(first_ring, second_ring):
    """Only corresponding wall edges; no central vertex or triangle fan."""
    a,b=np.asarray(first_ring,int),np.asarray(second_ring,int)
    if a.ndim!=1 or a.shape!=b.shape or len(a)<3:raise TubeGeometryError('rim rings must correspond')
    return np.array([[y,x,u] for x,y,u in zip(a,np.roll(a,-1),b)]+[[y,u,v] for y,u,v in zip(np.roll(a,-1),b,np.roll(b,-1))],dtype=np.int64)


def finite_wall_cells(fields, faces):
    """Exact extrema of the straight triangular-prism Jacobian polynomial.

    Outward-oriented tissue and inward interfaces require strictly negative
    determinant at every barycentric corner and throughout the wall depth.
    This local gate is additional to global intersections/containment checks.
    """
    p=np.asarray(fields,float)
    if p.ndim!=3 or p.shape[2]!=3 or len(p)<2 or not np.isfinite(p).all():raise TubeGeometryError('invalid shared wall fields')
    f=_faces(faces,p.shape[1]);failed=[];maxima=[]
    for a,b in zip(p[:-1],p[1:]):
        aa,bb=a[f],b[f];delta=bb-aa
        u,v=aa[:,1]-aa[:,0],aa[:,2]-aa[:,0]
        du,dv=delta[:,1]-delta[:,0],delta[:,2]-delta[:,0]
        c0=np.einsum('fi,fvi->fv',np.cross(u,v),delta)
        c1=np.einsum('fi,fvi->fv',np.cross(du,v)+np.cross(u,dv),delta)
        c2=np.einsum('fi,fvi->fv',np.cross(du,dv),delta)
        maximum=np.maximum(c0,c0+c1+c2);stationary=np.zeros_like(c0)
        np.divide(-c1,2*c2,out=stationary,where=c2!=0)
        valid=(c2!=0)&(stationary>0)&(stationary<1)
        extrema=c0+c1*stationary+c2*stationary**2
        maximum[valid]=np.maximum(maximum[valid],extrema[valid]);m=maximum.max(axis=1)
        failed.append(np.flatnonzero(m>=0).tolist());maxima.append(float(m.max()))
    return dict(finite_wall_cells_pass=not any(failed),invalid_faces_by_interval=failed,
        invalid_counts_by_interval=[len(x) for x in failed],maximum_jacobian_by_interval=maxima,
        whole_mesh_intersections_checked=False,anatomical_acceptance=False)
