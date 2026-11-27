"""Boundary-conforming plane splitting of a closed, possibly concave cell.

No anatomical partition is inferred. This is a carrier/LOD preparation helper,
not a reason to call a coarse volume a microscopic lobule. Original controls
remain literal; intersections use canonical shared edges. Concave plane loops
are ear-clipped, never closed with an unverified organ-centre fan. Nested loops
are rejected rather than silently filling a genuine hole.
"""
from dataclasses import dataclass
from collections import Counter, defaultdict
import numpy as np
from .convex_surface_channel import earcut, area


class ClosedCellClipError(ValueError): pass


class SharedCutBank:
    def __init__(self, positions):
        self.source = np.asarray(positions,dtype=np.float64).copy()
        self.positions = self.source.tolist()
        self.cache = {}
        self.intersections = []
        self.planes = {}

    def intersection(self,a,b,normal,offset,plane_id):
        a,b=sorted((int(a),int(b)));key=(a,b,plane_id)
        if key not in self.cache:
            aa,bb=np.asarray(self.positions[a]),np.asarray(self.positions[b])
            da,db=aa@normal-offset,bb@normal-offset
            fraction=da/(da-db)
            if not 0 < fraction < 1:raise ClosedCellClipError('Plane intersection must be on the source edge')
            point=aa+(bb-aa)*fraction
            idx=len(self.positions);self.positions.append(point.tolist());self.cache[key]=idx
            self.intersections.append({'plane':plane_id,'edge':[a,b],'fraction':float(fraction),'control':idx})
        return self.cache[key]


@dataclass(frozen=True)
class ClosedCellSplit:
    negative_faces:np.ndarray
    positive_faces:np.ndarray
    interface_faces:np.ndarray
    negative_ancestry:np.ndarray
    positive_ancestry:np.ndarray
    boundary_loops:tuple
    evidence:dict


def _boundary(faces):
    edges=defaultdict(list)
    for t in faces:
        for a,b in zip(t,np.roll(t,-1)):edges[tuple(sorted((int(a),int(b))))].append((int(a),int(b)))
    if any(len(e)>2 for e in edges.values()):raise ClosedCellClipError('Nonmanifold input or split')
    return [uses[0] for uses in edges.values() if len(uses)==1]


def _inside_polygon(point,poly):
    x,y=point;inside=False
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        if (a[1]>y)!=(b[1]>y):
            crossing=a[0]+(y-a[1])*(b[0]-a[0])/(b[1]-a[1])
            if x<crossing:inside=not inside
    return inside


def _volume(vertices,faces):
    shift=vertices[np.unique(faces)].mean(0);p=vertices[faces]-shift
    return float(np.einsum('ij,ij->i',p[:,0],np.cross(p[:,1],p[:,2])).sum()/6)


def split_closed_cell(bank,faces,normal,offset,plane_id,tolerance=1e-12,*,allow_holes=False):
    """Return two complementary closed boundaries and one shared interface.

    Faces must have consistent outward orientation and a closed boundary.
    The caller owns anatomical labels and must propagate ``bank.intersections``
    to any *other* stored surface sharing the intersected edges. Thus this
    helper does not silently split one member of a multi-tissue assembly.
    """
    point_count=len(bank.positions);intersection_count=len(bank.intersections)
    keys=set(bank.cache);planes=dict(bank.planes)
    try:return _split_closed_cell(bank,faces,normal,offset,plane_id,tolerance,allow_holes)
    except Exception:
        # Rejection must not leave cuts in a caller's shared assembly bank.
        del bank.positions[point_count:];del bank.intersections[intersection_count:]
        for key in bank.cache.keys()-keys:del bank.cache[key]
        bank.planes=planes
        raise


def _split_closed_cell(bank,faces,normal,offset,plane_id,tolerance,allow_holes):
    source_count=len(bank.positions);before=np.asarray(bank.positions,dtype=float).copy()
    faces=np.asarray(faces,dtype=np.int64)
    if _boundary(faces):raise ClosedCellClipError('Input cell is open')
    normal=np.asarray(normal,dtype=float).copy();length=np.linalg.norm(normal)
    if not np.isfinite(length) or length <= 0 or not np.isfinite(offset):
        raise ClosedCellClipError('Plane coefficients must be finite and nonzero')
    normal/=length;offset=float(offset)/length
    plane=(tuple(normal.tolist()),offset)
    if plane_id in bank.planes and bank.planes[plane_id] != plane:
        raise ClosedCellClipError('A shared plane ID cannot name different coefficients')
    bank.planes[plane_id]=plane
    # Local translation avoids cancellation at metre-scale body offsets.
    origin=before[np.unique(faces)].mean(0)
    distance=(before-origin)@normal-(float(offset)-origin@normal)
    if np.min(distance[np.unique(faces)])>=-tolerance or np.max(distance[np.unique(faces)])<=tolerance:
        raise ClosedCellClipError('Plane does not divide the cell')
    negative,positive,na,pa=[],[],[],[]
    for fi,tri in enumerate(faces):
        polys=[]
        for side in (-1,1):
            poly=[]
            for a,b in zip(tri,np.roll(tri,-1)):
                da,db=distance[a],distance[b]
                if side*da>=-tolerance:poly.append(int(a))
                if da*db<0 and abs(da)>tolerance and abs(db)>tolerance:
                    poly.append(bank.intersection(a,b,normal,offset,plane_id))
            clean=[]
            for i in poly:
                if not clean or clean[-1]!=i:clean.append(i)
            if len(clean)>1 and clean[0]==clean[-1]:clean.pop()
            # Oppositely owned copies of a shared face must choose the same
            # quadrilateral diagonal after clipping, regardless of winding.
            if clean:
                first=clean.index(min(clean));clean=clean[first:]+clean[:first]
            out=negative if side==-1 else positive;owners=na if side==-1 else pa
            for k in range(1,len(clean)-1):
                t=[clean[0],clean[k],clean[k+1]]
                if len(set(t))==3:out.append(t);owners.append(fi)
    negative=np.asarray(negative,dtype=np.int64);positive=np.asarray(positive,dtype=np.int64)
    boundary=_boundary(negative)
    follow=defaultdict(list);incoming=Counter()
    for a,b in boundary:follow[a].append(b);incoming[b]+=1
    if not boundary or any(len(bs)!=1 or incoming[a]!=1 for a,bs in follow.items()):
        raise ClosedCellClipError('Cut boundary is not a collection of simple directed loops')
    todo=set(follow);loops=[]
    while todo:
        start=min(todo);loop=[];current=start
        while current in todo:
            todo.remove(current);loop.append(current);current=follow[current][0]
        if current!=start or len(loop)<3:raise ClosedCellClipError('Incomplete cut loop')
        loops.append(loop[::-1]) # Interface faces point out of the negative cell.
    vertices=np.asarray(bank.positions,dtype=float)
    reference=np.eye(3)[int(np.argmin(abs(normal)))];u=np.cross(reference,normal);u/=np.linalg.norm(u);v=np.cross(normal,u)
    controls=set(i for loop in loops for i in loop)
    coordinates={i:np.array([(vertices[i]-origin)@u,(vertices[i]-origin)@v]) for i in controls}
    for i in controls:
        if abs((vertices[i]-origin)@normal-(offset-origin@normal))>2*tolerance:
            raise ClosedCellClipError('An unaccounted input boundary survives the split')
    projected=[np.array([coordinates[i] for i in loop]) for loop in loops]
    if allow_holes:
        # An annular tissue cut preserves the lumen. No new controls are
        # introduced: all cap boundaries are literal carrier intersections.
        import manifold3d
        depths=[sum(_inside_polygon(poly[0],other) for j,other in enumerate(projected) if i!=j)
                for i,poly in enumerate(projected)]
        for poly,depth in zip(projected,depths):
            if area(poly)*(-1 if depth%2 else 1)<=0:
                raise ClosedCellClipError('Cut loop winding disagrees with hole nesting')
        indices=np.concatenate([np.asarray(loop,dtype=np.int64) for loop in loops])
        triangles=np.asarray(manifold3d.triangulate(projected,epsilon=tolerance,allow_convex=False),dtype=np.int64)
        if not len(triangles):raise ClosedCellClipError('Hole triangulation returned no faces')
        caps=indices[triangles]
        cap_xy=np.concatenate(projected)[triangles]
        a,b=cap_xy[:,1]-cap_xy[:,0],cap_xy[:,2]-cap_xy[:,0]
        signed=(a[:,0]*b[:,1]-a[:,1]*b[:,0])/2
        expected_area=sum(area(poly) for poly in projected)
        if np.any(signed<=0) or abs(float(signed.sum())-expected_area)>max(1e-16,expected_area*1e-9):
            raise ClosedCellClipError('Hole triangulation loses or reverses tissue area')
    else:
        for i,poly in enumerate(projected):
            if area(poly)<=0:raise ClosedCellClipError('Unexpected plane-loop orientation or hole')
            if any(_inside_polygon(poly[0],other) for j,other in enumerate(projected) if i!=j):
                raise ClosedCellClipError('Nested cut loops require explicit hole triangulation')
        caps=np.concatenate([np.asarray(earcut(loop,coordinates),dtype=np.int64) for loop in loops])
    negative=np.concatenate([negative,caps]);positive=np.concatenate([positive,caps[:,[0,2,1]]])
    if _boundary(negative) or _boundary(positive):raise ClosedCellClipError('Output cell is not closed')
    for child in (negative,positive):
        directed=Counter((int(a),int(b)) for t in child for a,b in zip(t,np.roll(t,-1)))
        if any(directed[e]!=directed[e[::-1]] for e in directed):raise ClosedCellClipError('Output winding is inconsistent')
    source_volume=_volume(before,faces);nv,pv=_volume(vertices,negative),_volume(vertices,positive)
    if min(nv,pv)<=0:raise ClosedCellClipError('A child has nonpositive volume')
    if abs(nv+pv-source_volume)>max(1e-16,abs(source_volume)*1e-9):raise ClosedCellClipError('Cell volume is not conserved')
    # Convex source-triangle fragments conserve each parent's literal area;
    # the new interface cancels and contributes no extra whole-organ surface.
    original=before[faces];old_area=np.linalg.norm(np.cross(original[:,1]-original[:,0],original[:,2]-original[:,0]),axis=1)/2
    all_faces=np.concatenate([negative[:-len(caps)],positive[:-len(caps)]]);parents=np.asarray(na+pa)
    tri=vertices[all_faces];new_area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    errors=np.bincount(parents,weights=new_area,minlength=len(faces))-old_area
    if np.max(abs(errors))>1e-14:raise ClosedCellClipError('A source face is lost or duplicated')
    evidence={'sourceControlsExact':bool(np.array_equal(vertices[:source_count],before)),
              'sourceVolumeMetersCubed':source_volume,'negativeVolumeMetersCubed':nv,'positiveVolumeMetersCubed':pv,
              'volumeDifferenceMetersCubed':nv+pv-source_volume,'maximumSourceFaceAreaErrorMetersSquared':float(np.max(abs(errors))),
              'interfaceTriangles':len(caps),'simpleCutLoops':len(loops),'oneLiteralSharedInterface':True,
              'explicitHoleTriangulation':bool(allow_holes),
              'unverifiedAnatomicalOwnership':True}
    return ClosedCellSplit(negative,positive,caps,np.r_[na,np.full(len(caps),-1)],np.r_[pa,np.full(len(caps),-1)],tuple(tuple(loop) for loop in loops),evidence)
