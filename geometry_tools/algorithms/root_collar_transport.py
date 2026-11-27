"""Source-defined root collar between a fixed attachment and a rigid assembly.

Only source triangles with both root and distal controls are subdivided.
Shared edge/level keys retain exact common controls and face ancestry. A
station's complete inner and outer profiles receive one proper transform.
This does not establish anatomy, tissue thickness, discrete-cell validity,
packing or physiological patency; callers must verify their actual meshes.
"""
import numpy as np
from scipy.spatial.transform import Rotation
from .annular_pipe_envelopes import boundary_loops


def subdivide_root_collar(vertices, faces, root_vertices, levels=(.25,.5,.75)):
    p=np.asarray(vertices,float);f=np.asarray(faces)
    root=np.asarray(root_vertices)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite 3D source controls required')
    if f.ndim!=2 or f.shape[1]!=3 or not np.issubdtype(f.dtype,np.integer)or not len(f)or f.min()<0 or f.max()>=len(p):raise ValueError('Valid indexed source triangles required')
    if root.ndim!=1 or not np.issubdtype(root.dtype,np.integer)or not len(root)or len(set(root))!=len(root)or root.min()<0 or root.max()>=len(p):raise ValueError('Distinct source root indices required')
    levels=np.asarray(levels,float)
    if levels.ndim!=1 or not np.isfinite(levels).all()or np.any(levels<=0)or np.any(levels>=1)or np.any(np.diff(levels)<=0):raise ValueError('Strictly ordered interior stations required')
    station=np.ones(len(p));station[root]=0.
    mixed=np.flatnonzero(np.ptp(station[f],axis=1)>0)
    root_faces=np.flatnonzero(np.all(station[f]==0,axis=1))
    if not len(mixed)or not len(root_faces):raise ValueError('A source attachment and a nonempty transition collar are required')
    points=p.tolist();values=station.tolist();supports=[(i,)for i in range(len(p))]
    parents=[(i,i)for i in range(len(p))];fractions=[0.]*len(p);cache={};output=[];ancestry=[]
    def intersection(a,b,value):
        if values[a]==value:return a
        if values[b]==value:return b
        support=sorted(set(supports[a])|set(supports[b]))
        if len(support)!=2:raise ValueError('A cut intersection left its actual source edge')
        x,y=support;key=(x,y,float(value))
        if key not in cache:
            amount=(value-station[x])/(station[y]-station[x])
            if not 0<amount<1:raise ValueError('Invalid source-edge station interpolation')
            cache[key]=len(points);points.append((p[x]+amount*(p[y]-p[x])).tolist());values.append(float(value));supports.append((x,y));parents.append((x,y));fractions.append(float(amount))
        return cache[key]
    def clip(polygon,value,above):
        result=[]
        for a,b in zip(polygon,polygon[1:]+polygon[:1]):
            ina=values[a]>=value if above else values[a]<=value
            inb=values[b]>=value if above else values[b]<=value
            if ina:result.append(a)
            if ina!=inb:result.append(intersection(a,b,value))
        clean=[]
        for i in result:
            if not clean or i!=clean[-1]:clean.append(i)
        if len(clean)>1 and clean[0]==clean[-1]:clean.pop()
        return clean
    cuts=np.r_[0.,levels,1.]
    mixed_set=set(map(int,mixed))
    for face_id,tri in enumerate(f):
        if face_id not in mixed_set:output.append(tri.tolist());ancestry.append(face_id);continue
        for lo,hi in zip(cuts,cuts[1:]):
            polygon=clip(clip(tri.tolist(),float(lo),True),float(hi),False)
            for k in range(1,len(polygon)-1):output.append([polygon[0],polygon[k],polygon[k+1]]);ancestry.append(face_id)
    q=np.asarray(points,float);ff=np.asarray(output,int);s=np.asarray(values,float);owner=np.asarray(ancestry,int)
    t=q[ff]
    if np.any(np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)==0):raise ValueError('A subdivided source triangle collapsed')
    region=np.isin(owner,np.r_[mixed,root_faces]);profiles=[]
    for level in cuts:
        part=ff[region&np.all(s[ff]<=level,axis=1)]
        loops=boundary_loops(part)
        if len(loops)!=2:raise ValueError('The source collar does not have paired inner/outer station profiles')
        if not all(np.all(s[ring]==level)for ring in loops):raise ValueError('A station profile is not a common source contour')
        profiles.append(dict(station=float(level),rings=loops))
    return dict(vertices=q,faces=ff,station=s,source_face_ids=owner,source_root_face_ids=root_faces,
        source_collar_face_ids=mixed,root_region_faces=np.flatnonzero(region),profiles=profiles,
        interpolation_edges=np.asarray(parents,int),interpolation_fractions=np.asarray(fractions,float),
        original_vertex_count=len(p),levels=cuts)


def transport_root_collar(vertices, station, rotation, translation, center):
    p=np.asarray(vertices,float);s=np.asarray(station,float);r=np.asarray(rotation,float);t=np.asarray(translation,float);c=np.asarray(center,float)
    if p.ndim!=2 or p.shape[1]!=3 or s.shape!=(len(p),)or not np.isfinite(p).all()or not np.isfinite(s).all()or np.any(s<0)or np.any(s>1):raise ValueError('Finite controls with explicit [0,1] stations required')
    if r.shape!=(3,3)or not np.isfinite(r).all()or not np.allclose(r.T@r,np.eye(3),rtol=0,atol=1e-10)or not np.isclose(np.linalg.det(r),1.,rtol=0,atol=1e-10):raise ValueError('A proper rigid distal rotation is required')
    if t.shape!=(3,)or c.shape!=(3,)or not np.isfinite(t).all()or not np.isfinite(c).all():raise ValueError('Finite metre vectors required')
    w=s*s*s*(10+s*(-15+6*s));rotvec=Rotation.from_matrix(r).as_rotvec()
    out=Rotation.from_rotvec(w[:,None]*rotvec).apply(p-c)+c+w[:,None]*t
    out[s==0]=p[s==0]
    return out
