"""Closed classifiers for a communicating, role-labelled tubular graph.

This is construction geometry, not an anatomical acceptance certificate.
Cap faces are classification data. Callers must remove classifier closures from
anatomical walls and preserve actual source port rings when attaching parents.
Boolean forks subtract overlapping inner blood spaces rather than leaving
intersecting solid tubes or a membrane at their branch junction.
"""
from __future__ import annotations
import heapq
import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames
from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains,near_boundary_distance


class InteriorRouter:
    """Finite-triangle-constrained A* initializer, followed by checked curves."""
    def __init__(self,triangles,pitch=.0006,clearance=.00022):
        self.triangles=np.asarray(triangles,float);self.pitch=pitch
        low=np.floor(self.triangles.min((0,1))/pitch).astype(int)
        high=np.ceil(self.triangles.max((0,1))/pitch).astype(int)
        index=np.stack(np.meshgrid(*(np.arange(a,b+1) for a,b in zip(low,high)),indexing='ij'),-1).reshape(-1,3)
        points=index*pitch;inside,ambiguous=ray_contains(points,self.triangles)
        # Ambiguous boundary rays do not become trusted interior route nodes.
        distance=near_boundary_distance(points,self.triangles,.002)
        valid=inside&~ambiguous&(distance>clearance)
        self.index=index[valid];self.points=points[valid];self.distance=distance[valid]
        self.by_index={tuple(p):i for i,p in enumerate(self.index)};self.tree=cKDTree(self.points)
        self.steps=np.array([(x,y,z) for x in [-1,0,1] for y in [-1,0,1] for z in [-1,0,1] if (x,y,z)!=(0,0,0)])
        self.step_length=np.linalg.norm(self.steps,axis=1)

    @classmethod
    def from_cache(cls,triangles,cache,pitch=.0006):
        self=cls.__new__(cls);self.triangles=np.asarray(triangles,float);self.pitch=pitch
        z=np.load(cache);self.index=z['index'];self.points=z['points'];self.distance=z['distance']
        self.by_index={tuple(p):i for i,p in enumerate(self.index)};self.tree=cKDTree(self.points)
        self.steps=np.array([(x,y,z) for x in [-1,0,1] for y in [-1,0,1] for z in [-1,0,1] if (x,y,z)!=(0,0,0)])
        self.step_length=np.linalg.norm(self.steps,axis=1);return self

    def route(self,a,b,forbidden=None,start_direction=None,departure_distance=.0015,curve_forbidden=None):
        a,b=np.asarray(a,float),np.asarray(b,float)
        permitted=np.ones(len(self.points),bool) if forbidden is None else ~np.asarray(forbidden,bool)
        # A graph may use a conservative edge-length safety inflation. Curves
        # are checked against their own corridor mask when supplied; keeping
        # the default preserves the previous graph/curve constraint behavior.
        curve_permitted=permitted if curve_forbidden is None else ~np.asarray(curve_forbidden,bool)
        if curve_permitted.shape!=permitted.shape:raise ValueError('Curve corridor mask must match the cached grid')
        pool=np.flatnonzero(permitted);assert len(pool)
        start=pool[np.linalg.norm(self.points[pool]-a,axis=1).argmin()]
        stop=pool[np.linalg.norm(self.points[pool]-b,axis=1).argmin()]
        queue=[(0.,int(start))];cost={int(start):0.};parent={};closed=set()
        while queue:
            _,i=heapq.heappop(queue)
            if i in closed:continue
            closed.add(i)
            if i==stop:break
            for step,length in zip(self.steps,self.step_length):
                j=self.by_index.get(tuple(self.index[i]+step))
                if j is None or not permitted[j]:continue
                # A large clearance penalty selects a central course rather
                # than chasing the shortest shallow subcapsular diagonal.
                value=cost[i]+length*(1+.00022/min(self.distance[j],.002))
                if value>=cost.get(j,np.inf):continue
                cost[j]=value;parent[j]=i
                heuristic=np.linalg.norm(self.points[j]-self.points[stop])/self.pitch
                heapq.heappush(queue,(value+heuristic,j))
        if int(stop) not in cost:raise ValueError('No constrained interior route')
        ids=[int(stop)]
        while ids[-1]!=start:ids.append(parent[ids[-1]])
        raw=np.r_[a[None],self.points[ids[::-1]],b[None]]
        if start_direction is not None:
            far=np.flatnonzero(np.linalg.norm(raw-a,axis=1)>departure_distance)
            if len(far):raw=np.r_[a[None],raw[far[0]:]]
        # Exact endpoints stay fixed. Round only newly constructed controls.
        raw=raw[np.r_[True,np.linalg.norm(np.diff(raw,axis=0),axis=1)>1e-10]]
        rounded=gaussian_filter1d(raw,.8,axis=0,mode='nearest');rounded[0]=a;rounded[-1]=b
        for points in [rounded,raw]:
            s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
            if len(s)<3:return points
            t=np.linspace(0,s[-1],max(6,int(np.ceil(s[-1]/.00018))))
            bc='natural' if start_direction is None else ((1,np.asarray(start_direction,float)),(2,np.zeros(3)))
            curve=CubicSpline(s,points,bc_type=bc)(t)
            inside,ambiguous=ray_contains(curve,self.triangles)
            distances=near_boundary_distance(curve,self.triangles,.0003)
            if np.all(inside&~ambiguous&(distances>.00016)):
                if (forbidden is None and curve_forbidden is None) or np.all(curve_permitted[self.tree.query(curve)[1]]):return curve
        raise ValueError('Rounded route does not retain its constrained clearance')

    def exclude_segments(self,segments,radii,margin):
        forbidden=np.zeros(len(self.points),bool)
        for (a,b),r in zip(np.asarray(segments,float),np.asarray(radii,float)):
            delta=b-a;length=np.linalg.norm(delta)
            ids=np.array(self.tree.query_ball_point((a+b)/2,length/2+r+margin),int)
            if not len(ids):continue
            t=np.clip((self.points[ids]-a)@delta/(delta@delta),0,1)
            forbidden[ids[np.linalg.norm(self.points[ids]-a-t[:,None]*delta,axis=1)<=r+margin]]=True
        return forbidden


def tube_classifier(path,radii,ring_size=12,start_ring=None,end_ring=None,extend_joints=True):
    """Closed swept polygon volume and explicit closure-face IDs.

    Exact input rings can be carried literally at a source interface. New
    interior graph joints overlap by 1.15 local radii, so union has an actual
    aperture; a coincidence of rim centers never substitutes for a join.
    """
    path=np.asarray(path,float).copy();radii=np.broadcast_to(np.asarray(radii,float),len(path)).copy()
    tangent=np.gradient(path,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    if extend_joints:
        if start_ring is None:path[0]-=tangent[0]*radii[0]*1.15
        if end_ring is None:path[-1]+=tangent[-1]*radii[-1]*1.15
    reference=None if start_ring is None else np.asarray(start_ring)[0]-path[0]
    frames=rotation_minimizing_frames(path,tangent,initial_reference=reference)
    angle=np.arange(ring_size)*2*np.pi/ring_size
    circle=np.c_[np.cos(angle),np.sin(angle),np.zeros(ring_size)]
    if start_ring is not None:
        # Source knots may be nonuniform or run in the opposite direction.
        # Sweep the literal polygon instead of mismatching it to a new circle.
        local=(np.asarray(start_ring,float)-path[0])@frames[0]
        circle=local/radii[0]
    sections=np.einsum('kij,nj->kni',frames,circle)*radii[:,None,None]+path[:,None,:]
    if start_ring is not None:
        ring=np.asarray(start_ring,float);assert len(ring)==ring_size;sections[0]=ring
    if end_ring is not None:
        ring=np.asarray(end_ring,float);assert len(ring)==ring_size;sections[-1]=ring
    v=sections.reshape(-1,3);f=[]
    for k in range(len(path)-1):
        for j in range(ring_size):
            a=k*ring_size+j;b=k*ring_size+(j+1)%ring_size;c=b+ring_size;d=a+ring_size
            f.extend([[a,b,c],[a,c,d]])
    caps=[]
    for k,reverse in [(0,True),(len(path)-1,False)]:
        # A fan about a separate point is classifier-only, never wall tissue.
        center=len(v);v=np.r_[v,sections[k].mean(0)[None]];faces=[]
        for j in range(ring_size):
            tri=[center,k*ring_size+j,k*ring_size+(j+1)%ring_size]
            faces.append(tri[::-1] if reverse else tri)
        caps.append(np.arange(len(f),len(f)+ring_size));f.extend(faces)
    f=np.array(f,np.uint32)
    volume=np.einsum('ij,ij->i',v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]])).sum()/6
    if volume<0:f=f[:,::-1]
    return v,f,caps


def boolean_union(tubes,origin):
    """Float64 millimetre-local CSG with original per-triangle provenance."""
    import manifold3d as m
    origin=np.asarray(origin,float);solids=[];provenance=[];offset=0
    for role,v,f,caps in tubes:
        ids=np.arange(offset,offset+len(f),dtype=np.uint64);offset+=len(f)
        solid=m.Manifold(m.Mesh64((np.asarray(v,float)-origin)*1000,np.asarray(f,np.uint64),face_id=ids))
        if solid.status()!=m.Error.NoError:raise ValueError((role,str(solid.status())))
        solids.append(solid);provenance.append({'role':role,'firstSourceFaceId':int(ids[0]),'sourceFaceCount':len(f),'classifierSourceCapFaceIds':[ids[c].astype(int).tolist() for c in caps]})
    joined=m.Manifold.batch_boolean(solids,m.OpType.Add)
    if joined.status()!=m.Error.NoError:raise ValueError(str(joined.status()))
    mesh=joined.to_mesh64();vertices=np.asarray(mesh.vert_properties)[:,:3]/1000+origin
    return vertices,np.asarray(mesh.tri_verts,np.uint32),np.asarray(mesh.face_id,np.uint64),provenance,joined
