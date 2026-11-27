"""Persistent triangle BVHs and PostgreSQL 3D broad-phase queries."""
from __future__ import annotations
from pathlib import Path
import heapq,time,os
import numpy as np
from .common import CACHE,VERSION,identity,digest_file,atomic_npz,owned_path
from .distance import point_triangle,triangle_pairs


def box_distance_squared(low_a,high_a,low_b,high_b):
    delta=np.maximum(np.maximum(np.asarray(low_a)-high_b,np.asarray(low_b)-high_a),0)
    return np.einsum('...i,...i->...',delta,delta)


class BVH:
    def __init__(self,mesh,leaf_size=12,arrays=None):
        if leaf_size<2 or leaf_size>128:raise ValueError('Leaf size must be between 2 and 128')
        if not len(mesh.faces):raise ValueError('A triangle BVH requires a nonempty mesh')
        self.mesh=mesh;self.leaf_size=leaf_size;self.triangles=mesh.triangles
        if arrays is not None:
            for key in ['low','high','children','starts','counts','order']:setattr(self,key,np.asarray(arrays[key]))
            if sorted(self.order.tolist())!=list(range(len(mesh.faces))):raise ValueError('Cached BVH does not cover every source triangle exactly once')
            return
        t=self.triangles;lo=t.min(1);hi=t.max(1);centers=(lo+hi)*.5;low=[];high=[];children=[];starts=[];counts=[];order=[]
        def node(ids):
            index=len(low);low.append(lo[ids].min(0));high.append(hi[ids].max(0));children.append([-1,-1]);starts.append(-1);counts.append(len(ids))
            if len(ids)<=leaf_size:
                starts[index]=len(order);order.extend(ids.tolist())
            else:
                axis=int(np.argmax(np.ptp(centers[ids],axis=0)));ids=ids[np.argsort(centers[ids,axis],kind='stable')];mid=len(ids)//2
                children[index]=[node(ids[:mid]),node(ids[mid:])]
            return index
        node(np.arange(len(t)));self.low=np.array(low);self.high=np.array(high);self.children=np.array(children,dtype='<i8');self.starts=np.array(starts,dtype='<i8');self.counts=np.array(counts,dtype='<i8');self.order=np.array(order,dtype='<i8')

    def save(self,path):
        return atomic_npz(path,geometry_sha=np.array(self.mesh.geometry_sha),leaf_size=np.array(self.leaf_size),**{k:getattr(self,k)for k in ['low','high','children','starts','counts','order']})

    @classmethod
    def load(cls,mesh,path):
        with np.load(path,allow_pickle=False)as z:
            if str(z['geometry_sha'])!=mesh.geometry_sha:raise ValueError('A cached BVH belongs to different actual triangles')
            return cls(mesh,int(z['leaf_size']),{k:z[k].copy()for k in ['low','high','children','starts','counts','order']})

    def leaf(self,node):return self.children[node,0]<0
    def ids(self,node):return self.order[self.starts[node]:self.starts[node]+self.counts[node]]

    def query_bounds(self,low,high,*,distance_m=0.):
        """Return every literal triangle box within a finite query box.

        The BVH is immutable. A bounded deformation can conservatively query
        this source index with an explicit displacement margin.
        """
        low=np.asarray(low,float);high=np.asarray(high,float)
        if low.shape!=(3,)or high.shape!=(3,)or not np.isfinite(np.r_[low,high]).all()or np.any(low>high):raise ValueError('Ordered finite metre bounds required')
        if not np.isfinite(distance_m)or distance_m<0:raise ValueError('Finite nonnegative query margin required')
        stack=[0];found=[];squared=distance_m*distance_m
        while stack:
            node=stack.pop()
            if box_distance_squared(low,high,self.low[node],self.high[node])>squared:continue
            if self.leaf(node):
                ids=self.ids(node);t=self.triangles[ids]
                found.extend(ids[box_distance_squared(low,high,t.min(1),t.max(1))<=squared].tolist())
            else:stack.extend(map(int,self.children[node]))
        return np.asarray(sorted(found),dtype=np.int64)

    def closest_point(self,point,*,seed_face=None):
        point=np.asarray(point,float)
        if point.shape!=(3,)or not np.isfinite(point).all():raise ValueError('A finite metre point is required')
        queue=[(float(box_distance_squared(point,point,self.low[0],self.high[0])),0)];best=np.inf;result=None;tested=0
        if seed_face is not None:
            if not isinstance(seed_face,(int,np.integer))or not 0<=seed_face<len(self.triangles):raise ValueError('Seed must identify a literal source triangle')
            d,q,w=point_triangle(point[None],self.triangles[seed_face][None]);best=float(d[0]);tested=1
            result=dict(point=q[0],face_id=int(seed_face),barycentric=w[0])
        while queue:
            lower,node=heapq.heappop(queue)
            if lower>best:break
            if self.leaf(node):
                ids=self.ids(node);d,q,w=point_triangle(np.broadcast_to(point,(len(ids),3)),self.triangles[ids]);tested+=len(ids);j=int(np.argmin(d))
                if d[j]<best:best=float(d[j]);result=dict(point=q[j],face_id=int(ids[j]),barycentric=w[j])
                if best==0:break
            else:
                for child in self.children[node]:
                    lower=float(box_distance_squared(point,point,self.low[child],self.high[child]))
                    if lower<=best:heapq.heappush(queue,(lower,int(child)))
        return {**result,'distance_m':float(np.sqrt(best)),'tested_triangles':tested}

    def closest_points(self,points,*,backend='bvh'):
        """Batch queries, optionally native-seeded; every result gets BVH search."""
        points=np.asarray(points,float)
        if points.ndim!=2 or points.shape[1]!=3 or not np.isfinite(points).all():raise ValueError('Finite metre points with shape [N,3] required')
        if backend not in ['auto','open3d','bvh']:raise ValueError('Distance backend must be auto, open3d or bvh')
        if backend=='open3d'or backend=='auto'and len(points)>=32:
            from .open3d_backend import SurfaceIndex,available
            if backend=='open3d'or available():
                try:
                    if not hasattr(self,'_native_surface'):self._native_surface=SurfaceIndex(self.mesh)
                    return self._native_surface.closest_points(points,certify=True,bvh=self)
                except ValueError:
                    # Thin triangles may collapse in the native Float32 copy.
                    # The Float64 index still covers every original triangle.
                    if backend=='open3d':raise
        rows=[self.closest_point(p)for p in points]
        arrays=dict(distance_m=np.array([r['distance_m']for r in rows]),
                    points=np.array([r['point']for r in rows]).reshape(-1,3),
                    face_ids=np.array([r['face_id']for r in rows],dtype=np.int64),
                    barycentric=np.array([r['barycentric']for r in rows]).reshape(-1,3),
                    tested_triangles=np.array([r['tested_triangles']for r in rows],dtype=np.int64))
        return arrays,dict(backend='bvh',geometry_sha256=self.mesh.geometry_sha,frame=self.mesh.frame,
                           queries=len(points),nearest_triangle_verified=True,float64_source_witnesses=True,
                           distance_is_upper_bound_only=False,anatomical_acceptance=False)

    def candidate_batches(self,other,distance_m,*,selfcheck=False,batch_size=2048,root_pair=None):
        if distance_m<0 or not np.isfinite(distance_m):raise ValueError('A finite nonnegative proximity distance is required')
        if selfcheck and self.mesh.geometry_sha!=other.mesh.geometry_sha:raise ValueError('Self query requires the same exact triangle sequence')
        if root_pair is None:root_pair=(0,0)
        if len(root_pair)!=2 or any(not isinstance(i,(int,np.integer))for i in root_pair)or not 0<=root_pair[0]<len(self.low)or not 0<=root_pair[1]<len(other.low):raise ValueError('A pair of existing BVH subtree nodes is required')
        squared=distance_m*distance_m;stack=[tuple(map(int,root_pair))];buffer=[]
        while stack:
            a,b=stack.pop()
            if box_distance_squared(self.low[a],self.high[a],other.low[b],other.high[b])>squared:continue
            if self.leaf(a)and other.leaf(b):
                ai,bi=np.meshgrid(self.ids(a),other.ids(b),indexing='ij');pairs=np.c_[ai.ravel(),bi.ravel()]
                if selfcheck and a==b:pairs=pairs[pairs[:,0]<pairs[:,1]]
                if selfcheck and a!=b:pairs=np.sort(pairs,axis=1)
                if not len(pairs):continue
                ta,tb=self.triangles[pairs[:,0]],other.triangles[pairs[:,1]]
                take=box_distance_squared(ta.min(1),ta.max(1),tb.min(1),tb.max(1))<=squared;pairs=pairs[take]
                buffer.extend(pairs.tolist())
                while len(buffer)>=batch_size:
                    yield np.asarray(buffer[:batch_size],dtype='<i8');del buffer[:batch_size]
            elif selfcheck and a==b:
                left,right=self.children[a];stack.extend([(int(left),int(left)),(int(left),int(right)),(int(right),int(right))])
            elif not self.leaf(a)and(other.leaf(b)or self.counts[a]>=other.counts[b]):
                stack.extend((int(child),b)for child in self.children[a])
            else:stack.extend((a,int(child))for child in other.children[b])
        if buffer:yield np.asarray(buffer,dtype='<i8')

    def closest_mesh(self,other):
        queue=[(float(box_distance_squared(self.low[0],self.high[0],other.low[0],other.high[0])),0,0)];best=np.inf;result=None;tested=0
        while queue:
            lower,a,b=heapq.heappop(queue)
            if lower>best:break
            if self.leaf(a)and other.leaf(b):
                ai,bi=np.meshgrid(self.ids(a),other.ids(b),indexing='ij');pairs=np.c_[ai.ravel(),bi.ravel()];r=triangle_pairs(self.triangles[pairs[:,0]],other.triangles[pairs[:,1]]);tested+=len(pairs);i=int(np.argmin(r['distance_m']));d=float(r['distance_m'][i])**2
                if d<best:best=d;result={k:v[i]for k,v in r.items()};result.update(first_face_id=int(pairs[i,0]),second_face_id=int(pairs[i,1]))
                if best==0:break
            else:
                choices=[(int(child),b)for child in self.children[a]]if not self.leaf(a)and(other.leaf(b)or self.counts[a]>=other.counts[b])else[(a,int(child))for child in other.children[b]]
                for x,y in choices:
                    lower=float(box_distance_squared(self.low[x],self.high[x],other.low[y],other.high[y]))
                    if lower<=best:heapq.heappush(queue,(lower,x,y))
        result['tested_triangle_pairs']=tested
        return result


def cached_bvh(db,mesh,leaf_size=12):
    parameters=dict(leaf_size=leaf_size,numpy=np.__version__,implementation_sha256=digest_file(__file__))
    key=identity(VERSION,'bvh',mesh.geometry_sha,parameters);record=db.cache_get(key)
    if record:return BVH.load(mesh,record['path']),dict(cache_hit=True,key=key)
    db.ensure_geometry(mesh)
    bvh=BVH(mesh,leaf_size);path=CACHE/'bvh'/(key+'.npz')
    if path.exists():
        owned_path(path);os.replace(path,path.with_suffix('.corrupt-'+str(time.time_ns())+'.npz'))
    bvh.save(path);db.cache_put(key=key,geometry_sha=mesh.geometry_sha,kind='triangle-bvh',parameters=parameters,path=path,summary=dict(triangles=len(mesh.faces),nodes=len(bvh.low)))
    return bvh,dict(cache_hit=False,key=key)


def neighbors(db,reference,*,model=None,distance_m=.002,include_hidden=False,limit=1000):
    row=db.mesh(reference,model)
    if not np.isfinite(distance_m)or distance_m<0:raise ValueError('Nonnegative finite clearance required')
    if not 1<=limit<=100000:raise ValueError('Invalid neighbour limit')
    with db.connect()as c:
        if row['bbox']is None:raise ValueError('The source has no actual triangle bounds')
        b=c.execute('SELECT ARRAY[cube_ll_coord(bbox,1),cube_ll_coord(bbox,2),cube_ll_coord(bbox,3)]AS lo,ARRAY[cube_ur_coord(bbox,1),cube_ur_coord(bbox,2),cube_ur_coord(bbox,3)]AS hi FROM atlas_geom.meshes WHERE id=%s',(row['id'],)).fetchone()
        lo=(np.array(b['lo'])-distance_m).tolist();hi=(np.array(b['hi'])+distance_m).tolist()
        hidden=''if include_hidden else'AND NOT m.hidden AND NOT m.suppressed'
        result=c.execute('''SELECT m.id,m.part_id,m.name,m.system,m.hidden,m.suppressed,m.geometry_sha256,
                           m.bbox <-> s.bbox AS bbox_distance_lower_bound_m FROM atlas_geom.meshes m
                           JOIN atlas_geom.current_catalogues a ON a.catalogue_sha256=m.catalogue_sha256
                           CROSS JOIN atlas_geom.meshes s WHERE s.id=%s AND m.id<>s.id AND m.frame_id=s.frame_id
                           AND m.error IS NULL AND m.bbox&&cube(%s::float8[],%s::float8[]) '''+hidden+'''
                           ORDER BY bbox_distance_lower_bound_m,m.id LIMIT %s''',(row['id'],lo,hi,limit+1)).fetchall()
    return dict(source_id=row['id'],neighbors=result[:limit],truncated=len(result)>limit,broad_phase_only=True,anatomical_attachment_inferred=False)
