"""Partition existing shared cells with one common set of plane intersections.

Each new cap is built once and used by the two adjoining tissues with opposite
ownership. Inner contours are retained as holes, never filled with a disk.
"""
from collections import defaultdict,Counter
import numpy as np
import manifold3d as md
from .source_face_complement import conform_surface_boundary_edges


def partition(positions,surfaces,cells,bands,*,origin,axis=1):
    origin=np.asarray(origin,float);axis=int(axis)
    if axis not in [0,1,2]:raise ValueError('A coordinate axis is required')
    for a,b in zip(bands,bands[1:]):
        if abs(a[2]-b[1])>1e-12:raise ValueError('Bands must share literal boundaries')
    points=[];lookup={}
    def add(p):
        local=(np.asarray(p)-origin).astype('<f4');key=local.tobytes()
        if key not in lookup:lookup[key]=len(points);points.append(local.astype(float)+origin)
        return lookup[key]
    original=np.array([add(p)for p in positions],int)
    faces={n:original[np.asarray(f,int)]for n,f in surfaces.items()}
    planes=[float(np.float32(bands[i][2]-origin[axis]))+origin[axis]for i in range(len(bands)-1)]
    limits=[-np.inf,*planes,np.inf];edge_points={};source_edges={};initial_count=len(points)
    def crossing(a,b,plane):
        # Reuse the original edge for every plane. Re-intersecting a rounded
        # point from a previous parallel plane can otherwise produce two
        # different Float32 positions on the very same shared cut.
        if a<initial_count and b<initial_count:edge=tuple(sorted((a,b)))
        elif a in source_edges and b in source_edges[a]:edge=source_edges[a]
        elif b in source_edges and a in source_edges[b]:edge=source_edges[b]
        elif source_edges.get(a)==source_edges.get(b)and a in source_edges:edge=source_edges[a]
        else:raise ValueError('A parallel cut left its original triangle edge')
        key=(*edge,plane)
        if key not in edge_points:
            p,q=np.array(points[edge[0]]),np.array(points[edge[1]]);t=(plane-p[axis])/(q[axis]-p[axis]);r=p+t*(q-p);r[axis]=plane
            edge_points[key]=add(r)
            if edge_points[key]>=initial_count:source_edges[edge_points[key]]=edge
        return edge_points[key]
    def clip(poly,plane,above):
        if not np.isfinite(plane):return poly
        result=[]
        for a,b in zip(poly,poly[1:]+poly[:1]):
            da=points[a][axis]-plane;db=points[b][axis]-plane
            ia=da>=0 if above else da<=0;ib=db>=0 if above else db<=0
            if ia:result.append(a)
            if ia!=ib:result.append(crossing(a,b,plane))
        return list(dict.fromkeys(result))
    out={};uses={n+'@'+band[0]:[]for n in cells for band in bands};owners=defaultdict(list)
    for n,uu in cells.items():
        for u in uu:owners[u['surface']].append((n,u))
    source_points=np.asarray(points)
    for sid,ff in faces.items():
        # Most facets lie wholly in one slab. Classify them in bulk; clip
        # only genuine boundary facets. This avoids visiting every triangle
        # in Python for every segment of a long anatomical assembly.
        ordinates=source_points[ff,axis]
        minimum=ordinates.min(axis=1);maximum=ordinates.max(axis=1)
        triangles=source_points[ff]
        valid=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)>0
        for k,(band,_,_)in enumerate(bands):
            inside=valid&(minimum>=limits[k])&(maximum<=limits[k+1])
            intersects=valid&(maximum>=limits[k])&(minimum<=limits[k+1])&~inside
            emitted=[]
            for tri in ff[intersects]:
                poly=clip(list(map(int,tri)),limits[k],True)
                if len(poly)<3:continue
                poly=clip(poly,limits[k+1],False)
                if len(poly)<3:continue
                for j in range(1,len(poly)-1):emitted.append([poly[0],poly[j],poly[j+1]])
            if not emitted and not inside.any():continue
            clipped=np.asarray(emitted,int).reshape(-1,3)
            if len(clipped):
                q=np.asarray([points[i]for i in clipped.ravel()]).reshape(-1,3,3)
                clipped=clipped[np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0]
            f=np.concatenate([ff[inside],clipped])
            if not len(f):continue
            key=sid+'@'+band;out[key]=f
            for name,u in owners[sid]:
                neighbor=u.get('neighbor');uses[name+'@'+band].append(dict(surface=key,side=u['side'],neighbor=neighbor+'@'+band if neighbor else None))
    # Find the actual missing cap contours of the upper slab, then reverse
    # its border orientation to close that slab. The lower slab uses the
    # identical chart with the opposite sign.
    axes=[i for i in range(3)if i!=axis];cap_count=0
    v=np.asarray(points)
    for k,plane in enumerate(planes):
        lower=bands[k][0];upper=bands[k+1][0]
        for name in cells:
            uu=uses[name+'@'+upper];fs=[]
            for u in uu:
                f=out[u['surface']];fs.append(f if u['side']==1 else f[:,[0,2,1]])
            f=np.concatenate(fs);ee=np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]
            counts=Counter(tuple(sorted(e))for e in ee)
            edges=[(int(b),int(a))for a,b in ee if counts[tuple(sorted((a,b)))]==1 and points[a][axis]==plane and points[b][axis]==plane]
            if not edges:continue
            nxt={a:b for a,b in edges}
            if len(nxt)!=len(edges)or set(nxt)!=set(nxt.values()):raise ValueError('Cut contour is branched or incomplete: '+name)
            todo=set(nxt);loops=[]
            while todo:
                a=min(todo);loop=[];i=a
                while i in todo:todo.remove(i);loop.append(i);i=nxt[i]
                if i!=a:raise ValueError('Cut contour does not close')
                loops.append(np.array(loop,int))
            polys=[v[loop][:,axes]*1000 for loop in loops]
            if sum(float(np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)))for p in polys)<0:
                loops=[loop[::-1]for loop in loops];polys=[p[::-1]for p in polys]
            ids=np.concatenate(loops);cap=ids[np.asarray(md.triangulate(polys,epsilon=0.,allow_convex=False),int)]
            q=v[cap];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);cap[normal[:,axis]>0]=cap[normal[:,axis]>0][:,[0,2,1]]
            key='cut:'+str(k)+':'+name;out[key]=cap
            uses[name+'@'+upper].append(dict(surface=key,side=1,neighbor=name+'@'+lower))
            uses[name+'@'+lower].append(dict(surface=key,side=-1,neighbor=name+'@'+upper));cap_count+=1
    v=np.asarray(points)
    # A triangulator may omit a collinear boundary control. Insert existing
    # controls back in its boundary; no cap area or anatomy is inferred here.
    out,proof=conform_surface_boundary_edges(v,out,uses,tolerance=1e-8)
    return dict(positions=v,surfaces=out,cells=uses,evidence=dict(
        shared_cut_caps=cap_count,common_edge_intersections=len(edge_points),
        encoded_position_bank=True,cap_holes_preserved=True,boundary_conformity=proof,
        maximum_input_encoding_error_m=float(np.linalg.norm(v[original]-positions,axis=1).max())))
