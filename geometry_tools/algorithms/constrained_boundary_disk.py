"""Embedded computational disks on an existing nonplanar contour.

Candidate triangles come from a 3D Delaunay boundary triangulation; an interval
dynamic program minimizes their total area subject to real obstacle tests.
This uses the search-space idea described in CGAL's primary hole-filling manual
(https://doc.cgal.org/6.1/Polygon_mesh_processing/index.html), not CGAL's full
dihedral/area quality ordering. It neither constructs anatomical openings nor
accepts a model. Every caller must check its whole closed result independently.
"""
from collections import defaultdict
from itertools import combinations
import numpy as np
from scipy.spatial import Delaunay,QhullError,cKDTree
from scripts.ta98_triangle_intersections import sweep_hits


class BoundaryDiskError(ValueError):
    def __init__(self,message,details):super().__init__(message);self.details=details


def oriented_key(t):return min(tuple(map(int,np.roll(t,k)))for k in range(3))


def cancel_opposed_faces(faces):
    f=np.asarray(faces,int);owners=defaultdict(list)
    for i,t in enumerate(f):owners[tuple(sorted(map(int,t)))].append(i)
    remove=set();pairs=[]
    for ids in owners.values():
        if len(ids)==1:continue
        if len(ids)!=2 or oriented_key(f[ids[0]])!=oriented_key(f[ids[1]][[0,2,1]]):
            raise BoundaryDiskError('Coincident faces are not one opposed boundary pair',dict(faceIds=ids))
        remove.update(ids);pairs.append(ids)
    return np.delete(f,sorted(remove),axis=0),dict(opposedPairsCancelled=len(pairs),faceIds=pairs,coordinatesChanged=False)


def constrained_boundary_disk(positions,ring,wall_faces,other_caps=(),max_retries=6,forbidden_edges=()):
    v=np.asarray(positions,float);ring=np.asarray(ring,int)[::-1];wall=np.asarray(wall_faces,int);other=np.asarray(other_caps,int).reshape(-1,3)
    if len(set(map(int,ring)))!=len(ring)or len(ring)<3:raise ValueError('A simple literal contour is required')
    p=v[ring];center=p.mean(0);_,sing,axes=np.linalg.svd(p-center,full_matrices=False)
    try:
        delaunay=Delaunay(p);candidate_sets={tuple(sorted(int(q)for q in t if q!=skip))for t in delaunay.simplices for skip in t};dimension=3
    except QhullError:
        # This changes no coordinate; a planar contour has a 2D search space.
        delaunay=Delaunay((p-center)@axes[:2].T);candidate_sets={tuple(sorted(map(int,t)))for t in delaunay.simplices};dimension=2
    wall_tri=v[wall];wc=wall_tri.mean(1);radius=float(np.linalg.norm(wall_tri-wc[:,None],axis=2).max());tree=cKDTree(wc)
    literal={tuple(sorted(map(int,t))):i for i,t in enumerate(wall)};valid={};area={};probes=[];blocked=set()
    forbidden={tuple(sorted(map(int,e)))for e in forbidden_edges}
    def eligible(key):
        if key in blocked:return False
        if key in valid:return valid[key]
        ids=np.array(key,int);f=ring[ids];tri=v[f];a=float(np.linalg.norm(np.cross(tri[1]-tri[0],tri[2]-tri[0]))/2);area[key]=a
        if any(tuple(sorted((int(x),int(y))))in forbidden for x,y in zip(f,np.roll(f,-1))):valid[key]=False;return False
        if a==0:valid[key]=False;return False
        duplicate=literal.get(tuple(sorted(map(int,f))))
        if duplicate is not None and oriented_key(f)!=oriented_key(wall[duplicate][[0,2,1]]):valid[key]=False;return False
        c=tri.mean(0);r=float(np.linalg.norm(tri-c,axis=1).max());near=np.array(tree.query_ball_point(c,r+radius+1e-7),int)
        if duplicate is not None:near=near[near!=duplicate]
        hit=sweep_hits(tri[None],wall_tri[near],False,max_examples=1)if len(near)else dict(properOrCoplanarIntersections=0)
        caphit=sweep_hits(tri[None],v[other],False,max_examples=1)if len(other)else dict(properOrCoplanarIntersections=0)
        clean=hit['properOrCoplanarIntersections']==0 and caphit['properOrCoplanarIntersections']==0;valid[key]=clean
        probes.append(dict(triangle=list(key),wallCrossings=hit['properOrCoplanarIntersections'],priorCapCrossings=caphit['properOrCoplanarIntersections'],literalOpposedWallFace=duplicate is not None))
        return clean
    intervals=defaultdict(list)
    for key in candidate_sets:
        a,b,c=key;intervals[(a,c)].append(b)
    n=len(ring);retries=[];full_search=False;attempt=0
    while attempt<max_retries:
        costs=np.full((n,n),np.inf);split={}
        for i in range(n-1):costs[i,i+1]=0.
        for width in range(2,n):
            for i in range(n-width):
                j=i+width
                for k in intervals.get((i,j),()):
                    if not np.isfinite(costs[i,k]+costs[k,j]):continue
                    key=(i,k,j)
                    if not eligible(key):continue
                    value=costs[i,k]+costs[k,j]+area[key]
                    if value<costs[i,j]:costs[i,j]=value;split[(i,j)]=k
        if not np.isfinite(costs[0,n-1]):
            if not full_search:
                # A valid nonplanar patch need not use only Delaunay faces.
                # Preserve the same coordinates/guards; broaden the finite
                # combinatorial search, as in the documented hole-filling
                # fallback. Geometry tests remain lazy and cached.
                candidate_sets=set(combinations(range(n),3));intervals=defaultdict(list)
                for a,b,c in sorted(candidate_sets):intervals[(a,c)].append(b)
                full_search=True;continue
            raise BoundaryDiskError('No embedded disk in the complete boundary-triple search space',dict(dimension=dimension,controls=n,candidateTriangles=len(candidate_sets),trianglesTested=len(valid),probes=probes,retries=retries))
        triangles=[]
        def collect(i,j):
            if j-i<2:return
            k=split[(i,j)];triangles.append((i,k,j));collect(i,k);collect(k,j)
        collect(0,n-1);local=np.array(triangles,int);f=ring[local];hit=sweep_hits(v[f],v[f],True,max_examples=20)
        if hit['properOrCoplanarIntersections']==0:
            assert len(f)==n-2
            return f,dict(method='Obstacle-constrained Delaunay boundary disk',searchDimension=dimension,controls=n,candidateTriangles=len(candidate_sets),trianglesTested=len(valid),
                          completeBoundaryTripleFallback=full_search,forbiddenSourceEdges=len(forbidden),areaMetersSquared=float(costs[0,n-1]),selfIntersections=hit,probes=probes,retries=retries,newControls=0,anatomicalDisk=False)
        retries.append(dict(attempt=attempt,selfIntersections=hit))
        # Reject one crossed candidate and solve the same finite search again.
        a,b=hit['examples'][0];bad=max([tuple(map(int,local[a])),tuple(map(int,local[b]))],key=lambda key:area[key]);blocked.add(bad);attempt+=1
    raise BoundaryDiskError('Embedded disk retry bound exhausted',dict(controls=n,retries=retries,probes=probes))
