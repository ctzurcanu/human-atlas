"""Bounded diagonal flips at failed coat patches, preserving every control."""
from pathlib import Path
import argparse,json,time
from collections import defaultdict
import numpy as np
from scipy.spatial import ConvexHull
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

ap=argparse.ArgumentParser();ap.add_argument('carrier');ap.add_argument('--fields',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();source=Path(args.carrier).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
z=np.load(source/'surface.npz');v=z['vertices'];f=z['faces'].copy();d=json.loads((Path(args.fields)/'field-check.json').read_text())
bad=set(map(int,np.concatenate([np.asarray(x,int)for x in d['prisms']['invalid_faces_by_interval']])));edits=[]
def quality(ff):
    q=v[ff];area=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)
    return 2*np.sqrt(3)*area/np.sum((q-np.roll(q,1,axis=1))**2,axis=(1,2))
for iteration in range(24):
    adjacency=defaultdict(list)
    for i,tri in enumerate(f):
        for a,b in zip(tri,np.roll(tri,-1)):adjacency[tuple(sorted((int(a),int(b))))].append(i)
    best=None
    for first in sorted(bad):
        tri=f[first]
        for aa,bb in zip(tri,np.roll(tri,-1)):
            edge=tuple(sorted((int(aa),int(bb))));owners=adjacency[edge]
            if len(owners)!=2:continue
            second=next(i for i in owners if i!=first);a,b=int(aa),int(bb);c=int(next(i for i in tri if i not in edge));other=f[second]
            dd=[int(i)for i in other if i not in edge]
            if len(dd)!=1:continue
            d0=dd[0]
            if c==d0 or tuple(sorted((c,d0)))in adjacency:continue
            old=f[[first,second]];new=np.array([[c,a,d0],[c,d0,b]],int)
            # Select orientation by the unchanged patch's projected normal.
            normal=np.cross(v[old[:,1]]-v[old[:,0]],v[old[:,2]]-v[old[:,0]]).sum(0)
            nnew=np.cross(v[new[:,1]]-v[new[:,0]],v[new[:,2]]-v[new[:,0]])
            if np.all(nnew@normal<0):new=new[:,[0,2,1]]
            if np.any(np.cross(v[new[:,1]]-v[new[:,0]],v[new[:,2]]-v[new[:,0]])@normal<=0):continue
            ids=np.unique(np.r_[old.ravel(),new.ravel()]);points=v[ids];center=points.mean(0);_,_,basis=np.linalg.svd(points-center,full_matrices=False)
            uv=(points-center)@basis[:2].T;lookup={int(n):i for i,n in enumerate(ids)}
            try:hull=ConvexHull(uv)
            except Exception:continue
            if len(hull.vertices)!=4:continue
            areas=[]
            for faces in [old,new]:
                q=uv[[[lookup[int(i)]for i in t]for t in faces]]
                areas.append(float(abs((q[:,1,0]-q[:,0,0])*(q[:,2,1]-q[:,0,1])-(q[:,1,1]-q[:,0,1])*(q[:,2,0]-q[:,0,0])).sum()/2))
            if any(abs(area-hull.volume)>hull.volume*1e-8 for area in areas):continue
            # Both triangulations cover the same convex projected patch.
            # Its total height range bounds their surface separation.
            bound=float(np.ptp((points-center)@basis[2]))
            if bound>.00008:continue
            before=float(quality(old).min());after=float(quality(new).min());gain=after-before
            if gain<=max(.005,before*.1):continue
            score=gain/max(before,.001)
            if best is None or score>best[0]:best=(score,first,second,new,bound,before,after)
    if best is None:break
    _,first,second,new,bound,before,after=best;old=f[[first,second]].copy();f[first],f[second]=new
    edits.append(dict(faces=[first,second],before=old.tolist(),after=new.tolist(),projected_patch_surface_distance_bound_m=bound,
        minimum_quality_before=before,minimum_quality_after=after))
    bad.discard(first);bad.discard(second)
    if not bad:break
mesh=Mesh(v,f,'model:local-female-ta98:assembled-metres');proof=audit(mesh,require_closed=False,max_pairs=250000,max_tests=2000000)
atomic_npz(out/'surface.npz',vertices=v,faces=f,frame=np.array(mesh.frame))
metadata=json.loads((source/'surface.json').read_text());metadata.update(source_carrier_sha256=digest_file(source/'surface.npz'),
    carrier_face_repairs=edits,every_control_literal=True,audit=proof,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'surface.json',metadata)
assert proof['passed'],proof['failures']
print(json.dumps(dict(flips=len(edits),unresolved_seed_faces=sorted(bad),self_crossings=proof['self_intersections']['defect_pairs'],
    maximum_surface_change_bound_mm=max((e['projected_patch_surface_distance_bound_m']for e in edits),default=0)*1000,
    elapsed_seconds=time.monotonic()-start)),flush=True)
