"""Flip a local source diagonal without moving any anatomical control.

Only an oriented, convex, almost planar four-control patch can be changed.
This is a triangulation proposal, not anatomical or neighbor acceptance.
"""
import time
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics


def repair_diagonals(mesh,*,maximum_patch_thickness_m=.00025,wall_seconds=3.,maximum_tests=64,allow_partial=False,protected_face_ids=()):
 start=time.monotonic();initial=diagnostics(mesh,include_pairs=True);current=mesh;current_d=initial;changes=[];tested=0
 initial_pairs={tuple(map(int,p))for p in initial['source_face_pairs']}
 assert initial['source_face_pairs_complete']
 protected=set(map(int,protected_face_ids))
 if any(i<0 or i>=len(mesh.faces)for i in protected):raise ValueError('Protected source face outside mesh')
 while current_d['self_intersection_pairs'] and tested<maximum_tests:
  if time.monotonic()-start>wall_seconds:raise TimeoutError('Source diagonal repair deadline')
  f=current.faces;target=set(np.array(current_d['source_face_pairs']).ravel());edges={}
  for i,face in enumerate(f):
   for u,v in zip(face,np.roll(face,-1)):edges.setdefault(tuple(sorted((int(u),int(v)))),[]).append((i,int(u),int(v)))
  accepted=None
  for edge,inc in sorted(edges.items()):
   if len(inc)!=2 or not any(i in target for i,_,_ in inc):continue
   (i,u,v),(j,x,y)=inc
   if i in protected or j in protected:continue
   if (x,y)!=(v,u):continue
   c=next(int(k)for k in f[i]if k not in edge);d=next(int(k)for k in f[j]if k not in edge)
   if c==d or tuple(sorted((c,d)))in edges:continue
   ids=np.array([v,c,u,d]);q=current.vertices[ids];tri=current.triangles[[i,j]];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal=fn.sum(0);norm=np.linalg.norm(normal)
   if norm<1e-20:continue
   normal/=norm;thickness=float(np.ptp(q@normal))
   if thickness>maximum_patch_thickness_m:continue
   turns=np.cross(np.roll(q,-1,axis=0)-q,np.roll(q,-2,axis=0)-np.roll(q,-1,axis=0))@normal
   if np.any(turns<=1e-16):continue
   replacement=np.array([[c,d,v],[d,c,u]]);newtri=current.vertices[replacement];newfn=np.cross(newtri[:,1]-newtri[:,0],newtri[:,2]-newtri[:,0])
   if np.any(newfn@normal<=0)or np.linalg.norm(newfn,axis=1).min()<.15*np.linalg.norm(fn,axis=1).min():continue
   faces=f.copy();faces[[i,j]]=replacement;proposal=Mesh(current.vertices,faces,current.frame,current.attributes,current.source);tested+=1;diag=diagnostics(proposal,include_pairs=True)
   pairs={tuple(map(int,p))for p in diag['source_face_pairs']}
   if diag['source_face_pairs_complete']and pairs<=initial_pairs and diag['self_intersection_pairs']<current_d['self_intersection_pairs']and diag['nonmanifold_edges']==current_d['nonmanifold_edges']and diag['nonmanifold_vertices']==current_d['nonmanifold_vertices']:
    accepted=(proposal,diag,dict(faceIds=[i,j],oldFaces=f[[i,j]].tolist(),newFaces=replacement.tolist(),sourcePatchThickness_m=thickness));break
   if tested>=maximum_tests or time.monotonic()-start>wall_seconds:break
  if accepted is None:break
  current,current_d,proof=accepted;changes.append(proof)
 if time.monotonic()-start>wall_seconds:raise TimeoutError('Source diagonal repair deadline; partial output excluded')
 if current_d['self_intersection_pairs']and not(allow_partial and changes):raise ValueError(f'No complete local diagonal repair: {initial["self_intersection_pairs"]} -> {current_d["self_intersection_pairs"]}; tested {tested}')
 return current,dict(complete=True,initialCrossings=initial['self_intersection_pairs'],finalCrossings=current_d['self_intersection_pairs'],changedPatches=changes,tested=tested,protectedSourceFaceIds=sorted(protected),partialRepairAllowed=bool(allow_partial),allIntrinsicCrossingsResolved=current_d['self_intersection_pairs']==0,newIntrinsicCrossingPairs=0,coordinatesLiteral=True,attributesLiteral=True,maximumPatchThickness_m=maximum_patch_thickness_m,elapsedSeconds=time.monotonic()-start,neighborAcceptancePending=True)
