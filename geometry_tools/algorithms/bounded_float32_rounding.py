"""Resolve diagnosed storage crossings by bounded neighboring F32 words.

This is coordinate rounding, not general mesh fitting. Literal protected
source controls stay fixed, face orientation is preserved and the immutable
BVH guard rejects every new global crossing. Full final audits remain required.
"""
import itertools,time
import numpy as np
from .local_crossing_guard import LocalCrossingGuard


def round_diagnosed_crossings(positions,faces,defects,*,origin,protected=(),maximum_displacement_m=1e-7,seconds=20.,maximum_moves=32,maximum_pair_tests=2000000,avoid_new_literal_aliases=True,anchor_positions=None,maximum_displacement_by_vertex=None):
 start=time.monotonic();p=np.asarray(positions,float);f=np.asarray(faces,int);o=np.asarray(origin,float);fixed=set(map(int,protected));guard=LocalCrossingGuard(p,f,defects,maximum_displacement_m=maximum_displacement_m,maximum_pair_tests=maximum_pair_tests);reference=np.cross(p[f[:,1]]-p[f[:,0]],p[f[:,2]]-p[f[:,0]]);directions=sorted([q for q in itertools.product([-1,0,1],repeat=3)if any(q)],key=lambda q:sum(v*v for v in q));moves=[];reason=None
 occupied={tuple(q):set()for q in p}
 for i,q in enumerate(p):occupied[tuple(q)].add(i)
 anchor=p if anchor_positions is None else np.asarray(anchor_positions,float)
 if anchor.shape!=p.shape or not np.isfinite(anchor).all()or np.linalg.norm(p-anchor,axis=1).max()>maximum_displacement_m:raise ValueError('Valid immutable rounding anchor required')
 limits=np.full(len(p),maximum_displacement_m)if maximum_displacement_by_vertex is None else np.asarray(maximum_displacement_by_vertex,float)
 if limits.shape!=(len(p),)or not np.isfinite(limits).all()or np.any(limits<0)or np.any(limits>maximum_displacement_m):raise ValueError('Finite per-control bounds within the global limit required')
 epochs=np.zeros(len(p),np.int64);failed={};neighbors={};skipped=0
 for iteration in range(maximum_moves):
  if not guard.remaining:break
  if time.monotonic()-start>=seconds:reason='finite rounding deadline';break
  vertices=np.unique(f[np.array(sorted(guard.remaining)).ravel()]);accepted=False
  for vertex in vertices:
   vertex=int(vertex)
   if vertex in fixed:continue
   if vertex not in neighbors:
    pair_faces=guard.pairs_for([vertex]);incident0=guard.incident[guard.starts[vertex]:guard.starts[vertex+1]]
    neighbors[vertex]=np.unique(f[np.unique(np.r_[pair_faces.ravel(),incident0])])
   epoch=int(epochs[neighbors[vertex]].max(initial=0))
   if failed.get(vertex)==epoch:skipped+=1;continue
   old=guard.positions[vertex].copy();local=(old-o).astype('<f4');minus=np.nextafter(local,np.float32(-np.inf));plus=np.nextafter(local,np.float32(np.inf));incident=guard.incident[guard.starts[vertex]:guard.starts[vertex+1]]
   for direction in directions:
    if time.monotonic()-start>=seconds:reason='finite rounding deadline';break
    q=local.copy()
    for axis,sign in enumerate(direction):
     if sign:q[axis]=plus[axis]if sign>0 else minus[axis]
    value=q.astype(float)+o
    if avoid_new_literal_aliases and occupied.get(tuple(value),set())-{vertex}:continue
    if np.linalg.norm(value-anchor[vertex])>limits[vertex] or np.linalg.norm(value-p[vertex])>maximum_displacement_m:continue
    triangles=guard.positions[f[incident]].copy();triangles[f[incident]==vertex]=value;normal=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    if np.any(np.einsum('ij,ij->i',normal,reference[incident])<=0):continue
    proposal=guard.check_move([vertex],value[None])
    if proposal is not None and len(proposal.remaining)<len(guard.remaining):
     before=len(guard.remaining);guard.commit(proposal);epochs[vertex]=guard.revision;occupied[tuple(old)].discard(vertex);occupied.setdefault(tuple(value),set()).add(vertex);moves.append(dict(vertex=vertex,from_m=old.tolist(),to_m=value.tolist(),move_m=float(np.linalg.norm(value-old)),remaining_before=before,remaining_after=len(guard.remaining)));accepted=True;break
   if not accepted and not reason:failed[vertex]=epoch
   if accepted or reason:break
  if not accepted:
   if reason is None:reason='no single-word move improves the complete defect worklist'
   break
 return dict(positions=guard.positions,evidence=dict(passed=not guard.remaining,moves=moves,stop_reason=reason,maximum_allowed_displacement_m=maximum_displacement_m,maximum_actual_displacement_m=float(np.linalg.norm(guard.positions-p,axis=1).max()),maximum_total_anchor_displacement_m=float(np.linalg.norm(guard.positions-anchor,axis=1).max()),unchanged_failed_vertices_skipped=skipped,protected_controls_literal=bool(all(np.array_equal(guard.positions[i],p[i])for i in fixed)),elapsed_seconds=time.monotonic()-start,**guard.evidence()))
