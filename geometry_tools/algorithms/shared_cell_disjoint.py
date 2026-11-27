"""Independent disjointness certificate for exactly opposed shared facets.

Requires already audited, positive, closed oriented cells. Check every surface
contact and every nonshared face centroid in both directions. Unlike a volume
threshold, this rejects crossing, same-facing overlap and nested occupancy.
An uncertain/on-boundary winding number fails rather than authorizing a mesh.
"""
import numpy as np
from geometry_tools.contact import pair_contacts

def shared_cell_disjoint(first,second,*,max_pairs=250000,max_tests=2000000):
 import igl
 result,contacts=pair_contacts(first,second,0.,max_pairs=max_pairs,max_tests=max_tests)
 if not contacts['complete_within_pair_budget']or contacts['crossing_face_pairs']:
  return dict(passed=False,reason='crossing-or-incomplete-contacts',contacts=contacts)
 vertices,inverse=np.unique(np.r_[first.vertices,second.vertices],axis=0,return_inverse=True)
 a=inverse[:len(first.vertices)][first.faces];b=inverse[len(first.vertices):][second.faces]
 def oriented(t):
  q=tuple(map(int,t));return min(q,q[1:]+q[:1],q[2:]+q[:2])
 lookup={tuple(sorted(t)):i for i,t in enumerate(b)};ashared=set();bshared=set()
 for i,t in enumerate(a):
  j=lookup.get(tuple(sorted(t)))
  if j is None:continue
  if oriented(t)!=oriented(b[j,[0,2,1]]):return dict(passed=False,reason='same-facing-shared-facet')
  ashared.add(i);bshared.add(j)
 if not ashared:return dict(passed=False,reason='no-explicit-opposed-interface')
 for index in np.flatnonzero(result['relation']=='coplanar-overlap'):
  i,j=map(int,result['pairs'][index])
  if i not in ashared or lookup[tuple(sorted(a[i]))]!=j:
   return dict(passed=False,reason='nonidentical-coplanar-overlap',faces=[i,j])
 origin=(np.minimum(first.bounds[0],second.bounds[0])+np.maximum(first.bounds[1],second.bounds[1]))/2
 winding=[]
 for mesh,other,shared in [(first,second,ashared),(second,first,bshared)]:
  ids=np.asarray([i for i in range(len(mesh.faces))if i not in shared],int)
  if not len(ids):return dict(passed=False,reason='no-independent-exterior-witness')
  points=(mesh.triangles[ids]-origin).mean(1);w=igl.winding_number(other.vertices-origin,other.faces,points)
  row=dict(nonshared_faces_checked=len(ids),maximum_absolute_winding=float(abs(w).max()),inside=int(np.sum(abs(w)>.5)),uncertain=int(np.sum((abs(w)>1e-7)&(abs(w)<=.5))))
  winding.append(row)
  if not np.isfinite(w).all()or row['maximum_absolute_winding']>1e-7:
   return dict(passed=False,reason='contained-or-uncertain-nonshared-boundary',winding=winding,contacts=contacts)
 return dict(passed=True,method='complete-surface-contacts-and-bidirectional-generalized-winding',opposed_shared_faces=len(ashared),all_nonshared_faces_checked=True,winding=winding,contacts=contacts,anatomical_acceptance=False)
