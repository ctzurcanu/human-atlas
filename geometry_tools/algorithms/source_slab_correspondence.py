"""Literal source-face correspondence through an explicitly clipped tube slab.

Connectivity outside the two planes cannot join otherwise separate strands.
True forks inside the slab remain explicit many-to-one/one-to-many components.
No point proximity, semantic ordering or nearest-centre assignment is used.
"""
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .scalar_surface_clip import clip_scalar_surface


def source_slab_correspondence(vertices,faces,low_y,high_y,first_face_sets,second_face_sets):
 v,f=np.asarray(vertices,float),np.asarray(faces,int)
 if not np.isfinite([low_y,high_y]).all()or low_y>=high_y:raise ValueError('Two ordered literal Y planes required')
 active=np.flatnonzero((v[f,1].max(1)>=low_y)&(v[f,1].min(1)<=high_y))
 if not len(active):raise ValueError('No actual source triangles meet the slab')
 cut=clip_scalar_surface(v,f[active],v[:,1],low_y);v1,f1=cut['positions'],cut['positive_faces'];parent1=active[cut['positive_parents']]
 if not len(f1):raise ValueError('No source surface above first plane')
 cut=clip_scalar_surface(v1,f1,v1[:,1],high_y);v2,f2=cut['positions'],cut['negative_faces'];parents=parent1[cut['negative_parents']]
 if not len(f2):raise ValueError('No source surface within slab')
 _,inv=np.unique(v2,axis=0,return_inverse=True);physical=inv[f2];owners=defaultdict(list)
 for fi,tri in enumerate(physical):
  for a,b in zip(tri,np.roll(tri,-1)):owners[tuple(sorted((int(a),int(b))))].append(fi)
 if any(len(values)>2 for values in owners.values()):raise ValueError('Nonmanifold literal source slab has unresolved ownership')
 edges=np.array([values for values in owners.values()if len(values)==2],int).reshape(-1,2);_,labels=connected_components(coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(f2),len(f2))).tocsr(),directed=False)
 def component_of(face_sets):
  answer=[]
  for values in face_sets:
   values=np.asarray(values,int)
   if not len(values)or values.min()<0 or values.max()>=len(f):raise ValueError('Actual nonempty source-face ownership required at each end')
   hits=np.flatnonzero(np.isin(parents,values));found=np.unique(labels[hits])
   if len(found)!=1:raise ValueError('Section does not map to one literal slab component')
   answer.append(int(found[0]))
  return answer
 first,second=component_of(first_face_sets),component_of(second_face_sets);links=[];components=[]
 for component in sorted(set(first)|set(second)):
  aa=[i for i,value in enumerate(first)if value==component];bb=[i for i,value in enumerate(second)if value==component];components.append(dict(component=component,firstProfiles=aa,secondProfiles=bb,sourceFaces=np.unique(parents[labels==component]).tolist(),explicitForkOrMerge=len(aa)>1 or len(bb)>1))
  if aa and bb:links.extend([i,j,1./(len(aa)*len(bb))]for i in aa for j in bb)
 return dict(links=links,components=components,literalPlaneInterval_m=[float(low_y),float(high_y)],sourceFacesPreserved=True,connectivityOutsideSlabExcluded=True,nearestCentreAssignmentUsed=False,anatomicalAcceptance=False)
