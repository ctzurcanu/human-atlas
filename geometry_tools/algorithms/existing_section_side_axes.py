"""Local separation proposals oriented by existing profile offsets, without a global side flip."""
import numpy as np
def existing_section_side_axes(first,second,first_stations,second_stations,pairs,*,minimum_offset_m=1e-8):
 pairs=np.asarray(pairs,int)
 if pairs.ndim!=2 or pairs.shape[1]!=2 or minimum_offset_m<=0:raise ValueError('Declared source face pairs and positive side-resolution threshold required')
 def fields(mesh,station):
  station=np.asarray(station);assert station.shape==(len(mesh.vertices),) and np.all(station>=0)
  groups={g:np.flatnonzero(station==g) for g in np.unique(station)};centers={g:mesh.vertices[ids].mean(0) for g,ids in groups.items()};center=np.zeros_like(mesh.vertices);tangent=np.zeros_like(mesh.vertices)
  for group,ids in groups.items():
   _,singular,axes=np.linalg.svd(mesh.vertices[ids]-centers[group],full_matrices=False)
   if len(ids)>=3 and singular[1]>1e-12:axis=axes[-1]
   else:
    adjacent=np.unique(station[mesh.faces[np.any(np.isin(mesh.faces,ids),axis=1)]]);offsets=np.array([centers[g]-centers[group] for g in adjacent if g!=group])
    if not len(offsets):raise ValueError('Unresolved source chart has no adjacent course')
    _,_,course=np.linalg.svd(offsets,full_matrices=False);axis=course[0]
   center[ids]=centers[group];tangent[ids]=axis
  return center,tangent
 ca,ta=fields(first,first_stations);cb,tb=fields(second,second_stations);fa,fb=first.faces[pairs[:,0]],second.faces[pairs[:,1]];delta=ca[fa].mean(1)-cb[fb].mean(1);u=ta[fa[:,0]];v=tb[fb[:,0]];v*=np.where(np.einsum('ij,ij->i',u,v)<0,-1.,1.)[:,None]
 aligned=u+v;aligned/=np.maximum(np.linalg.norm(aligned,axis=1),1e-20)[:,None];radial=delta-np.einsum('ij,ij->i',delta,aligned)[:,None]*aligned;cross=np.cross(u,v);cross_len=np.linalg.norm(cross,axis=1);cross/=np.maximum(cross_len,1e-20)[:,None];side=np.einsum('ij,ij->i',cross,delta)
 use_cross=(cross_len>.1)&(abs(side)>minimum_offset_m);axes=np.where(use_cross[:,None],cross*np.where(side<0,-1.,1.)[:,None],radial);length=np.linalg.norm(axes,axis=1);resolved=length>minimum_offset_m;axes=np.divide(axes,length[:,None],out=np.zeros_like(axes),where=resolved[:,None])
 return axes,dict(complete=True,resolvedMask=resolved,sourcePairs=len(pairs),resolvedPairs=int(resolved.sum()),unresolvedPairs=int((~resolved).sum()),minimumExistingOffset_m=minimum_offset_m,globalPreferredAxisUsed=False,anatomicalSideInferred=False)
