"""Rigid tube profiles with a harmonic displacement on intervening source cores.

Every original face stays owned. An optional physical identity vector must
come from source provenance; it never authorizes welding nearby geometry.
Raw alias offsets remain geometrically present and rotate with their profile.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

def incident_profile_core_closure(vertex_count, faces, profiles, moving_profile_ids):
 """Whole non-profile source patches incident to explicitly moved sections.

 Connectivity uses source control indices only. No spatial cutoff, coordinate
 welding or ownership inference is permitted. All other section controls stay
 outside the closure, including sections that bound the same branch patch.
 """
 n=int(vertex_count);f=np.asarray(faces,int);moving_ids=np.asarray(moving_profile_ids,int)
 if n<1 or f.ndim!=2 or f.shape[1]!=3 or np.any(f<0)or np.any(f>=n):raise ValueError('Valid source triangle controls required')
 if moving_ids.ndim!=1 or np.any(moving_ids<0)or np.any(moving_ids>=len(profiles)):raise ValueError('Explicit moving section indices required')
 owner=np.full(n,-1,int)
 for k,raw in enumerate(profiles):
  p=np.asarray(raw,int)
  if p.ndim!=1 or len(p)<3 or len(np.unique(p))!=len(p)or np.any(p<0)or np.any(p>=n)or np.any(owner[p]>=0):raise ValueError('Disjoint original source section controls required')
  owner[p]=k
 edges=np.unique(np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
 adj=coo_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n)).tocsr()
 free=np.flatnonzero(owner<0);mask=np.zeros(n,bool);evidence=[]
 if len(free):
  count,labels=connected_components(adj[free][:,free],directed=False)
  moving=np.isin(owner,moving_ids)&(owner>=0)
  for k in range(count):
   p=free[labels==k]
   if not adj[p][:,moving].nnz:continue
   mask[p]=True
   neighbours=adj[p].indices;bound=np.unique(owner[neighbours]);bound=bound[bound>=0]
   evidence.append(dict(coreComponent=int(k),controlCount=len(p),boundaryProfileIds=bound.tolist()))
 return dict(core_mask=mask,components=evidence,near_points_welded=False)

def profile_translation_weights(vertices,faces,profiles,*,physical_ids=None):
 """Exact linear harmonic-core influence for rigid section translations.

 Returns a raw-control by section matrix. Unanchored source components stay
 literal (zero rows). This preserves source aliases; it does not assign tissue
 ownership or certify the eventual geometry.
 """
 v=np.asarray(vertices,float);f=np.asarray(faces,int)
 ids=np.asarray(physical_ids,int)if physical_ids is not None else np.unique(v,axis=0,return_inverse=True)[1]
 if ids.shape!=(len(v),)or (ids<0).any():raise ValueError('One original identity per raw point required')
 _,inv=np.unique(ids,return_inverse=True);n=int(inv.max())+1
 field=np.zeros((n,len(profiles)));owned=np.full(n,-1,int)
 for k,p in enumerate(profiles):
  p=inv[np.asarray(p,int)]
  if len(set(p))!=len(p)or np.any(owned[p]>=0):raise ValueError('Profiles need disjoint physical controls')
  owned[p]=k;field[p,k]=1.
 ef=inv[f];edges=np.unique(np.sort(ef[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0);edges=edges[edges[:,0]!=edges[:,1]]
 adj=coo_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n)).tocsr()
 _,components=connected_components(adj,directed=False);fixed=owned>=0
 for ci in np.unique(components):
  mask=components==ci
  if not np.any(owned[mask]>=0):fixed[mask]=True
 boundary=np.flatnonzero(fixed);free=np.flatnonzero(~fixed);L=diags(np.asarray(adj.sum(1)).ravel())-adj
 if len(free):field[free]=spsolve(L[free][:,free].tocsc(),-(L[free][:,boundary]@field[boundary]))
 if not np.isfinite(field).all():raise ValueError('Invalid harmonic influence')
 return field[inv]

def transport_profile_cores(vertices,faces,profiles,target_centres,rotations,*,fixed_vertices=(),physical_ids=None,
                            prescribed_vertices=(),prescribed_positions=None,
                            core_method='harmonic',core_maximum_steps=32,core_wall_seconds=2.):
 v=np.asarray(vertices,float);f=np.asarray(faces,int);targets=np.asarray(target_centres,float);R=np.asarray(rotations,float)
 if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or targets.shape!=(len(profiles),3)or R.shape!=(len(profiles),3,3):raise ValueError('Source triangles and one proper rigid pose per profile required')
 if not np.isfinite(np.r_[v.ravel(),targets.ravel(),R.ravel()]).all()or not np.allclose(R.transpose(0,2,1)@R,np.eye(3),atol=1e-11,rtol=0)or np.any(np.linalg.det(R)<1-1e-11):raise ValueError('Finite inputs and proper profile rotations required')
 if physical_ids is None:_,representative,inv=np.unique(v,axis=0,return_index=True,return_inverse=True)
 else:
  ids=np.asarray(physical_ids,int)
  if ids.shape!=(len(v),)or (ids<0).any():raise ValueError('Explicit original source identity is required for each raw point')
  _,representative,inv=np.unique(ids,return_index=True,return_inverse=True)
 canonical=v[representative];n=len(canonical);fixed=np.unique(inv[np.asarray(fixed_vertices,int)]);fixed_mask=np.zeros(n,bool);fixed_mask[fixed]=True;boundary=set(map(int,fixed));d=np.zeros_like(canonical);profile_of=np.full(n,-1,int);profile_physical=[]
 for k,raw in enumerate(profiles):
  raw=np.asarray(raw,int);p=inv[raw]
  if len(p)<3 or len(set(p))!=len(p)or np.any(profile_of[p]>=0):raise ValueError('Profiles must have disjoint original physical point ownership')
  profile_of[p]=k;centre=canonical[p].mean(0);q=(canonical[p]-centre)@R[k].T+targets[k];delta=q-canonical[p]
  if fixed_mask[p].any()and not(np.array_equal(R[k],np.eye(3))and np.array_equal(targets[k],centre)):raise ValueError('A protected source point conflicts with the profile pose')
  d[p]=delta;boundary.update(map(int,p));profile_physical.append(p)
 prescribed=np.asarray(prescribed_vertices,int)
 if prescribed.size:
  positions=np.asarray(prescribed_positions,float)
  if prescribed.ndim!=1 or positions.shape!=(len(prescribed),3)or np.any(prescribed<0)or np.any(prescribed>=len(v))or not np.isfinite(positions).all():raise ValueError('Finite explicit source-point boundary controls required')
  seen={}
  for raw,target in zip(prescribed,positions):
   key=int(inv[raw]);value=target-v[raw]
   if key in seen and not np.allclose(value,seen[key],atol=1e-12,rtol=0):raise ValueError('Conflicting prescribed source aliases')
   if profile_of[key]>=0 or fixed_mask[key]:
    if not np.allclose(d[key],value,atol=1e-12,rtol=0):raise ValueError('Prescribed boundary conflicts with a rigid profile or fixed source point')
   seen[key]=value;d[key]=value;boundary.add(key)
 ef=inv[f];edges=np.unique(np.sort(ef[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0);edges=edges[edges[:,0]!=edges[:,1]];adj=coo_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n)).tocsr();degree=np.asarray(adj.sum(1)).ravel();_,components=connected_components(adj,directed=False);unanchored=[]
 for ci in np.unique(components):
  ids=np.flatnonzero(components==ci)
  if not any(int(i)in boundary for i in ids):boundary.update(map(int,ids));unanchored.append(int(ci))
 boundary=np.array(sorted(boundary),int);free=np.setdiff1d(np.arange(n),boundary);L=diags(degree)-adj;residual=0.
 if len(free):
  H=L[free][:,free].tocsc();rhs=-(L[free][:,boundary]@d[boundary]);d[free]=spsolve(H,rhs)
  if not np.isfinite(d).all():raise ValueError('Unresolved harmonic core system')
  residual=float(np.max(np.abs(H@d[free]-rhs)))
 arap_evidence=None
 if core_method=='arap':
  from .arap_core import arap_core
  posed,arap_evidence=arap_core(canonical,edges,boundary,canonical[boundary]+d[boundary],canonical+d,maximum_steps=core_maximum_steps,wall_seconds=core_wall_seconds)
  d=posed-canonical
 elif core_method!='harmonic':raise ValueError('Explicit harmonic or arap core method required')
 q=v+d[inv];raw_profile=profile_of[inv];profile_mask=raw_profile>=0
 # Retain each actual packed source offset instead of welding aliases to the
 # first representative. This makes every raw section congruent to its source.
 k=raw_profile[profile_mask];alias_offset=v[profile_mask]-canonical[inv[profile_mask]];q[profile_mask]+=np.einsum('nij,nj->ni',R[k]-np.eye(3),alias_offset);fixed_raw=np.isin(inv,fixed);q[fixed_raw]=v[fixed_raw];rawR=np.repeat(np.eye(3)[None],len(v),axis=0);rawR[profile_mask]=R[k]
 return dict(vertices=q,raw_profile_ids=raw_profile,profile_rotations=rawR,physical_ids=inv,physical_displacements=d,core_free_physical_ids=free,fixed_raw_point_ids=np.flatnonzero(fixed_raw),prescribed_raw_point_ids=prescribed,unanchored_components_retained_literal=unanchored,harmonic_equation_maximum_residual_m=residual,harmonic_residual_applies_to_initialization_only=core_method=='arap',core_method=core_method,arap_evidence=arap_evidence,source_faces_unchanged=True,near_points_welded=False)
