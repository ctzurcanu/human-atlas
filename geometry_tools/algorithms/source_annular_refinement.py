"""Refine explicit source annuli along a declared common course displacement.

The original noncircular profile points and annulus end coordinates remain.
Only declared exporter quads are replaced; junction/core faces stay owned.
This constructs an outer source surface, not a finite wall or lumen proof.
"""
from collections import Counter
import numpy as np

def refine_annuli(source_vertices,target_vertices,faces,annuli,quad_pairs,
                  selected_strip_ids,displacement,*,maximum_station_m=.002,
                  attributes=None,maximum_new_vertices=100000,
                  proper_profile_poses=False,maximum_source_loft_area_error=.05,
                  profile_pose=None,skip_short_strips=False,
                  proper_source_profile_strips=(),station_limits_m=None,
                  profile_clearance=None,profile_pose_context=False):
 v=np.asarray(source_vertices,float);q=np.asarray(target_vertices,float);f=np.asarray(faces,int);pairs=np.asarray(quad_pairs,int);attrs={k:np.asarray(x)for k,x in (attributes or{}).items()}
 if v.shape!=q.shape or v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or not np.isfinite(np.r_[v.ravel(),q.ravel()]).all():raise ValueError('Finite corresponding source and posed vertices required')
 if not 0<maximum_station_m<=.01:raise ValueError('Specify a positive station step of at most10mm')
 if not np.isfinite(maximum_source_loft_area_error)or not 0<maximum_source_loft_area_error<=.05:raise ValueError('Source-profile area tolerance must be positive and at most5%')
 if profile_pose is not None and not proper_profile_poses:raise ValueError('An explicit station pose requires proper profile transport')
 if profile_clearance is not None and not proper_profile_poses:raise ValueError('Profile clearance requires rigid target profile poses')
 if pairs.ndim!=2 or pairs.shape[1]!=2 or len(np.unique(pairs))!=pairs.size:raise ValueError('Literal disjoint exporter triangle pairs required')
 if any(len(x)!=len(v)for x in attrs.values()):raise ValueError('Source attributes must address original vertices')
 physical=np.asarray(annuli['raw_to_physical']);rep=np.asarray(annuli['representative_raw_ids']);pf=physical[f];stripids=list(map(int,selected_strip_ids))
 if len(set(stripids))!=len(stripids):raise ValueError('A strip may be refined only once')
 source_pose_strips=set(map(int,proper_source_profile_strips))
 if not source_pose_strips<=set(stripids):raise ValueError('Proper source profile transport must address selected source strips')
 station_limits={}if station_limits_m is None else dict(station_limits_m)
 if not set(station_limits)<=set(stripids)or any(not np.isfinite(x)or x<=0 or x>maximum_station_m for x in station_limits.values()):raise ValueError('Per-strip station limits must refine selected source strips')
 source=list(v.copy());posed=list(q.copy());newattrs={k:list(x.copy())for k,x in attrs.items()};replace=set();newfaces=[];parent_faces=[];profiles=[];proof=[]
 pose_context=None;source_profile_context=None
 def add_point(a,b,t,extra,j):
  sv=(1-t)*v[rep[a]]+t*v[rep[b]]
  if source_profile_context is not None:
   lookup,profiles=source_profile_context;sv=profiles[j,lookup[int(a)]]
  if pose_context is None:tv=(1-t)*q[rep[a]]+t*q[rep[b]]+extra
  else:
   rotation,source_centre,target_centre=pose_context;tv=(sv-source_centre[j])@rotation[j].T+target_centre[j]+extra
  index=len(source);source.append(sv);posed.append(tv)
  if len(source)-len(v)>maximum_new_vertices:raise ValueError('Refinement exceeds the declared vertex budget')
  return index
 for stripid in stripids:
  strip=annuli['strips'][stripid];owned=set(map(int,strip['source_face_ids']));rings=[np.asarray(annuli['rings'][i],int)for i in strip['rings']];A,B=map(set,rings);ra,rb=map(lambda x:rep[x],rings);ca,cb=v[ra].mean(0),v[rb].mean(0);mapping={};quads=[]
  for ids in pairs:
   if not owned.intersection(map(int,ids)):continue
   if not set(map(int,ids))<=owned:raise ValueError('Partial exporter quad ownership is forbidden')
   pp=pf[ids];edges=Counter(tuple(sorted((int(a),int(b))))for tri in pp for a,b in zip(tri,np.roll(tri,-1)));rim=[e for e,n in edges.items()if n==1];long=[]
   for a,b in rim:
    if a in B and b in A:a,b=b,a
    if a in A and b in B:
     if a in mapping and mapping[a]!=b:raise ValueError('Source longitudinal edge correspondence is ambiguous')
     mapping[a]=b;long.append((a,b))
   if len(long)!=2:raise ValueError('Each source quad needs two actual longitudinal edges')
   quads.append((ids.copy(),long))
  if set(mapping)!=A or set(mapping.values())!=B or len(quads)!=len(A)or len(A)!=len(B):raise ValueError('An entire one-to-one source annulus is required')
  if replace.intersection(owned):raise ValueError('Refinement duplicates source face ownership')
  longlength=max(np.linalg.norm(v[rep[a]]-v[rep[b]])for a,b in mapping.items());steps=max(1,int(np.ceil(max(longlength,np.linalg.norm(cb-ca))/station_limits.get(stripid,maximum_station_m))))
  if steps==1:
   if skip_short_strips:continue
   raise ValueError('Selected strip already meets the station limit')
  replace.update(owned)
  ts=np.linspace(0,1,steps+1);centers=(1-ts[:,None])*ca+ts[:,None]*cb;source_profile_context=None
  if stripid in source_pose_strips:
   from scipy.spatial.transform import Rotation,Slerp
   ordered=list(map(int,rings[0]));a=v[rep[ordered]]-ca;b=v[[rep[mapping[x]]for x in ordered]]-cb;U,s,Vt=np.linalg.svd(a.T@b);R=Vt.T@np.diag([1.,1.,np.linalg.det(Vt.T@U.T)])@U.T;assert np.allclose(R.T@R,np.eye(3),atol=1e-12)and np.linalg.det(R)>0;aligned=b@R;shape=(1-ts[:,None,None])*a+ts[:,None,None]*aligned;rotations=Slerp([0.,1.],Rotation.from_matrix(np.array([np.eye(3),R])))(ts).as_matrix();source_profiles=np.einsum('nij,nkj->nki',rotations,shape)+centers[:,None];source_profile_context=({a:i for i,a in enumerate(ordered)},source_profiles)
  ds=np.asarray(displacement(centers),float);end=np.asarray(displacement(np.array([ca,cb])),float)
  if ds.shape!=centers.shape or end.shape!=(2,3)or not np.isfinite(ds).all():raise ValueError('Declared displacement must return finite metre vectors')
  bulge=ds-((1-ts[:,None])*end[0]+ts[:,None]*end[1]);bulge[[0,-1]]=0;point_indices={};pose_context=None
  if proper_profile_poses:
   from scipy.spatial.transform import Rotation,Slerp
   rotations=[];target_centres=[]
   for ids in (ra,rb):
    x,y=v[ids],q[ids];xc,yc=x.mean(0),y.mean(0);U,s,Vt=np.linalg.svd((x-xc).T@(y-yc));R=Vt.T@np.diag([1.,1.,np.linalg.det(Vt.T@U.T)])@U.T
    if np.max(np.linalg.norm((x-xc)@R.T+yc-y,axis=1))>5e-8:raise ValueError('An original annulus end is not a rigid source profile')
    rotations.append(R);target_centres.append(yc)
   if profile_pose is None:
    matrices=Slerp([0.,1.],Rotation.from_matrix(np.asarray(rotations)))(ts).as_matrix();target_centres=(1-ts[:,None])*target_centres[0]+ts[:,None]*target_centres[1]
   else:
    context=dict(strip_id=stripid,source_face_ids=sorted(owned))
    field_R,field_t=map(lambda x:np.asarray(x,float),profile_pose(centers,context)if profile_pose_context else profile_pose(centers))
    if field_R.shape!=(len(ts),3,3)or field_t.shape!=centers.shape or not np.isfinite(field_R).all()or not np.isfinite(field_t).all()or not np.allclose(field_R.transpose(0,2,1)@field_R,np.eye(3),atol=2e-12)or np.any(np.linalg.det(field_R)<1-2e-12):raise ValueError('Station pose field must contain finite proper rotations without scale')
    correction=np.array([rotations[0]@field_R[0].T,rotations[1]@field_R[-1].T]);matrices=Slerp([0.,1.],Rotation.from_matrix(correction))(ts).as_matrix()@field_R;field_centres=np.einsum('nij,nj->ni',field_R,centers)+field_t;delta=np.asarray(target_centres)-field_centres[[0,-1]];target_centres=field_centres+(1-ts[:,None])*delta[0]+ts[:,None]*delta[1];bulge[:]=0
   clearance_records={}
   if profile_clearance is not None:
    for j in range(1,steps):
     source_profile=(1-ts[j])*v[rep[rings[0]]]+ts[j]*v[[rep[mapping[int(a)]]for a in rings[0]]]
     if source_profile_context is not None:source_profile=source_profile_context[1][j]
     target_profile=(source_profile-centers[j])@matrices[j].T+target_centres[j]+bulge[j];shift,evidence=profile_clearance(target_profile,dict(strip_id=stripid,source_face_ids=sorted(owned),t=float(ts[j])));shift=np.asarray(shift,float)
     if shift.shape!=(3,)or not np.isfinite(shift).all()or np.linalg.norm(shift)>.0015:raise ValueError('Profile clearance exceeded the1.5mm proper translation bound')
     target_centres[j]+=shift;clearance_records[j]=dict(translation_m=shift.tolist(),evidence=evidence)
   pose_context=(matrices,centers,target_centres)
  for a,b in mapping.items():
   for j in range(1,steps):point_indices[a,j]=add_point(a,b,ts[j],bulge[j],j)
  def area(points):
   points=points-points.mean(0);return .5*np.linalg.norm(np.cross(points,np.roll(points,-1,axis=0)).sum(0))
  end_areas=[area(v[ra]),area(v[rb])]
  for j in range(1,steps):
   ids=[point_indices[int(a),j]for a in rings[0]];expected=(1-ts[j])*end_areas[0]+ts[j]*end_areas[1];relative=float(area(np.asarray(source)[ids])/expected-1)
   if abs(relative)>maximum_source_loft_area_error:raise ValueError('Intermediate source-profile area exceeds the declared end-area interpolation tolerance')
   profiles.append(dict(strip_id=stripid,t=float(ts[j]),raw_ids=ids,sourceLoftRelativeAreaError=relative,properTargetProfilePose=proper_profile_poses,properSourceProfilePose=stripid in source_pose_strips,profileClearance=clearance_records.get(j)if proper_profile_poses else None))
  for ids,long in quads:
   raw=f[ids].ravel();raw_by={int(p):int(raw[np.flatnonzero(physical[raw]==p)[0]])for p in np.unique(physical[raw])};(a,b),(c,d)=long
   # Fixed exporter winding is determined once from the actual source quad.
   initial=np.array([v[rep[a]],v[rep[c]],v[rep[d]]]);ref=np.cross(v[f[ids]][:,1]-v[f[ids]][:,0],v[f[ids]][:,2]-v[f[ids]][:,0]).sum(0);forward=np.dot(np.cross(initial[1]-initial[0],initial[2]-initial[0]),ref)>0
   def vertex(p,j):
    if j==0:return raw_by[p]
    if j==steps:return raw_by[mapping[p]]
    return point_indices[p,j]
   for j in range(steps):
    aa,cc,dd,bb=vertex(a,j),vertex(c,j),vertex(c,j+1),vertex(a,j+1);tri=[[aa,cc,dd],[aa,dd,bb]]
    if not forward:tri=[x[::-1]for x in tri]
    newfaces.extend(tri);parent_faces.extend([ids.tolist(),ids.tolist()])
   # Original corner attributes interpolate on their own chart. Geometry
   # aliases are shared, while original UV seams remain explicit.
   for x,y in long:
    for j in range(1,steps):
     k=point_indices[x,j]
     for name,values in attrs.items():
      while len(newattrs[name])<=k:newattrs[name].append(None)
      if name=='normalSigns':value=np.asarray(1,dtype=values.dtype)
      elif name=='normals':value=np.array([0,0,32767],dtype=values.dtype)
      else:value=((1-ts[j])*values[raw_by[x]]+ts[j]*values[raw_by[y]]).astype(values.dtype)
      previous=newattrs[name][k]
      if previous is not None and name not in ('normals','normalSigns')and not np.array_equal(previous,value):
       # A shared raw edge cannot silently merge distinct attribute charts.
       raise ValueError('Source longitudinal edge has an attribute seam; use explicit corner charts')
      newattrs[name][k]=value
  proof.append(dict(strip_id=stripid,source_face_ids=sorted(owned),stations=steps+1,source_profile_ids=list(map(int,strip['rings'])),original_ends_unchanged=True))
 source=np.asarray(source);posed=np.asarray(posed);kept=np.asarray([i for i in range(len(f))if i not in replace],int);outfaces=np.vstack([f[kept],np.asarray(newfaces,int).reshape(-1,3)]);resultattrs={k:np.asarray(values,dtype=attrs[k].dtype)for k,values in newattrs.items()}
 # Actual resulting incident areas define smooth shared normals. Keep all
 # unmodified source normal words and signs outside the refined charts.
 if 'normals'in resultattrs:
  signs=resultattrs.get('normalSigns',np.ones(len(source),dtype='i1'));oldtri=q[f];oldcross=np.cross(oldtri[:,1]-oldtri[:,0],oldtri[:,2]-oldtri[:,0]);original_effective=attrs['normals'].astype(float)*attrs.get('normalSigns',np.ones(len(v)))[:,None];orientation=1 if np.median(np.einsum('ij,ij->i',oldcross,original_effective[f].mean(1)))>=0 else-1
  _,inv=np.unique(posed,axis=0,return_inverse=True);sums=np.zeros((inv.max()+1,3));t=posed[outfaces];cross=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])*orientation;np.add.at(sums,inv[outfaces].ravel(),np.repeat(cross,3,axis=0));length=np.linalg.norm(sums,axis=1);affected=np.unique(inv[np.unique(np.asarray(newfaces).ravel())]);mask=np.isin(inv,affected);assert np.all(length[inv[mask]]>0);norm=sums[inv[mask]]/length[inv[mask],None];resultattrs['normals'][mask]=np.rint(np.clip(norm*signs[mask,None],-1,1)*32767).astype(attrs['normals'].dtype)
 return dict(source_vertices=source,vertices=posed,faces=outfaces,attributes=resultattrs,
             kept_source_face_ids=kept,replaced_source_face_ids=np.asarray(sorted(replace)),
             new_face_parent_pairs=np.asarray(parent_faces,int).reshape(-1,2),new_profiles=profiles,strips=proof,
             original_vertex_count=len(v),new_vertices=len(source)-len(v),new_triangles=len(outfaces)-len(f),
             source_faces_partitioned_once=True,finite_wall_established=False,anatomical_acceptance=False)
