"""Resolve declared overlapping finite assemblies without moving their union.

The caller chooses ownership and anatomical semantics. Shared charts inside
either assembly remain reciprocal; this operation does not infer fusion or
declare a coaptation surface to be an anatomical attachment.
"""
from collections import defaultdict
import numpy as np
from .finite_tissue_boolean import CommonPositionBank,subtract_occupied,barycentric
from .shared_triangle_edges import reconcile_boundaries,refine_surfaces

def subtract_assembly(positions,surfaces,cells,host_names,cutter_names,*,tolerance=2e-12,exact_source_replay=True):
 p=np.asarray(positions,float);host_names=set(host_names);cutter_names=set(cutter_names)
 if not host_names or not cutter_names or host_names&cutter_names:raise ValueError('Declare two distinct nonempty assembly groups')
 if not(host_names|cutter_names)<=set(cells):raise ValueError('Unknown declared tissue owner')
 owners=defaultdict(list)
 for name,uses in cells.items():
  for u in uses:owners[u['surface']].append((name,u))
 for sid,uses in owners.items():
  for group in [host_names,cutter_names]:
   contained=[name in group for name,u in uses]
   if any(contained)and not all(contained):raise ValueError('A declared assembly crosses an existing reciprocal interface: '+sid)
 host_surfaces={sid:fs for sid,fs in surfaces.items()if owners[sid][0][0]in host_names};host_uses={name:cells[name]for name in cells if name in host_names};host_faces={name:np.concatenate([host_surfaces[u['surface']]if u['side']==1 else host_surfaces[u['surface']][:,[0,2,1]]for u in uses])for name,uses in host_uses.items()};host_parents={name:[(u['surface'],i,u['side'])for u in uses for i in range(len(host_surfaces[u['surface']]))]for name,uses in host_uses.items()};host_rows={sid:dict(id=sid,kind='outer'if len(owners[sid])==1 else'shared',array=sid)for sid in host_surfaces};cutter_faces=[];cutter_owners=[];cutter_parent=[];internal={}
 for sid,fs in surfaces.items():
  uses=owners[sid]
  if uses[0][0]not in cutter_names:continue
  if len(uses)==2:internal[sid]=fs;continue
  name,u=uses[0];cutter_faces.extend(fs if u['side']==1 else fs[:,[0,2,1]]);cutter_owners.extend([name]*len(fs));cutter_parent.extend((sid,i,u['side'])for i in range(len(fs)))
 cf=np.array(cutter_faces,int);result=subtract_occupied(p,host_surfaces,host_rows,host_faces,host_parents,host_uses,p,cf,occupied_face_owners=np.array(cutter_owners,dtype=object),partition_free_outer=True,exact_source_replay=exact_source_replay,intersection_tolerance=tolerance,reconcile_shared_subdivisions=True);q=result['positions'];output=dict(result['surfaces']);uses={name:[dict(u)for u in ulist]for name,ulist in result['cells'].items()};origins={};covered=defaultdict(list)
 for sid,fs in output.items():
  row=next(r for r in result['surfaceRows']if r['id']==sid);pa=result['charts'][sid]['parentFaces'];family=row.get('parentSurface',sid)
  if family!='occupied-outer':origins[sid]=[(family,int(i),1)for i in pa];continue
  sign=row['parentNormalSign'];origins[sid]=[(cutter_parent[int(i)][0],cutter_parent[int(i)][1],sign*cutter_parent[int(i)][2])for i in pa]
  for face,index in zip(fs,pa):
   source_sid,ordinal,side=cutter_parent[int(index)];covered[(source_sid,ordinal)].append(face if sign*side==1 else face[[0,2,1]])
 for name in cutter_names:uses[name]=[dict(u)for u in cells[name]if u['surface']in internal]
 for u in result['vascularContacts']:uses[u['vascularOwner']].append(dict(surface=u['surface'],side=u['side'],neighbor=u['neighbor']))
 for u in result['vascularFreeSurfaces']:uses[u['vascularOwner']].append(dict(surface=u['surface'],side=u['side']))
 for sid,fs in internal.items():output[sid]=fs;origins[sid]=[(sid,i,1)for i in range(len(fs))]
 rest=set(cells)-host_names-cutter_names
 for name in rest:
  uses[name]=[dict(u)for u in cells[name]]
  for u in cells[name]:
   sid=u['surface'];output[sid]=surfaces[sid];origins[sid]=[(sid,i,1)for i in range(len(surfaces[sid]))]
 plans=defaultdict(set);deviation=0.
 for sid,fs in surfaces.items():
  if sid in internal or owners[sid][0][0]not in cutter_names:continue
  for ordinal,face in enumerate(fs):
   fragments=np.array(covered[(sid,ordinal)],int).reshape(-1,3)
   if not len(fragments):raise ValueError('Cutter source face disappeared')
   deviation=max(deviation,reconcile_boundaries(q,fragments,face[None],plans,tolerance=tolerance))
 bank=CommonPositionBank(q,0.);output,ancestry,changed=refine_surfaces(bank,output,plans);fields={};maximum_plane=0.;maximum_orientation_error=0.
 for sid,fs in output.items():
  refs=[origins[sid][int(i)]for i in ancestry[sid]];weights=[]
  for face,(source_sid,parent,side)in zip(fs,refs):
   tri=p[surfaces[source_sid][parent]];w=barycentric(q[face],tri);weights.append(w);maximum_plane=max(maximum_plane,float(np.max(np.linalg.norm(w@tri-q[face],axis=1))));fn=np.cross(*(q[face][1:]-q[face][0]));source_fn=np.cross(*(tri[1:]-tri[0]))*side
   if not fn@source_fn>0:maximum_orientation_error+=1
  fields[sid]=dict(source_surfaces=[v[0]for v in refs],source_faces=np.array([v[1]for v in refs],int),normal_signs=np.array([v[2]for v in refs],int),barycentrics=np.array(weights))
 if maximum_orientation_error or maximum_plane>tolerance*4:raise ValueError('Refined ownership fragment leaves its declared parent chart')
 return dict(positions=q,surfaces=output,cells=uses,charts=fields,evidence=dict(boolean=result['evidence'],edge_subdivision_count=len(plans),refined_faces=changed,maximum_collinear_deviation_m=deviation,maximum_parent_plane_distance_m=maximum_plane,source_controls_literal=np.array_equal(q[:len(p)],p),ownership_is_explicit_not_anatomical_inference=True))
