"""Replace an explicitly certified source sleeve while retaining its branches.

The caller declares the two proper similarities. Source mouths remain literal
under those poses, a source-derived cross-section interpolates between them,
and an annular strip leaves the channel open. Packing and emitted-surface
validity are independent acceptance gates.
"""
import numpy as np
from numpy.polynomial import polynomial as poly
from .source_profile_axis import source_profile_graph
from .source_root_sections import source_annular_stem
from .port_hermite_course import port_hermite_course
from .tube_sections import unit,rotation_minimizing_frames
from .tube_rim_bridge import bridge_boundaries
from .tube_continuity import boundary_loops
from geometry_tools.model import Mesh


def _similarity(matrix):
 m=np.asarray(matrix,float)
 if m.shape!=(4,4)or not np.isfinite(m).all()or not np.array_equal(m[3],[0,0,0,1]):raise ValueError('Finite homogeneous similarity required')
 s=np.cbrt(np.linalg.det(m[:3,:3]));R=m[:3,:3]/s
 if not s>0 or not np.allclose(R.T@R,np.eye(3),atol=1e-11,rtol=0)or not np.isclose(np.linalg.det(R),1.,atol=1e-11):raise ValueError('Proper positive uniform similarity required')
 return m,s,R


def _sample_loop(points,count):
 p=np.asarray(points,float);closed=np.r_[p,p[:1]];distance=np.r_[0,np.cumsum(np.linalg.norm(np.diff(closed,axis=0),axis=1))]
 if np.any(np.diff(distance)<=0):raise ValueError('Source loop has repeated neighboring controls')
 if count==len(p):return p.copy()
 if count%len(p)==0:
  subdivisions=count//len(p);return np.array([a+(b-a)*j/subdivisions for a,b in zip(p,np.roll(p,-1,axis=0))for j in range(subdivisions)])
 return np.array([np.interp(np.arange(count)*distance[-1]/count,distance,closed[:,j])for j in range(3)]).T


def replace_source_stem(mesh,root_port,root_matrix,rest_matrix,*,spacing_m=.0015,
                        minimum_radius_ratio=1.25,maximum_length_ratio=3.):
 rootM,s0,R0=_similarity(root_matrix);restM,s1,R1=_similarity(rest_matrix)
 if not 0<spacing_m<=.003 or minimum_radius_ratio<1.:raise ValueError('Finite explicit spacing and non-pinching radius ratio required')
 p,f=mesh.vertices,mesh.faces;g=source_profile_graph(p,f,triangle_pairing=np.arange(len(f)//2*2).reshape(-1,2));root=g['boundary_raw'][root_port];stem=source_annular_stem(p,f,root,g['profiles']);cut=np.array(stem['cutRawIds'],int);removed=np.array(stem['sourceFaceIds'],int);retained=np.setdiff1d(np.arange(len(f)),removed);c0=p[root].mean(0);c1=p[cut].mean(0);old_axis=unit(c1-c0)
 def normal(raw):
  local=p[raw]-p[raw].mean(0);n=unit(np.cross(local,np.roll(local,-1,axis=0)).sum(0))
  return n if n@old_axis>0 else-n
 n0,n1=normal(root),normal(cut);u0=p[root[0]]-c0;u0=unit(u0-(u0@n0)*n0);u1=p[cut[0]]-c1;u1=unit(u1-(u1@n1)*n1);F0=np.column_stack((u0,np.cross(n0,u0),n0));F1=np.column_stack((u1,np.cross(n1,u1),n1));local0=(p[root]-c0)@F0*s0;local1=(p[cut]-c1)@F1*s1
 def signed_area(local):
  other=np.roll(local,-1,axis=0);return float((local[:,0]*other[:,1]-local[:,1]*other[:,0]).sum()/2)
 def area(local):return abs(signed_area(local))
 original_sections=[root]+[np.array(x['rawIds'],int)for x in stem['sourceSectionProfiles']];centres=np.array([p[raw].mean(0)for raw in original_sections]);distance=np.r_[0,np.cumsum(np.linalg.norm(np.diff(centres,axis=0),axis=1))]
 if np.any(np.diff(distance)<=1e-10):raise ValueError('Source sleeve contains coincident section centres')
 station=distance/distance[-1];normals=np.array([normal(raw)for raw in original_sections]);source_frames=rotation_minimizing_frames(centres,normals,initial_reference=u0);F0,F1=source_frames[[0,-1]];u0,u1=F0[:,0],F1[:,0];templates=[];source_areas=[];phase_choices=[]
 for k,(raw,centre,frame)in enumerate(zip(original_sections,centres,source_frames)):
  local=(p[raw]-centre)@frame;blend_scale=float(np.exp((1-station[k])*np.log(s0)+station[k]*np.log(s1)));local*=blend_scale;sampled=_sample_loop(local,len(root))
  if templates and np.sign(signed_area(sampled))!=np.sign(signed_area(templates[0])):sampled=sampled[::-1]
  choices=[];reference=templates[-1]if templates else sampled;reference_direction=unit(np.c_[reference[:,:2],np.zeros(len(root))])
  for phase in range(len(root)):
   q=np.roll(sampled,phase,axis=0);direction=unit(np.c_[q[:,:2],np.zeros(len(root))]);choices.append((float(np.linalg.norm(reference_direction-direction,axis=1).sum()),phase,q))
  _,phase,sampled=min(choices,key=lambda x:(x[0],x[1]));templates.append(sampled);source_areas.append(area(local));phase_choices.append(phase)
 templates=np.array(templates);A0,A1=source_areas[0],source_areas[-1]
 if min(A0,A1)<1e-12:raise ValueError('Source mouth has no stable projected flow area')
 radius=max(np.linalg.norm(x[:,:2],axis=1).max()for x in templates);start=c0@rootM[:3,:3].T+rootM[:3,3];end=c1@restM[:3,:3].T+restM[:3,3]
 coefficients,course=port_hermite_course(start,R0@n0,end,R1@n1,minimum_radius_m=minimum_radius_ratio*radius,maximum_length_ratio=maximum_length_ratio)
 if not course['accepted']:raise ValueError('Source-calibre stem course rejected: '+str(course))
 count=max(4,int(np.ceil(course['length_m']/spacing_m))+1)
 if count>160:raise ValueError('Source stem exceeds the finite station budget')
 ts=np.linspace(0,1,count);path=poly.polyval(ts,coefficients).T;tangents=unit(poly.polyval(ts,np.array([poly.polyder(x)for x in coefficients.T]).T).T);frames=rotation_minimizing_frames(path,tangents,initial_reference=R0@u0,final_reference=R1@u1)
 vertices=(p@restM[:3,:3].T+restM[:3,3]).tolist();root_xyz=p[root]@rootM[:3,:3].T+rootM[:3,3]
 for raw,value in zip(root,root_xyz):vertices[int(raw)]=value.tolist()
 faces=[];rings=[root];areas=[A0];radial=[];corrections=[]
 for k,t in enumerate(ts[1:-1],1):
  local=np.array([np.interp(t,station,templates[:,j,axis])for j in range(len(root))for axis in range(3)]).reshape(len(root),3);wanted=float(np.interp(t,station,np.sqrt(source_areas))**2);correction=np.sqrt(wanted/area(local));local[:,:2]*=correction
  if not .90<correction<1.15:raise ValueError('Interpolated source contour requires excessive area correction')
  q=path[k]+local@frames[k].T;ids=np.arange(len(vertices),len(vertices)+len(q));vertices.extend(q.tolist());previous=rings[-1]
  for i,j in zip(range(len(root)),np.roll(np.arange(len(root)),-1)):
   a,b,c,d=int(previous[i]),int(previous[j]),int(ids[j]),int(ids[i]);faces.extend([[a,b,c],[a,c,d]])
  rings.append(ids);areas.append(wanted);corrections.append(float(correction));radial.extend((q-path[k]).tolist())
 v=np.array(vertices);last,bridge=bridge_boundaries(v,rings[-1][::-1],cut);faces.extend(last.tolist());new_faces=np.array(faces,int);all_faces=np.r_[f[retained],new_faces]
 attrs={};old_count=len(p);added=len(v)-old_count
 for key,value in mesh.attributes.items():
  if key=='normals':
   signed=value.astype(float)*mesh.attributes['normalSigns'][:,None];signed=signed@R1.T;signed[root]=(value[root].astype(float)*mesh.attributes['normalSigns'][root,None])@R0.T;signed=unit(signed)
   tri=v[all_faces];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);sums=np.zeros_like(v)
   oldfn=np.cross(mesh.triangles[:,1]-mesh.triangles[:,0],mesh.triangles[:,2]-mesh.triangles[:,0]);oldn=value.astype(float)*mesh.attributes['normalSigns'][:,None];orientation=float(np.sign(np.median(np.einsum('ij,ij->i',oldfn,oldn[f].mean(1)))))
   if not orientation:raise ValueError('Source normal orientation is ambiguous')
   for i in range(3):np.add.at(sums,all_faces[:,i],fn*orientation)
   newn=unit(sums[old_count:]);attrs[key]=np.rint(np.r_[signed,newn]*32767).astype('<i2')
  elif key=='normalSigns':attrs[key]=np.ones(len(v),dtype='i1')
  else:attrs[key]=np.concatenate([value,np.repeat(value[root[:1]],added,axis=0)])
 candidate=Mesh(v,all_faces,mesh.frame,attrs,mesh.source)
 # Source endpoint and retained source-face geometry are asserted literally.
 np.testing.assert_array_equal(candidate.vertices[root],root_xyz);np.testing.assert_array_equal(candidate.triangles[:len(retained)],(p@restM[:3,:3].T+restM[:3,3])[f[retained]])
 _,rep,inv=np.unique(v,axis=0,return_index=True,return_inverse=True);loops=boundary_loops(inv[all_faces]);original_loops=boundary_loops(np.unique(p,axis=0,return_inverse=True)[1][f]);assert len(loops)==len(original_loops)
 proof=dict(sourceRootRawIds=root.tolist(),sourceCutRawIds=cut.tolist(),removedSourceFaceIds=removed.tolist(),retainedSourceFaceIds=retained.tolist(),retainedSourceFacesLiteralUnderDeclaredPose=True,sourceNeckAnnulus=stem['annularAudit'],course=course,stations=count,
  newFaces=len(new_faces),newVertices=added,constructedFacesStart=len(retained),finalBridge=bridge,rootMatrix=rootM.tolist(),restMatrix=restM.tolist(),templatePhase=phase,sourceMouthAreas_m2=[A0,A1],areaCorrectionRange=[min(corrections),max(corrections)],sectionAreas_m2=areas+[A1],
  lumenDisksCreated=0,finiteVesselWallEstablished=False,discreteAndNeighborAcceptancePending=True,anatomicalAcceptance=False,
  sourceSectionParameters=station.tolist(),sourceSectionProjectedAreas_m2=source_areas,sourceTemplatePhases=phase_choices)
 return candidate,proof
