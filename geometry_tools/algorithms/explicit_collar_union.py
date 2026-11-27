"""Boolean union of isolated source collars, preserving their actual port edges.

Temporary computation caps close each input cell. The three distal caps are
identified by lineage and removed from the emitted carrier. Shared-root caps
must disappear in the union; otherwise the junction is rejected.
"""
import time
import numpy as np
import manifold3d as mf
from scipy.spatial import cKDTree
from .three_way_source_junction import triple_ring_parts
from .source_root_sections import recover_root_sections
from .tube_continuity import boundary_loops,continuity_audit
from geometry_tools.model import Mesh


def union_explicit_collars(mesh,*,source_layers=1):
 tick=time.monotonic();root,pieces=triple_ring_parts(mesh);removed=[];rings=[];cells=[];rootcaps=set();portcaps=set();channels={};values=[mesh.vertices];at=3
 _,physical_rep,physical=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True)
 if 'normals'in mesh.attributes:
  normal=mesh.attributes['normals'].astype(float)*mesh.attributes['normalSigns'][:,None];normal/=np.linalg.norm(normal,axis=1)[:,None];values.append(normal);channels['normals']=(at,at+3);at+=3
 for key in ['uvs','colors']:
  if key in mesh.attributes:
   q=mesh.attributes[key].astype(float);values.append(q);channels[key]=(at,at+q.shape[1]);at+=q.shape[1]
 properties=np.column_stack(values);condition_origin=mesh.vertices[root].mean(0);condition_scale=.001
 for branch,ids in enumerate(pieces):
  sections=recover_root_sections(mesh.vertices,mesh.faces[ids],root,maximum_layers=source_layers)
  if len(sections)!=source_layers:raise ValueError('Declared complete source collar depth is unavailable')
  sourceids=ids[np.concatenate([x['sourceFaceIds']for x in sections])];removed.extend(sourceids.tolist());ring=np.array(sections[-1]['rawIds']);f=mesh.faces[sourceids];ports=[physical_rep[loop]for loop in boundary_loops(physical[f])];props=properties.copy();faces=f.tolist();faceids=sourceids.tolist();source_ring_coordinates={tuple(x)for x in mesh.vertices[ring]};source_root_coordinates={tuple(x)for x in mesh.vertices[root]};actual_ring=None
  for loop in ports:
   coordinates={tuple(x)for x in mesh.vertices[loop]};isroot=coordinates==source_root_coordinates
   if not isroot and coordinates!=source_ring_coordinates:raise ValueError('Unexpected source collar boundary')
   if not isroot:actual_ring=loop
   capid=1000000+branch*2+int(not isroot);(rootcaps if isroot else portcaps).add(capid);centre=len(props);props=np.r_[props,(np.r_[condition_origin,properties[loop,3:].mean(0)]if isroot else props[loop].mean(0))[None]]
   for a,b in zip(loop,np.roll(loop,-1)):faces.append([centre,int(b),int(a)]);faceids.append(capid)
  faces=np.array(faces,int);used,inv=np.unique(faces,return_inverse=True);p=props[used].copy();p[:,:3]=(p[:,:3]-condition_origin)/condition_scale;ff=inv.reshape(-1,3);tri=p[ff,:3];origin=tri.mean((0,1));local=tri-origin;volume=np.einsum('ij,ij->i',local[:,0],np.cross(local[:,1],local[:,2])).sum()/6
  if volume<0:ff=ff[:,[0,2,1]]
  # Exact-coordinate seams are topology aliases, while each UV/normal property
  # corner remains distinct. Never approximate-weld or discard its properties.
  _,representatives,inverse=np.unique(p[:,:3],axis=0,return_index=True,return_inverse=True);canonical=representatives[inverse];seams=np.flatnonzero(canonical!=np.arange(len(p)))
  raw=mf.Mesh64(np.ascontiguousarray(p),np.ascontiguousarray(ff,dtype=np.uint64),merge_from_vert=seams.astype(np.uint64),merge_to_vert=canonical[seams].astype(np.uint64),face_id=np.array(faceids,dtype=np.uint64),tolerance=0);solid=mf.Manifold(raw)
  if solid.status()!=mf.Error.NoError or solid.is_empty():raise ValueError('Original source collar does not form a valid computation solid: '+str(solid.status()))
  cells.append(solid);assert actual_ring is not None;rings.append(actual_ring)
 union=mf.Manifold.batch_boolean(cells,mf.OpType.Add)
 if union.status()!=mf.Error.NoError:raise ValueError('Local collar union failed: '+str(union.status()))
 out=union.to_mesh64();props=np.asarray(out.vert_properties).copy();props[:,:3]=props[:,:3]*condition_scale+condition_origin;faces=np.asarray(out.tri_verts,int);faceids=np.array(out.face_id,int);parents=np.arange(len(props))
 def find(i):
  while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
  return i
 for a,b in zip(out.merge_from_vert,out.merge_to_vert):parents[find(int(a))]=find(int(b))
 canonical=np.array([find(i)for i in range(len(props))]);merge_error=float(np.linalg.norm(props[:,:3]-props[canonical,:3],axis=1).max(initial=0))
 if merge_error>1e-11:raise ValueError('CSG property-seam position discrepancy exceeds the declared numerical limit')
 props[:,:3]=props[canonical,:3]
 if set(faceids)&rootcaps:raise ValueError('A shared-root computation cap survived the union')
 keep=~np.isin(faceids,list(portcaps));carrier=faces[keep];lineage=faceids[keep]
 if not len(carrier)or np.any(lineage>=1000000):raise ValueError('Computation caps may not enter the emitted anatomy')
 removed=np.unique(removed);retained=np.setdiff1d(np.arange(len(mesh.faces)),removed);target=rings;assert len(target)==3
 # Match only controls from the explicitly declared source collars. This is not
 # whole-body proximity welding. The Float64 union may round an original input.
 sourcecontrols=np.unique(mesh.faces[removed]);tree=cKDTree(mesh.vertices[sourcecontrols]);dist,near=tree.query(props[:,:3]);mapping=np.full(len(props),-1,int);literal=dist<=1e-11
 for i in np.flatnonzero(literal):
  choices=sourcecontrols[tree.query_ball_point(props[i,:3],1e-11)]
  compatible=[j for j in choices if np.allclose(properties[j,3:],props[i,3:],atol=1e-7,rtol=1e-7)]
  if compatible:mapping[i]=compatible[0]
 vertices=mesh.vertices.tolist();attributes={k:v.tolist()for k,v in mesh.attributes.items()};newraw=[]
 for i in np.unique(carrier):
  if mapping[i]>=0:continue
  mapping[i]=len(vertices);vertices.append(props[i,:3].tolist());newraw.append(int(i))
  for key in attributes:
   if key=='normalSigns':value=1
   elif key=='normals':
    lo,hi=channels[key];n=props[i,lo:hi];n=n/np.linalg.norm(n);value=np.rint(n*32767).astype('<i2').tolist()
   elif key in channels:
    lo,hi=channels[key];value=props[i,lo:hi].astype(mesh.attributes[key].dtype).tolist()
   else:raise ValueError('Unsupported source attribute must not be silently invented: '+key)
   attributes[key].append(value)
 new=mapping[carrier];_,physical=np.unique(np.array(vertices),axis=0,return_inverse=True);oldedges={(int(physical[a]),int(physical[b]))for loop in target for a,b in zip(loop,np.roll(loop,-1))}
 for sign in [1,-1]:
  q=new if sign==1 else new[:,[0,2,1]]
  try:ports=boundary_loops(physical[q])
  except ValueError:continue
  edges={(int(a),int(b))for loop in ports for a,b in zip(loop,np.roll(loop,-1))}
  if edges==oldedges:new=q;break
 else:
  raise ValueError('Union did not retain all three literal, oppositely wound source ports: '+str(dict(expected=[x.tolist()for x in target],actual=[x.tolist()for x in ports],missingEdges=sorted(oldedges-edges),extraEdges=sorted(edges-oldedges),maximumInputControlRounding_m=float(dist.min(initial=0)),unionFaces=len(faces),retainedCarrierFaces=len(carrier))))
 attrs={k:np.array(v,dtype=mesh.attributes[k].dtype)for k,v in attributes.items()};candidate=Mesh(np.array(vertices),np.r_[mesh.faces[retained],new],mesh.frame,attrs,mesh.source);audit=continuity_audit(candidate.vertices,new,expected_mouths=3,expected_components=1,expected_euler=-1)
 if not audit['topological_continuity_pass']:raise ValueError('Local union is not one manifold three-way carrier')
 proof=dict(complete=True,method='Conditioned Float64 Boolean union of three explicitly isolated source collars',conditioningOrigin_m=condition_origin.tolist(),conditioningScale_m=condition_scale,sourceLayers=source_layers,sourceSharedRawIds=root.tolist(),sourcePortRawIds=[x.tolist()for x in rings],removedSourceFaceIds=removed.tolist(),retainedSourceFaceIds=retained.tolist(),newFacesStart=len(retained),newFaceCount=len(new),newVertexCount=len(newraw),newFaceSourceLineage=lineage.tolist(),computationRootCapsRetained=0,emittedLumenDisks=0,sourcePortMaximumMatchedRounding_m=float(dist[literal].max(initial=0)),nativePropertySeamMergePairs=len(out.merge_from_vert),nativePropertySeamMaximumRounding_m=merge_error,sourceControlsAttributesRetained=True,carrierAudit=audit,finiteVesselWallEstablished=False,fullSurfacePackingPatencyAndAnatomicalAcceptancePending=True,elapsedSeconds=time.monotonic()-tick)
 return candidate,proof
