"""Replace one declared triple-ring pinch with a source-port three-way exterior.

Only certified source annular collars are removed. Original branch controls,
remote faces, open mouths and attributes remain unchanged. The local convex
carrier has three explicitly removed port caps; no lumen disk is emitted.
This is a local outer-envelope proposal, not a finite vascular wall.
"""
import itertools,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull,QhullError
from .source_root_sections import recover_root_sections
from .tube_continuity import boundary_loops,continuity_audit
from geometry_tools.model import Mesh


def triple_ring_parts(mesh):
 v,rep,inv=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True);f=inv[mesh.faces];es=np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1);ed,ei,count=np.unique(es,axis=0,return_inverse=True,return_counts=True);bad=ed[count>2]
 if len(bad)<3 or np.any(count[count>2]!=3):raise ValueError('One explicit three-sheet source ring is required')
 root=np.unique(bad);degree=np.bincount(bad.ravel(),minlength=len(v))
 if len(bad)!=len(root)or not np.all(degree[root]==2):raise ValueError('Nonmanifold edges are not one simple ring')
 graph=coo_matrix((np.ones(2*len(bad)),(np.r_[bad[:,0],bad[:,1]],np.r_[bad[:,1],bad[:,0]])),shape=(len(v),len(v))).tocsr()
 if connected_components(graph[root][:,root],directed=False)[0]!=1:raise ValueError('Multiple triple source rings are not inferred as one junction')
 order=np.argsort(ei);at=np.r_[0,np.cumsum(count)];pairs=[]
 for e in np.flatnonzero(count==2):pairs.append(order[at[e]:at[e+1]]//3)
 pairs=np.array(pairs,int);g=coo_matrix((np.ones(2*len(pairs)),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(f),len(f)));n,labels=connected_components(g,directed=False)
 # A named arterial mesh can contain other disconnected branches. Only the
 # three face groups incident to this certified ring belong to this surgery.
 # Unrelated groups must remain literal rather than preventing local repair.
 incident=np.unique(labels[np.any(np.isin(f,root),axis=1)])
 if len(incident)!=3:raise ValueError('Declared triple ring must separate exactly three incident source branches')
 pieces=[np.flatnonzero(labels==i)for i in incident]
 for edge in bad:
  owners=[int(labels[j])for j in np.flatnonzero(np.sum(np.isin(f,edge),axis=1)==2)]
  if sorted(owners)!=sorted(incident.tolist()):raise ValueError('Every triple-ring edge must belong to the same three source branches')
 return rep[root],pieces


def repair_three_way_junction(mesh,*,maximum_source_layers=3,wall_seconds=3.):
 start=time.monotonic();root,pieces=triple_ring_parts(mesh);sections=[]
 for ids in pieces:
  rows=recover_root_sections(mesh.vertices,mesh.faces[ids],root,maximum_layers=maximum_source_layers)
  if not rows:raise ValueError('Every source branch requires a certified complete collar')
  sections.append(rows)
 options=sorted(itertools.product(*(range(len(x))for x in sections)),key=lambda x:(sum(x),x));failures=[]
 for option in options:
  if time.monotonic()-start>wall_seconds:raise TimeoutError('Three-way source junction deadline')
  rings=[np.array(sections[i][k]['rawIds'])for i,k in enumerate(option)];removed=np.unique(np.concatenate([pieces[i][np.concatenate([np.array(row['sourceFaceIds'])for row in sections[i][:k+1]])]for i,k in enumerate(option)]));retained=np.setdiff1d(np.arange(len(mesh.faces)),removed)
  if not len(retained):continue
  raw=np.concatenate(rings)
  if len(np.unique(mesh.vertices[raw],axis=0))!=len(raw):continue
  try:hull=ConvexHull(mesh.vertices[raw])
  except QhullError:continue
  new=hull.simplices.copy();t=mesh.vertices[raw[new]];fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);flip=np.einsum('ij,ij->i',fn,hull.equations[:,:3])<0;new[flip]=new[flip][:,[0,2,1]];new=raw[new];cap=np.zeros(len(new),bool);cap_counts=[]
  for ring in rings:
   group=np.all(np.isin(new,ring),axis=1);cap_counts.append(int(group.sum()));cap|=group
   if not group.any()or not continuity_audit(mesh.vertices,new[group],expected_mouths=1,expected_components=1,expected_euler=1)['topological_continuity_pass']:break
  else:
   if all(cap_counts):pass
   else:continue
  if len(cap_counts)!=3 or not all(cap_counts):failures.append('Each source port must expose one actual disk in the hull');continue
  new=new[~cap]
  try:partial_loops=boundary_loops(new)
  except ValueError:continue
  if len(partial_loops)!=3:continue
  expansion={};valid=True
  for partial in partial_loops:
   match=[ring for ring in rings if set(partial)<=set(ring)]
   if len(match)!=1:valid=False;break
   ring=match[0];oriented=None
   for q in [ring,ring[::-1]]:
    where={int(v):j for j,v in enumerate(q)};advance=[(where[int(b)]-where[int(a)])%len(q)for a,b in zip(partial,np.roll(partial,-1))]
    if sum(advance)==len(q):oriented=q;break
   if oriented is None:valid=False;break
   where={int(v):j for j,v in enumerate(oriented)}
   for a,b in zip(partial,np.roll(partial,-1)):
    count=(where[int(b)]-where[int(a)])%len(oriented);path=[int(oriented[(where[int(a)]+j)%len(oriented)])for j in range(count+1)]
    if len(path)>2:expansion[int(a),int(b)]=path
  if not valid:continue
  expanded=[]
  for tri in new:
   hit=[(int(a),int(b),int(c))for a,b,c in zip(tri,np.roll(tri,-1),np.roll(tri,-2))if (int(a),int(b))in expansion]
   if not hit:expanded.append(tri.tolist());continue
   if len(hit)!=1:valid=False;break
   a,b,c=hit[0];path=expansion[a,b];expanded.extend([[u,v,c]for u,v in zip(path,path[1:])])
  if not valid:continue
  new=np.array(expanded,int)
  try:mouths=boundary_loops(mesh.faces[retained])
  except ValueError:continue
  port_loops=[]
  for ring in rings:
   matched=[x for x in mouths if set(x)==set(ring)]
   if len(matched)!=1:break
   port_loops.extend(matched)
  if len(port_loops)!=3:continue
  old_directed={(int(a),int(b))for loop in port_loops for a,b in zip(loop,np.roll(loop,-1))}
  for orientation in [1,-1]:
   proposed=new if orientation==1 else new[:,[0,2,1]]
   try:carrier_loops=boundary_loops(proposed)
   except ValueError:continue
   directed={(int(b),int(a))for loop in carrier_loops for a,b in zip(loop,np.roll(loop,-1))}
   if directed!=old_directed:continue
   cell=continuity_audit(mesh.vertices,proposed,expected_mouths=3,expected_components=1,expected_euler=-1)
   if not cell['topological_continuity_pass']:continue
   faces=np.r_[mesh.faces[retained],proposed];candidate=Mesh(mesh.vertices,faces,mesh.frame,mesh.attributes,mesh.source)
   source_counts=np.unique(np.sort(mesh.faces[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0,return_counts=True);open_edges=source_counts[0][source_counts[1]==1];actual=boundary_loops(faces);new_edges={tuple(sorted((int(a),int(b))))for loop in actual for a,b in zip(loop,np.roll(loop,-1))}
   if new_edges!=set(map(tuple,open_edges)):continue
   proof=dict(complete=True,sourceSharedRawIds=root.tolist(),sourceComponentFaces=[x.tolist()for x in pieces],sourcePortRawIds=[x.tolist()for x in rings],sourceLayers=[k+1 for k in option],removedSourceFaceIds=removed.tolist(),retainedSourceFaceIds=retained.tolist(),newFaceCount=len(proposed),newFacesStart=len(retained),sourcePortCapsExcluded=cap_counts,sourcePortContourEdgeExpansions={str(k):v for k,v in expansion.items()},lumenDisksCreated=0,sourceCoordinatesAttributesLiteral=True,externalSourceMouthEdgesLiteral=True,carrierAudit=cell,carrierHullVolumeUpperBound_m3=float(hull.volume),carrierBounds_m=[mesh.vertices[raw].min(0).tolist(),mesh.vertices[raw].max(0).tolist()],convexCarrier=not bool(expansion),insideOriginalSourcePortHull=True,finiteVesselWallEstablished=False,fullSurfaceNeighborAndAnatomicalAcceptancePending=True,elapsedSeconds=time.monotonic()-start)
   return candidate,proof
 raise ValueError('No source-port convex three-way exterior passes the declared boundaries: '+str(failures[:3]))
