"""Move a source vessel by complete section poses with explicit mouth anchors.

Source incidence defines every co-moving section and core. No cylinder is
substituted, no points are proximity-welded, and no surface validity is inferred
from the graph solve. Fixed and moving mouths are caller-declared interfaces.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
from scipy.spatial.transform import Rotation
from geometry_tools.model import Mesh
from .source_profile_axis import source_profile_graph,complete_source_sections
from .profile_network_arap import fit_profile_network
from .profile_graph_transport import transport_linked_profiles


def pose_anchored_sections(mesh,matrix,*,fixed_ports,moving_ports,
                           unanchored_components,wall_seconds=3.,iterations=20,core_method='rigid',recover_necks=False):
 m=np.asarray(matrix,float);fixed_ports=set(map(int,fixed_ports));moving_ports=set(map(int,moving_ports))
 if m.shape!=(4,4)or not np.isfinite(m).all()or not np.array_equal(m[3],[0,0,0,1]):raise ValueError('Explicit homogeneous proper similarity required')
 scale=np.cbrt(np.linalg.det(m[:3,:3]));R=m[:3,:3]/scale
 if not scale>0 or not np.allclose(R.T@R,np.eye(3),atol=1e-11,rtol=0)or not np.isclose(np.linalg.det(R),1.,atol=1e-11):raise ValueError('Positive uniform scale and proper rotation required')
 if not fixed_ports or not moving_ports or fixed_ports&moving_ports:raise ValueError('Disjoint explicit fixed and moving mouths required')
 if unanchored_components not in ['fixed','moving']:raise ValueError('Explicit orphan-component pose policy required')
 if core_method not in ['rigid','arap']:raise ValueError('Explicit rigid or ARAP junction-core method required')
 p=mesh.vertices;f=mesh.faces;physical,inv=np.unique(p,axis=0,return_inverse=True)
 guide=source_profile_graph(p,f,triangle_pairing=np.arange(len(f)//2*2).reshape(-1,2),ring_size=6)
 recovered=[]
 if recover_necks:
  from .source_root_sections import recover_root_sections
  known=[set(inv[raw])for raw in guide['profiles']+guide['boundary_raw']]
  for port_index,port_raw in enumerate(guide['boundary_raw']):
   for row in recover_root_sections(p,f,port_raw):
    ids=set(inv[row['rawIds']]);overlap=[old for old in known if old&ids]
    if overlap:
     if any(ids!=old for old in overlap):break
     continue
    guide['profiles'].append(np.array(row['rawIds'],int));known.append(ids);recovered.append(dict(row,port=port_index))
 ports=guide['boundary_raw'];all_ports=set(range(len(ports)))
 if fixed_ports|moving_ports!=all_ports:raise ValueError('Every source mouth requires an explicit pose policy')
 owner=np.full(len(physical),-1,int);groups=complete_source_sections(guide,inv)
 for i,raw in enumerate(groups):owner[inv[raw]]=i
 ef=inv[f];edges=np.unique(np.sort(ef[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
 free=np.flatnonzero(owner<0);local=np.full(len(owner),-1);local[free]=np.arange(len(free));live=edges[np.all(owner[edges]<0,axis=1)]
 graph=coo_matrix((np.ones(2*len(live)),(np.r_[local[live[:,0]],local[live[:,1]]],np.r_[local[live[:,1]],local[live[:,0]]])),shape=(len(free),len(free))).tocsr()
 count,label=connected_components(graph,directed=False)
 for i in range(count):owner[free[label==i]]=len(groups)+i
 raw_owner=owner[inv];n=int(owner.max())+1;centres=np.array([physical[owner==i].mean(0)for i in range(n)])
 links=np.unique(np.sort(owner[edges],axis=1),axis=0);links=links[links[:,0]!=links[:,1]]
 port_nodes=[]
 for raw in ports:
  nodes=np.unique(raw_owner[raw])
  if len(nodes)!=1:raise ValueError('Source mouth did not close into one complete rigid section')
  port_nodes.append(int(nodes[0]))
 fixed_nodes={port_nodes[i]for i in fixed_ports};moving_nodes={port_nodes[i]for i in moving_ports}
 if fixed_nodes&moving_nodes:raise ValueError('Intersecting source sections have incompatible interface poses')
 aa,zz=links.T;adj=coo_matrix((np.ones(2*len(links)),(np.r_[aa,zz],np.r_[zz,aa])),shape=(n,n)).tocsr();components,labels=connected_components(adj,directed=False)
 anchors=fixed_nodes|moving_nodes;orphans=[]
 for i in range(components):
  nodes=set(np.flatnonzero(labels==i))
  if not nodes&anchors:
   orphans.append(dict(nodes=sorted(map(int,nodes)),policy=unanchored_components))
   (moving_nodes if unanchored_components=='moving'else fixed_nodes).update(nodes)
 anchors=np.array(sorted(fixed_nodes|moving_nodes),int);free_nodes=np.setdiff1d(np.arange(n),anchors);weights=np.zeros(n);weights[list(moving_nodes)]=1.
 lap=diags(np.asarray(adj.sum(axis=1)).ravel())-adj
 if len(free_nodes):weights[free_nodes]=spsolve(lap[free_nodes][:,free_nodes].tocsc(),-lap[free_nodes][:,anchors]@weights[anchors])
 if np.any(weights< -1e-10)or np.any(weights>1+1e-10):raise ValueError('Invalid harmonic section pose weights')
 weights=np.clip(weights,0,1);rotvec=Rotation.from_matrix(R).as_rotvec();reference=Rotation.from_rotvec(weights[:,None]*rotvec).as_matrix();scales=np.exp(weights*np.log(scale))
 uniform=centres@m[:3,:3].T+m[:3,3];targets=centres+weights[:,None]*(uniform-centres)
 fit=fit_profile_network(centres,targets,links,fixed=anchors,reference_rotations=reference,
  target_weights=np.full(n,.25),rotation_regularization=.25,iterations=iterations,
  wall_seconds=wall_seconds,fixed_at_targets=True)
 rotations,frame_proof=transport_linked_profiles(centres,fit['centres'],fit['rotations'],links,protected=anchors)
 q=fit['centres'][raw_owner]+scales[raw_owner,None]*np.einsum('nij,nj->ni',rotations[raw_owner],p-centres[raw_owner])
 # Restore exact source and exact whole-similarity words at the declared ends.
 fixed_raw=np.flatnonzero(np.isin(raw_owner,list(fixed_nodes)));moving_raw=np.flatnonzero(np.isin(raw_owner,list(moving_nodes)))
 q[fixed_raw]=p[fixed_raw];q[moving_raw]=p[moving_raw]@m[:3,:3].T+m[:3,3]
 raw_rotations=rotations[raw_owner].copy();core_evidence=None
 if core_method=='arap':
  from .arap_core import arap_core
  _,representative=np.unique(p,axis=0,return_index=True)
  hard_raw=np.unique(np.r_[np.flatnonzero(raw_owner<len(groups)),fixed_raw,moving_raw]);hard=np.unique(inv[hard_raw])
  posed,core_evidence=arap_core(physical,edges,hard,q[representative[hard]],q[representative],maximum_steps=32,wall_seconds=wall_seconds)
  q=posed[inv];q[fixed_raw]=p[fixed_raw];q[moving_raw]=p[moving_raw]@m[:3,:3].T+m[:3,3]
  # Fit the actual remaining source stars for shading transport only. These
  # rotations never substitute for the independent discrete surface gate.
  i,j=edges.T;old=physical[i]-physical[j];new=posed[i]-posed[j];covariance=np.zeros((len(physical),3,3));outer=np.einsum('ni,nj->nij',new,old);np.add.at(covariance,i,outer);np.add.at(covariance,j,outer)
  u,_,vt=np.linalg.svd(covariance);D=np.broadcast_to(np.eye(3),covariance.shape).copy();D[:,-1,-1]=np.linalg.det(u@vt);local=u@D@vt;free_raw=np.flatnonzero(~np.isin(inv,hard));raw_rotations[free_raw]=local[inv[free_raw]]
 attrs={k:v.copy()for k,v in mesh.attributes.items()}
 if 'normals'in attrs:
  normal=attrs['normals'].astype(float)*attrs['normalSigns'][:,None];normal=np.einsum('nij,nj->ni',raw_rotations,normal);normal/=np.linalg.norm(normal,axis=1)[:,None]
  attrs['normals']=np.rint(normal*32767).astype('<i2');attrs['normalSigns']=np.ones(len(q),dtype='i1')
  for key in ['normals','normalSigns']:attrs[key][fixed_raw]=mesh.attributes[key][fixed_raw]
 profiles=[]
 for i,raw in enumerate(guide['profiles']+ports):
  k=raw_owner[raw[0]];old=p[raw];new=q[raw];wanted=np.linalg.norm(old[:,None]-old[None],axis=2)*scales[k];actual=np.linalg.norm(new[:,None]-new[None],axis=2)
  error=float(np.abs(wanted-actual).max());assert error<1e-10
  profiles.append(dict(index=i,rawIds=raw.tolist(),scale=float(scales[k]),maximumPairDistanceError_m=error))
 proof=dict(method='source-incidence complete rigid sections and branch cores; anchored proper frame transport',
  fixedPorts=sorted(fixed_ports),movingPorts=sorted(moving_ports),fixedRawIds=fixed_raw.tolist(),movingRawIds=moving_raw.tolist(),
  fixedSourceControlsLiteral=True,movingControlsEqualWholeSimilarity=True,sourceFacesAndAttributesRetained=True,
  sourceNormalFieldTransport='proper section rotation; fixed attachment normal words literal',
  minimumSectionScale=float(scales.min()),maximumSectionScale=float(scales.max()),profiles=profiles,
  frameTransport=frame_proof,unanchoredComponents=orphans,sectionNodes=n,
  sourceClinicalAdjacencyInferred=False,discreteAndNeighborAcceptancePending=True,anatomicalAcceptance=False,
  solverElapsedSeconds=fit['elapsed_seconds'])
 proof.update(coreMethod=core_method,coreEvidence=core_evidence)
 proof['additionalSourceNeckSections']=recovered
 return Mesh(q,f,mesh.frame,attrs,mesh.source),proof
