"""One bounded aortic collar candidate: proper ring transport, fixed mouth/far field."""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.spatial.transform import Rotation,Slerp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.tube_continuity import finite_wall_cells

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-section-transport-v1';D.mkdir(exist_ok=False)
cp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cp.read_text());assert digest_file(cp)=='a0d7b46e033d145617a8937e1aff489db90b7a2808aa76ec90192f155ad95518';p=next(p for p in a['parts']if p['id']=='ZA:Ascending aorta');reader=AtlasReader();source=reader.read_part(a,p,frame='model:local-ta98:assembled-metres');v,f=source.vertices,source.faces;branch=recover_station_rings(v,f)[0];rings,correspondence=align_station_correspondence(v,f,branch['rings']);centers=np.array([v[r].mean(0)for r in rings]);frames=[]
for i,r in enumerate(rings):
 q=v[r]-centers[i];_,_,vt=np.linalg.svd(q,full_matrices=False);n=vt[-1];tangent=centers[min(i+1,len(centers)-1)]-centers[max(i-1,0)]
 if n@tangent<0:n=-n
 u=q[0]-n*(q[0]@n);u/=np.linalg.norm(u);w=np.cross(n,u);frames.append(np.c_[u,w,n])
frames=np.array(frames);end=14;arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(centers[:end+1],axis=0),axis=1))];parameters=arc/arc[-1];slerp=Slerp([0,1],Rotation.from_matrix(frames[[0,end]]));target_frames=slerp(parameters).as_matrix();length=np.linalg.norm(centers[end]-centers[0]);start_tangent=frames[0,:,2]*length;end_tangent=frames[end,:,2]*length;new=v.copy();rotations={}
for i in range(1,end):
 t=parameters[i];center=(2*t**3-3*t*t+1)*centers[0]+(t**3-2*t*t+t)*start_tangent+(-2*t**3+3*t*t)*centers[end]+(t**3-t*t)*end_tangent;rotation=target_frames[i]@frames[i].T;r=rings[i];new[r]=(v[r]-centers[i])@rotation.T+center;rotations[i]=rotation
attributes={k:q.copy()for k,q in source.attributes.items()};normals=attributes['normals'].astype(float)/32767
for i,rotation in rotations.items():normals[rings[i]]=normals[rings[i]]@rotation.T
attributes['normals']=np.rint(normals*32767).astype('<i2');candidate=Mesh(new,f,source.frame,attributes);candidate.save(D/'candidate-outer.npz');source.save(D/'source-outer.npz');crossing,_=self_contact_summary(candidate,example_limit=1000);inner=new.copy();areas=[];spans=[]
for r in rings:
 before=v[r]-v[r].mean(0);after=new[r]-new[r].mean(0);areas.append(float(np.linalg.norm(np.cross(after,np.roll(after,-1,axis=0)).sum(0))/np.linalg.norm(np.cross(before,np.roll(before,-1,axis=0)).sum(0))));spans.append(float(np.linalg.norm(after,axis=1).max()/np.linalg.norm(before,axis=1).max()));radius=np.linalg.norm(after,axis=1);inner[r]=new[r]-.0009*after/radius[:,None]
prisms=finite_wall_cells(np.array([new,inner]),f);fixed=np.concatenate([rings[0],*rings[end:]]);assert np.array_equal(new[fixed],v[fixed]);assert np.array_equal(attributes['normals'][fixed],source.attributes['normals'][fixed]);assert np.array_equal(attributes['uvs'],source.attributes['uvs']);np.savez_compressed(D/'wall-fields.npz',fields=np.array([new,inner]),faces=f,frame=np.asarray(source.frame),rings=np.array(rings));report=dict(schema='male-aortic-proper-section-transport/v1',source_catalogue_sha256=digest_file(cp),source_geometry_sha256=source.geometry_sha,candidate_geometry_sha256=candidate.geometry_sha,changed_stations=list(range(1,end)),fixed_station_profiles=[0,*range(end,len(rings))],native_faces_and_uvs_literal=True,normal_words_properly_rotated_only_in_support=True,all_heart_and_other_vessels_literal=True,maximum_move_m=float(np.linalg.norm(new-v,axis=1).max()),section_area_ratio_range=[min(areas),max(areas)],section_max_radius_ratio_range=[min(spans),max(spans)],outer_self_intersections=crossing,paired_radial_wall_cells=prisms,thickness_m=.0009,thickness_status='Illustrative',passed=bool(crossing['defect_pairs']==0 and prisms['finite_wall_cells_pass']),coronary_ostia_unjoined=True,heart_and_arch_interfaces_unjoined=True,installable=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(D/'stage.json',report)
print(json.dumps({k:v for k,v in report.items()if k not in ['outer_self_intersections','paired_radial_wall_cells']}));print('outer_defects',crossing['defect_pairs'],'prism_invalid',prisms['invalid_counts_by_interval'])
