"""Bounded native-course aortic wall diagnostic; coronary and heart joins pending.

Preserve every existing outer control/corner word and native open end.
The 0.9 mm radial thickness is illustrative, not histological acceptance.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.tissue import layer_arrangement,normal_words,pack
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.algorithms.tube_continuity import boundary_loops

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-native-aortic-finite-wall-v1';D.mkdir(exist_ok=False)
cp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cp.read_text());source_sha=digest_file(cp);assert source_sha=='a0d7b46e033d145617a8937e1aff489db90b7a2808aa76ec90192f155ad95518'
reader=AtlasReader();rows={p['id']:p for p in a['parts']};frame='model:local-ta98:assembled-metres';pid='ZA:Ascending aorta';outer=reader.read_part(a,rows[pid],frame=frame);v,f=outer.vertices,outer.faces
source=recover_station_rings(v,f);assert len(source)==1;rings=source[0]['rings'];assert len(rings)==25 and len(rings[0])==12
centers=np.array([v[r].mean(0)for r in rings]);inner=v.copy();thickness=.0009
for ring,center in zip(rings,centers):
 radial=v[ring]-center;radius=np.linalg.norm(radial,axis=1);assert np.min(radius)>thickness*5
 inner[ring]=v[ring]-radial/radius[:,None]*thickness
fields=np.array([v,inner]);positions,surfaces,cells,construction=layer_arrangement(fields,f,['aortic-wall']);assert construction['finite_wall_cells_pass']
origins=[];review=[dict(mesh=str(D/'source-outer.npz'),name='Native aortic outer reference, literal',color=[.75,.35,.32,1])]
outer.save(D/'source-outer.npz')
for branch in ['ZA:Left coronary artery','ZA:Right coronary artery']:
 m=reader.read_part(a,rows[branch],frame=frame);vv,iv=np.unique(m.vertices,axis=0,return_inverse=True);ff=iv[m.faces];mouths=boundary_loops(ff);root=max(mouths,key=lambda r:vv[r,1].mean());center=vv[root].mean(0);distance,face,witness=igl.point_mesh_squared_distance(center[None],v,f)
 archive=D/(branch.split(':')[-1].replace(' ','-')+'.npz');m.save(archive)
 origins.append(dict(part=branch,source_geometry_sha256=m.geometry_sha,root_ring_controls=len(root),root_ring_m=vv[root].tolist(),root_center_m=center.tolist(),root_selection='Most superior native open end; anatomical origin course requires confirmation',distance_to_native_aortic_outer_m=float(np.sqrt(distance[0])),native_aortic_face=int(face[0]),closest_aortic_point_m=witness[0].tolist(),coronary_ostium_constructed=False,coronary_lumen_join_verified=False,blind_wall_installation_permitted=False))
 review.append(dict(mesh=str(archive),name=branch.replace('ZA:','')+' fixed native source',color=[.84,.08,.06,1]))
plan=dict(schema='male-native-aortic-wall-construction/v1',source_catalogue_sha256=source_sha,source_part=pid,source_bank_sha256=outer.source['bank_sha256'],native_source_vertices=len(v),native_source_faces=len(f),native_station_count=25,ring_controls=12,outer_surface_and_both_end_rings_literal=True,source_uvs_and_normal_words_literal=True,heart_surfaces_and_pose_literal=True,all_native_vessel_courses_fixed=True,wall_thickness_m=thickness,wall_thickness_status='Illustrative radial engineering field, not verified histological thickness',interfaces=[dict(name='ventricular-root',kind='unjoined',reason='No independently isolated proper paired host annulus established in current installed wall'),dict(name='arch-root',kind='unjoined',reason='Distal source arch has no paired inner wall established')],coronary_origin_interfaces=origins,anatomical_acceptance=False,installable=False,remaining_limits=['Coronary ostia must be constructed without occluding their entry before any installation','Heart-side and arch-side finite tissue interfaces remain unjoined','No artery histological coat separation claimed'])
atomic_json(D/'construction-plan.json',plan)
words=normal_words(positions.astype('<f4').astype(float),surfaces);words['interface:0']=outer.attributes['normals'][f]
uvs={key:outer.attributes['uvs'][f]if key.startswith('interface:')else np.zeros((len(ff),3,2),dtype='<f4')for key,ff in surfaces.items()}
result=pack(positions,surfaces,cells,{'aortic-wall':'TA98ARTERIAL:male:ascending-aorta:finite-wall'},frame,D/'packed',corner_words=words,corner_uvs=uvs,construction=construction,max_pairs=250000,max_tests=2000000)
assert result['passed'];candidate=Mesh.load(D/'packed/TA98ARTERIAL_male_ascending-aorta_finite-wall.npz')
assert np.array_equal(positions[surfaces['interface:0']],v[f]);assert np.array_equal(words['interface:0'],outer.attributes['normals'][f]);assert np.array_equal(uvs['interface:0'],outer.attributes['uvs'][f])
review[0]=dict(mesh=str(D/'packed/TA98ARTERIAL_male_ascending-aorta_finite-wall.npz'),name='Native-course finite aortic wall; coronary ostia still pending',color=[.75,.35,.32,1])
axis=centers[-1]-centers[0];axis/=np.linalg.norm(axis);mid=centers[len(centers)//2]
atomic_json(D/'review-spec.json',dict(frame=frame,meshes=review,resolution=[1100,900],time_limit_seconds=45,views=[dict(name='front',direction=[0,0,1]),dict(name='right',direction=[-1,0,0]),dict(name='aortic axial diagnostic',direction=axis.tolist(),section=dict(normal=axis.tolist(),point=mid.tolist(),clear_outer=True)),dict(name='aortic longitudinal diagnostic',direction=[0,0,1],section=dict(normal=[0,0,1],point=mid.tolist(),clear_outer=True))]))
report=dict(schema='male-native-aortic-finite-wall-stage/v1',source_catalogue_sha256=source_sha,bank_sha256=result['bank_sha256'],bank_bytes=result['total_bytes'],candidate_geometry_sha256=candidate.geometry_sha,closed_finite_wall_passed=result['passed'],finite_wall_cells_pass=construction['finite_wall_cells_pass'],native_outer_source_faces_preserved=len(f),native_uv_and_normal_corners_literal=True,annular_ends_only=True,no_lumen_disks=True,all_existing_geometry_unchanged=True,coronary_ostia_pending=2,root_interfaces_pending=2,installable=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T)
atomic_json(D/'stage.json',report);print(json.dumps(report))
