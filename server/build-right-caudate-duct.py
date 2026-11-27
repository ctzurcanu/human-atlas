"""A communicating finite right caudate branch, using the real right duct.

The outer and fluid unions are constructed together. No cylinder is left
inside another duct and no solid disk spans the communicating lumen.
"""
from pathlib import Path
import argparse,importlib.util,json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
from geometry_tools.audit import audit
from geometry_tools.open3d_backend import library
from geometry_tools.algorithms.tube_sections import centripetal_curve,rotation_minimizing_frames
from geometry_tools.algorithms.tube_continuity import loft_faces
from geometry_tools.tissue import normal_words
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['build','pack'],required=True);ap.add_argument('--corrected-cap-winding',action='store_true');ap.add_argument('--local-wall-union',action='store_true');args=ap.parse_args()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/('right-caudate-duct-v226'if args.local_wall_union else'right-caudate-duct-v224'if args.corrected_cap_winding else'right-caudate-duct-v222');out.mkdir(exist_ok=True);start=time.monotonic()
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cat.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';old=parts['TA98HEPATIC:male:biliary'];wall=reader.read_part(a,old,frame=frame,catalogue_path=cat)
def solid(m):
 v,inv=np.unique(m.vertices,axis=0,return_inverse=True);s=md.Manifold(md.Mesh64(np.ascontiguousarray(v),np.ascontiguousarray(inv[m.faces],dtype=np.uint64)));assert s.status()==md.Error.NoError;return s
def mesh(s):
 assert s.status()==md.Error.NoError
 q=s.to_mesh64();return Mesh(np.asarray(q.vert_properties)[:,:3],np.asarray(q.tri_verts,dtype=np.int64),frame)
def tube(path,radius):
 frames=rotation_minimizing_frames(path);angles=np.arange(20)*2*np.pi/20
 vertices=(path[:,None]+radius*(np.cos(angles)[None,:,None]*frames[:,None,:,0]+np.sin(angles)[None,:,None]*frames[:,None,:,1])).reshape(-1,3)
 faces=loft_faces(len(path),20).tolist();vertices=np.r_[vertices,path[[0,-1]]]
 for j in range(20):
  if args.corrected_cap_winding:faces.extend([[(j+1)%20,j,len(vertices)-2],[(len(path)-1)*20+j,(len(path)-1)*20+(j+1)%20,len(vertices)-1]])
  else:faces.extend([[j,(j+1)%20,len(vertices)-2],[(len(path)-1)*20+j,len(vertices)-1,(len(path)-1)*20+(j+1)%20]])
 m=Mesh(vertices,np.asarray(faces),frame)
 # Frame convention fixes the circle's orientation; use positive volume.
 p=m.triangles
 if np.einsum('ij,ij->i',p[:,0],np.cross(p[:,1],p[:,2])).sum()<0:m=Mesh(m.vertices,m.faces[:,[0,2,1]],frame)
 return solid(m)
if args.stage=='build':
 assert not(out/'wall.npz').exists()
 marker=a['liverSectionIntegration'];raw=reader.asset(marker['bankUrl']).read_bytes();d=next(x for x in marker['exclusions']if x['id']=='biliary');v=np.ndarray((d['vertexCount'],3),'<f4',buffer=raw,offset=d['positions']).astype(float)+d['origin'];f=np.ndarray((d['indexCount']//3,3),'<u4',buffer=raw,offset=d['indices']);outer=Mesh(v,f,frame);os,ws=solid(outer),solid(wall);fluid=os-ws
 assert fluid.status()==md.Error.NoError and fluid.volume()>0
 region=reader.read_part(a,parts['TA98HEPATIC:male:region:A05.8.01.064'],frame=frame,catalogue_path=cat)
 endpoint=np.array([-.022,1.263,.020]);q=BVH(region).closest_point(endpoint);point=np.asarray(q['point']);hit=BVH(outer).closest_point(point);surface=np.asarray(hit['point']);tri=outer.triangles[hit['face_id']];normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);normal/=np.linalg.norm(normal)
 if np.dot(normal,endpoint-surface)<0:raise ValueError('Right-duct outer witness faces away from the planned caudate territory')
 root=surface-normal*.0012;waypoints=np.array([root,surface+normal*.0015,(surface+endpoint)*.5,endpoint]);path,_=centripetal_curve(waypoints,chord_tolerance_m=.00001,max_step_m=.0005)
 oldfluid=mesh(fluid);o3d=library();scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor(oldfluid.vertices.astype('<f4')),o3d.core.Tensor(oldfluid.faces.astype('<u4')));assert scene.compute_occupancy(o3d.core.Tensor(root[None].astype('<f4')),nthreads=1).numpy()[0]>.5,'Branch inlet must start within the existing patent duct fluid'
 branch_outer=tube(path,.00065);inner_path=path.copy();inner_path[-1]+=(path[-1]-path[-2])/np.linalg.norm(path[-1]-path[-2])*.00015;branch_fluid=tube(inner_path,.00048)
 merged_outer=os+branch_outer;merged_fluid=fluid+branch_fluid
 # Preserve the old wall away from the actual new junction. Reconstructing
 # it from a separately rounded render classifier introduces false slivers.
 result=(ws-branch_fluid)+((branch_outer-fluid)-branch_fluid)if args.local_wall_union else merged_outer-merged_fluid
 assert len(merged_fluid.decompose())==len(fluid.decompose()),'Disconnected caudate lumen'
 for name,s in [('wall',result),('outer',merged_outer),('fluid',merged_fluid)]:mesh(s).save(out/(name+'.npz'))
 atomic_json(out/'construction.json',dict(source_wall_sha256=wall.geometry_sha,source_outer_sha256=outer.geometry_sha,root_m=root.tolist(),endpoint_m=endpoint.tolist(),waypoints=waypoints.tolist(),outer_radius_m=.00065,lumen_radius_m=.00048,lumen_components_before=len(fluid.decompose()),lumen_components_after=len(merged_fluid.decompose()),parent_term='A05.8.01.064',new_term='A05.8.01.068',lumen_unplugged=True,elapsed_seconds=time.monotonic()-start))
 print(json.dumps(dict(stage='build',wall_triangles=result.num_tri(),elapsed_seconds=time.monotonic()-start)),flush=True)
else:
 assert not(out/'network.json').exists();m=Mesh.load(out/'wall.npz');origin=np.array([0,1.25,0]);v=(m.vertices-origin).astype('<f4');ff,_,conform=conform_quantized_surfaces(m.vertices,{'wall':m.faces},origin);m=Mesh(v.astype(float)+origin,ff['wall'],frame);m.save(out/'emitted-wall.npz');atomic_json(out/'conform.json',conform)
 checks=audit(m,require_closed=True,max_pairs=1000000,max_tests=12000000);atomic_json(out/'wall-check.json',checks)
 if not checks['passed']:raise ValueError(('New right caudate duct failed',checks['failures']))
 atomic_json(out/'network.json',dict(passed=True,wall=str(out/'emitted-wall.npz'),outer=str(out/'outer.npz'),fluid=str(out/'fluid.npz'),construction=str(out/'construction.json'),audited_emitted_geometry_sha256=m.geometry_sha,elapsed_seconds=time.monotonic()-start))
 print(json.dumps(dict(stage='pack',passed=True,triangles=len(m.faces),elapsed_seconds=time.monotonic()-start)),flush=True)
