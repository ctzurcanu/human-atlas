"""Finite corresponding walls for an explicitly named native cardiac root."""
from pathlib import Path
from collections import deque
import argparse,json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,close_with_shared_disks,signed_volume
from geometry_tools.algorithms.tube_sections import path_tangents,arclength,sweep_sections,transport_spatial_sections
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--part-id',default='ZA:Superior vena cava',choices=['ZA:Superior vena cava','ZA:Ascending aorta','ZA:Pulmonary trunk']);ap.add_argument('--wall-mm',type=float,default=.6);ap.add_argument('--retain-native-end',action='store_true');ap.add_argument('--round-rings',action='store_true');ap.add_argument('--straighten',type=float,default=0.);ap.add_argument('--spatial-profiles',action='store_true');args=ap.parse_args();assert 0<=args.straighten<=1;start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());pid=args.part_id;term={'ZA:Superior vena cava':'A12.3.03.001','ZA:Ascending aorta':'A12.2.03.001','ZA:Pulmonary trunk':'A12.2.01.001'}[pid];record=next(p for p in a['parts']if p['id']==pid);reader=AtlasReader();source=reader.read_part(a,record,frame='model:local-ta98:assembled-metres',catalogue_path=path);source.save(out/'source.npz');p,inv=np.unique(source.vertices,axis=0,return_inverse=True);f=inv[source.faces];loops=boundary_loops(f);assert len(loops)==2 and len(loops[0])==len(loops[1]);adj=[set()for _ in p]
for t in f:
 for i,j in zip(t,np.roll(t,-1)):adj[i].add(int(j));adj[j].add(int(i))
d=np.full(len(p),-1);d[loops[0]]=0;queue=deque(map(int,loops[0]))
while queue:
 i=queue.popleft()
 for j in adj[i]:
  if d[j]<0:d[j]=d[i]+1;queue.append(j)
original_positions=p.copy();centers=np.array([p[d==k].mean(0)for k in range(d.max()+1)])
transport=None
if args.straighten:
 tubes=recover_station_rings(p,f);assert len(tubes)==1
 rings,correspondence=align_station_correspondence(p,f,tubes[0]['rings']);fields=p[np.array(rings)][None];old_centers=fields[0].mean(1);s=arclength(old_centers);line=old_centers[0]+(old_centers[-1]-old_centers[0])*s[:,None]/s[-1];target=old_centers+(line-old_centers)*args.straighten;target[[0,-1]]=old_centers[[0,-1]]
 swept,frames,transport=transport_spatial_sections(fields,target)if args.spatial_profiles else sweep_sections(fields,target,dtype=np.float64);atomic_json(out/'section-transport.json',dict(transport=transport,source_correspondence=correspondence,straightening_fraction=args.straighten));assert transport['sectionwise_isometry_pass']if args.spatial_profiles else transport['section_transport_pass'],transport
 p[np.array(rings).ravel()]=swept[0].reshape(-1,3)
elif args.round_rings:
 # Reorient each existing section by one proper rigid transform. Its centre,
 # radii and area stay literal; the folded source mouth no longer doubles back.
 tangents=path_tangents(centers)
 for k,tangent in enumerate(tangents):
  ids=np.flatnonzero(d==k);local=p[ids]-centers[k];_,_,vh=np.linalg.svd(local,full_matrices=False);normal=vh[-1]
  if normal@tangent<0:normal=-normal
  v=np.cross(normal,tangent);c=float(normal@tangent);skew=np.array([[0.,-v[2],v[1]],[v[2],0.,-v[0]],[-v[1],v[0],0.]])
  rotation=np.eye(3)+skew+skew@skew/(1+c);assert np.linalg.det(rotation)>.999999999
  p[ids]=local@rotation.T+centers[k]
inner=p.copy();thickness=args.wall_mm/1000;assert .0001<=thickness<=.001;centers=[]
for k in range(d.max()+1):
 ring=np.flatnonzero(d==k);assert len(ring)==len(loops[0]);c=p[ring].mean(0);radius=np.linalg.norm(p[ring]-c,axis=1);assert radius.min()>thickness*3;inner[ring]=c+(p[ring]-c)*(1-thickness/radius)[:,None];centers.append(c.tolist())
outer_v,outer_f,_=close_with_shared_disks(p,f);vv=list(inner);ff=list(f);caps=[]
for loop in loops:
 c=inner[loop].mean(0);_,_,vh=np.linalg.svd(p[loop]-c,full_matrices=False);direction=vh[-1];end=int(d[loop[0]]);neighbor=inner[d==(1 if end==0 else d.max()-1)].mean(0)
 if direction@(neighbor-c)>0:direction=-direction
 distance=max(0.,float(((p[loop]-c)@direction).max()-((inner[loop]-c)@direction).min()))+.0005;new=np.arange(len(vv),len(vv)+len(loop));vv.extend(inner[loop]+distance*direction);tip=len(vv);vv.append(c+distance*direction)
 for i,j,ii,jj in zip(loop,np.roll(loop,-1),new,np.roll(new,-1)):ff.extend([[j,i,ii],[j,ii,jj]]);caps.append([jj,ii,tip])
iv=np.array(vv);iff=np.r_[np.array(ff),np.array(caps)];origin=np.array([0.,1.375,.03125])
def solid(v,f):
 if signed_volume(v,f)<0:f=f[:,[0,2,1]]
 s=md.Manifold(md.Mesh64(v-origin,np.asarray(f,dtype='uint64'),tolerance=0.));assert str(s.status())=='Error.NoError';return s
fluid=Mesh.load(h/('blood-spaces-v130/fluid-boundary-4.npz'if pid=='ZA:Ascending aorta'else'blood-spaces-v130/fluid-boundary-3.npz'))
if args.retain_native_end:
 # Corresponding sheets already establish the finite tube. Direct annular
 # rims retain both nonplanar native mouths without a cap Boolean or slivers.
 n=len(p);rims=[]
 for loop in loops:
  for i,j in zip(loop,np.roll(loop,-1)):rims.extend([[j,i,i+n],[j,i+n,j+n]])
 mesh=Mesh(np.r_[p,inner],np.r_[f,f[:,[0,2,1]]+n,np.array(rims)],source.frame)
 body=solid(mesh.vertices,mesh.faces)
else:
 finite=solid(outer_v,outer_f)-solid(iv,iff)
 body=finite-solid(fluid.vertices,fluid.faces);z=body.to_mesh64()
 mesh=Mesh(z.vert_properties[:,:3]+origin,z.tri_verts,source.frame)
mesh.save(out/'finite-root.npz');checks=audit(mesh,require_closed=True,max_pairs=500000,max_tests=4000000);atomic_json(out/'audit.json',checks);atomic_npz(out/'source-walls.npz',outer_positions=p,inner_positions=inner,source_faces=f,station_indices=d,centers=np.array(centers));atomic_json(out/'construction.json',dict(passed=checks['passed'],source_part_id=pid,ta98_term=term,source_catalogue_sha256=digest_file(path),source_geometry_sha256=source.geometry_sha,wall_thickness_m=thickness,thickness_illustrative=True,source_outer_course_unchanged=not(args.round_rings or args.straighten),source_centers_retained=not bool(args.straighten),source_section_areas_retained=not args.spatial_profiles,source_spatial_profile_shapes_retained=True,source_end_centers_retained=True,straightening_fraction=args.straighten,section_transport=transport,maximum_source_section_rotation_displacement_m=float(np.linalg.norm(p-original_positions,axis=1).max()),native_axial_end_retained=args.retain_native_end,source_stations=len(centers),finite_volume_m3=float(body.volume()),triangles=len(mesh.faces),chamber_fluid_overlap_m3=float((body^solid(fluid.vertices,fluid.faces)).volume()),host_interface_pending=True,whole_lumen_checks_pending=True,elapsed_seconds=time.monotonic()-start,installed=False));print(json.dumps(dict(passed=checks['passed'],triangles=len(mesh.faces),failures=checks['failures'],elapsed_seconds=time.monotonic()-start)),flush=True)
raise SystemExit(0 if checks['passed']else 2)
