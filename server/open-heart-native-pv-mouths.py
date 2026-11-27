"""Open the actual native pulmonary-vein lumens through atrial tissue only.

Native source paths/outer vessels stay fixed. Temporary cutter end disks are
computational only and must not contribute a face inside myocardial tissue.
"""
from pathlib import Path
import argparse,sys,time,json
from collections import deque
import numpy as np
import manifold3d as md
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.annular_pipe_envelopes import face_components,boundary_loops,close_with_shared_disks,signed_volume
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='native-pv-mouths-v60');ap.add_argument('--extend-ends',action='store_true');ap.add_argument('--atrial-fluid-end',action='store_true');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());m=Mesh.load(b/'whole-wall-v25/whole-wall.npz');o=np.array([0.,1.375,.03125]);solid=md.Manifold(md.Mesh64(m.vertices-o,m.faces.astype('uint64')));reader=AtlasReader();sources=[];cutters=[];disks=[];ring_records=[];cap_checks=[]
source_scene=o3d.t.geometry.RaycastingScene(nthreads=1);source_scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces.astype('uint32')));source_bvh=BVH(m);gates=json.loads((b/'blood-space-v12/compartments.json').read_text())['gates'];atrial_seed=np.array([.018,1.405,-.012]);seed_clearance=float(source_scene.compute_signed_distance(o3d.core.Tensor((atrial_seed-o).astype('float32')[None]),nthreads=1,nsamples=3).numpy()[0]);assert seed_clearance>.0003
for part in a['parts']:
 if part['id']not in {'ZA:Left inferior pulmonary vein','ZA:Left superior pulmonary vein','ZA:Right inferior pulmonary vein','ZA:Right superior pulmonary vein'}:continue
 v=reader.read_part(a,part,frame=m.frame,catalogue_path=path);points,inv=np.unique(v.vertices,axis=0,return_inverse=True);f=inv[v.faces];nc,label=face_components(f);sources.append(dict(part_id=part['id'],geometry_sha256=v.geometry_sha,components=nc))
 for j in range(nc):
  fs=f[label==j];used,ii=np.unique(fs,return_inverse=True);p=points[used];fs=ii.reshape(-1,3);rings=boundary_loops(fs);assert len(rings)==2;adj=[set()for _ in p]
  for face in fs:
   for k in range(3):adj[face[k]].add(int(face[(k+1)%3]));adj[face[(k+1)%3]].add(int(face[k]))
  distance=np.full(len(p),-1);distance[rings[0]]=0;queue=deque(map(int,rings[0]))
  while queue:
   k=queue.popleft()
   for q in adj[k]:
    if distance[q]<0:distance[q]=distance[k]+1;queue.append(q)
  groups=[np.flatnonzero(distance==d)for d in range(distance.max()+1)];assert all(len(g)==len(rings[0])for g in groups);inner=p.copy()
  for g in groups:
   center=p[g].mean(0);inner[g]=center+.8*(p[g]-center)
  if not args.extend_ends:vv,ff,dd=close_with_shared_disks(inner,fs,expected_mouths=2)
  else:
   verts=list(inner);new_faces=list(fs);disk_faces=[]
   for loop in boundary_loops(fs):
    end=int(distance[loop[0]]);center=inner[loop].mean(0);other=inner[groups[1 if end==0 else-2]].mean(0);direction=center-other;direction/=np.linalg.norm(direction);cap0=inner[loop];done=False
    gate=next(g for g in gates if g['part_id']==part['id']);is_atrial_end=np.linalg.norm(center-np.array(gate['center']))<.0001
    for extension in np.arange(.002,.010001,.001):
     if args.atrial_fluid_end and is_atrial_end:
      scale=min(.3,.25*seed_clearance/np.linalg.norm(cap0-center,axis=1).max());cap=atrial_seed+scale*(cap0-center)
     else:
      if args.atrial_fluid_end:direction=center-atrial_seed;direction/=np.linalg.norm(direction)
      cap=cap0+direction*extension
     cap_center=cap.mean(0);cap_vertices=np.r_[cap,cap_center[None]];cf=np.array([[k,(k+1)%len(cap),len(cap)]for k in range(len(cap))]);cap_samples=np.r_[cap,cap_center[None],cap[cf[:,:2]].mean(1),cap_vertices[cf].mean(1)];signed=source_scene.compute_signed_distance(o3d.core.Tensor((cap_samples-o).astype('float32')),nthreads=1,nsamples=3).numpy()
     if signed.min()<.00002:continue
     _,touch=pair_contacts(Mesh(cap_vertices,cf,m.frame),m,0.,second_bvh=source_bvh,max_pairs=1000,max_tests=100000)
     if touch['returned_triangle_pairs']or touch['truncated']:continue
     done=True;break
    assert done,('Cutter cap remains in tissue',part['id'],j,end)
    new=np.arange(len(verts),len(verts)+len(loop));verts.extend(cap);tip=len(verts);verts.append(cap_center)
    for k in range(len(loop)):
     a0,b0=int(loop[k]),int(loop[(k+1)%len(loop)]);aa,bb=int(new[k]),int(new[(k+1)%len(loop)]);new_faces.extend([[b0,a0,aa],[b0,aa,bb]]);disk_faces.append([bb,aa,tip])
    cap_checks.append(dict(part_id=part['id'],component=j,end=end,extension_m=float(np.linalg.norm(cap_center-center)),computation_only_atrial_cap=args.atrial_fluid_end and is_atrial_end,all_samples_outside_myocardium=True,minimum_sample_clearance_m=float(signed.min()),exact_cap_tissue_contacts=touch['returned_triangle_pairs']))
   vv=np.array(verts);dd=np.array(disk_faces);ff=np.r_[np.array(new_faces),dd]
  if signed_volume(vv,ff)<0:ff=ff[:,[0,2,1]]
  cutter=md.Manifold(md.Mesh64(vv-o,ff.astype('uint64')));assert str(cutter.status())=='Error.NoError';cutters.append(cutter);disks.append(vv[dd]);ring_records.append(dict(part_id=part['id'],component=j,points_m=inner.tolist(),faces=fs.tolist(),inset_fraction=.8))
combined=md.Manifold.batch_boolean(cutters,md.OpType.Add);result=solid-combined;assert str(result.status())=='Error.NoError';removed=float((solid.volume()-result.volume())*1e6);assert 0<removed<8.,removed;z=result.to_mesh64();r=normalize_shared_boundary(z.vert_properties[:,:3]+o,{'wall':z.tri_verts},o);wall=Mesh(r['positions'],r['cells']['wall'],m.frame);wall.save(out/'whole-wall.npz');atomic_json(out/'encoding.json',r['evidence']);atomic_json(out/'source-courses.json',dict(sources=sources,rings=ring_records,cap_checks=cap_checks));print('Opened native PV lumen courses; removed tissue mL',removed,'faces',len(wall.faces),flush=True);proof=audit(wall,require_closed=True,max_pairs=3000000,max_tests=22000000);atomic_json(out/'emitted-audit.json',proof)
scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((wall.vertices-o).astype('float32')),o3d.core.Tensor(wall.faces.astype('uint32')));disc_centres=np.concatenate(disks).mean(1);dist=scene.compute_distance(o3d.core.Tensor((disc_centres-o).astype('float32')),nthreads=1).numpy();report=dict(passed=proof['passed'],source_catalogue_sha256=digest_file(path),source_geometry_sha256=m.geometry_sha,geometry_sha256=wall.geometry_sha,removed_myocardial_tissue_mL=removed,inset_fraction=.8,native_outer_vessels_moved=False,temporary_cap_centroid_minimum_distance_to_tissue_m=float(dist.min()),triangles=len(wall.faces),failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
