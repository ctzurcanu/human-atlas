"""Corresponding finite inner/outer wall charts for cardiac tissue peeling.

The 40 µm target is an illustrative construction allocation. Local facet
altitude and actual opposing-wall rays reduce it conservatively; this is not
a claim of measured endocardial/epicardial histological thickness.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.native_thyroid_capsule_field import incident_normal_field
from geometry_tools.algorithms.tube_continuity import finite_wall_cells
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';source=base/'whole-wall-v25';assert json.loads((source/'construction.json').read_text())['passed'];out=base/'wall-covering-v29';out.mkdir(exist_ok=False);m=Mesh.load(source/'whole-wall.npz');o=np.array([0.,1.375,.03125]);local=(m.vertices-o).astype('<f4');assert np.array_equal(local.astype(float)+o,m.vertices);delta,cones=incident_normal_field(local,m.faces,.00004,minimum_progress=.02,allow_feasible_lp_fallback=True);length=np.linalg.norm(delta,axis=1);direction=delta/length[:,None];scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor(local),o3d.core.Tensor(m.faces.astype('uint32')));epsilon=.0000003;rays=np.c_[local.astype(float)+direction*epsilon,direction].astype('float32');hit=scene.cast_rays(o3d.core.Tensor(rays),nthreads=1)['t_hit'].numpy();assert np.isfinite(hit).all()and np.all(hit>0),'No established opposing tissue wall at a carrier vertex'
q=m.triangles;heights=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)/np.linalg.norm(q-np.roll(q,1,axis=1),axis=2).max(1);local_height=np.full(len(local),np.inf)
for k in range(3):np.minimum.at(local_height,m.faces[:,k],heights)
fraction=np.minimum(1.,np.minimum(.1*(hit+epsilon),.1*local_height)/length);inner=(local.astype(float)+delta*fraction[:,None]).astype('<f4');fields=np.stack([local,inner]).astype(float)+o;construction=finite_wall_cells(fields,m.faces);atomic_npz(out/'fields.npz',fields=fields,faces=m.faces,frame=np.array(m.frame),geometry_origin=o);atomic_json(out/'cones.json',cones);atomic_json(out/'prisms.json',construction);print('Finite covering field saved',round(time.monotonic()-start,1),'s',flush=True)
core=Mesh(fields[1],m.faces,m.frame);shell=Mesh(fields.reshape(-1,3),np.r_[m.faces,m.faces[:,[0,2,1]]+len(local)],m.frame);core.save(out/'core.npz');shell.save(out/'covering.npz');checks={}
for name,mesh in [('core',core),('covering',shell)]:
 proof=audit(mesh,require_closed=True,max_pairs=2500000,max_tests=20000000);atomic_json(out/(name+'-audit.json'),proof);checks[name]=dict(passed=proof['passed'],failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs']);print(name,json.dumps(checks[name]),flush=True)
thickness=np.linalg.norm(fields[1]-fields[0],axis=1);report=dict(passed=all(x['passed']for x in checks.values())and construction['finite_wall_cells_pass'],source_sha256=digest_file(source/'whole-wall.npz'),source_geometry_sha256=m.geometry_sha,outer_controls_literal=True,nominal_target_m=.00004,thickness_quantiles_m=np.quantile(thickness,[0,.1,.5,.9,1]).tolist(),opposing_wall_ray_distance_quantiles_m=np.quantile(hit,[0,.1,.5,.9,1]).tolist(),checks=checks,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
