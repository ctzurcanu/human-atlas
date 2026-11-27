"""One relaxed linearized bound and literal port distances; never installs."""
from pathlib import Path
import sys,json,time,ast
import numpy as np
import igl
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.tube_continuity import boundary_loops,finite_wall_cells
from geometry_tools.contact import self_contact_summary

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-minimal-station-fit-v1/obstacle-diagnosis-v1';D.mkdir(exist_ok=False);source=Mesh.load(ROOT/'.local-models/geometry-workbench/heart/male-native-aortic-source-defects-v2/native-ascending-aorta.npz');v,f=source.vertices,source.faces;b=recover_station_rings(v,f)[0];rings,replay=align_station_correspondence(v,f,b['rings']);centers=np.array([v[r].mean(0)for r in rings]);radius=np.array([np.linalg.norm(v[r]-c,axis=1).max()for r,c in zip(rings,centers)]);base_inner=v.copy()
for r,c in zip(rings,centers):
 delta=v[r]-c;base_inner[r]=v[r]-.0009*delta/np.linalg.norm(delta,axis=1)[:,None]
moving=list(range(1,14));nvar=len(moving)*6;area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)*.0009;module=ast.parse((ROOT/'server/stage-male-aortic-minimal-station-fit.py').read_text());definitions=[node for node in module.body if isinstance(node,ast.FunctionDef)and node.name in ['fields','constraints']];exec(compile(ast.Module(body=definitions,type_ignores=[]),'reviewed-rigid-station-functions','exec'));x0=np.zeros(nvar);value=constraints(x0);jac=np.empty((len(value),nvar));step=1e-4
for j in range(nvar):
 x=x0.copy();x[j]=step;jac[:,j]=(constraints(x)-value)/step
rows=[]
for j,i in enumerate(moving):
 for delta in (v[rings[i]]-centers[i]):
  row=np.zeros((3,nvar));row[:,j*6:j*6+3]=np.eye(3)
  for axis in range(3):row[:,j*6+3+axis]=np.cross(np.eye(3)[axis],delta)/radius[i]
  rows.extend(row)
motion=np.array(rows);A=np.r_[np.c_[jac,np.zeros(len(jac))],np.c_[motion,-np.ones(len(motion))],np.c_[-motion,-np.ones(len(motion))]];rhs=np.r_[-1e-4-value,np.zeros(len(motion)*2)];solve=linprog(np.r_[np.zeros(nvar),1.],A_ub=A,b_ub=rhs,bounds=[(-10,10)]*nvar+[(0,10)],method='highs',options={'time_limit':10.});bound=dict(success=bool(solve.success),status=int(solve.status),message=str(solve.message),meaning='Minimum max per-control Cartesian component in the first-order positive-wall model; not a rigorous nonlinear lower bound')
if solve.success:
 x=solve.x[:-1];outer,inner=fields(x);maximum=float(np.linalg.norm(outer-v,axis=1).max());bound.update(linearized_component_budget_m=float(solve.x[-1]*.001),exact_proper_transform_maximum_move_m=maximum,within_user_one_mm_bound=bool(maximum<=.001));atomic_npz(D/'linearized-diagnostic-fields.npz',fields=np.array([outer,inner]),faces=f,parameters=x,frame=np.asarray(source.frame));test=finite_wall_cells(np.array([outer,inner]),f);cross,_=self_contact_summary(Mesh(outer,f,source.frame),example_limit=1000);bound.update(exact_source_self_crossings=cross['crossing_face_pairs'],exact_invalid_prisms=test['invalid_counts_by_interval']);
 if maximum<=.001 and cross['defect_pairs']==0 and test['finite_wall_cells_pass']:Mesh(outer,f,source.frame).save(D/'within-bound-candidate-outer.npz')
context_path=ROOT/'.local-models/geometry-workbench/jobs/b9c8595514cdc1b8671263ec425dbf7bbad8f059c484638ae4b24211c0b2d9dc/prepared-2/proposal.json';a=json.loads(context_path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();port=[];prox=v[rings[0]];far=v[rings[-1]]
for pid in ['TA98HEART:male:endocardium','TA98HEART:male:epicardium']:
 m=reader.read_part(a,parts[pid],frame=source.frame);distance,face,point=igl.point_mesh_squared_distance(prox,m.vertices,m.faces);mm=np.sqrt(distance)*1000;port.append(dict(part=pid,native_aortic_proximal_ring_control_distances_mm=mm.tolist(),minimum_mm=float(mm.min()),median_mm=float(np.median(mm)),maximum_mm=float(mm.max()),nearest_source_faces=face.tolist(),anatomical_pairing_inferred=False))
arch=reader.read_part(a,parts['ZA:Aortic arch'],frame=source.frame);vv,iv=np.unique(arch.vertices,axis=0,return_inverse=True);loops=boundary_loops(iv[arch.faces]);match=min(loops,key=lambda r:np.linalg.norm(vv[r].mean(0)-far.mean(0)));dist=np.linalg.norm(far[:,None]-vv[match][None],axis=2);port.append(dict(part='ZA:Aortic arch',matched_native_outer_control_count=len(match),native_aortic_distal_outer_ring_max_nearest_distance_mm=float(dist.min(1).max()*1000),native_aortic_distal_outer_ring_center_distance_mm=float(np.linalg.norm(far.mean(0)-vv[match].mean(0))*1000),paired_inner_wall_exists=False,anatomical_attachment_inferred=False));report=dict(schema='male-aortic-minimal-motion-obstacle/v1',source_phase_and_native_diagonals_exact=True,source_self_crossings=25,source_crossing_faces=26,folded_station_intervals=[0,1,2,8,9,10,11,12],bounded_attempt='../result.json',relaxed_linearized_bound=bound,literal_port_distances=port,heart_and_all_native_courses_unchanged=True,no_candidate_installed=True,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(D/'evidence.json',report);print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'],relaxed_bound=bound,port_summary=[{k:q[k]for k in q if k in ['part','minimum_mm','median_mm','maximum_mm','native_aortic_distal_outer_ring_max_nearest_distance_mm','native_aortic_distal_outer_ring_center_distance_mm']}for q in port])))
