"""One bounded <=1 mm rigid-station alternative; no radius or topology change."""
from pathlib import Path
import sys,time,json
import numpy as np
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.tube_continuity import finite_wall_cells
from geometry_tools.contact import self_contact_summary

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-minimal-station-fit-v1';D.mkdir(exist_ok=False);source=Mesh.load(ROOT/'.local-models/geometry-workbench/heart/male-native-aortic-source-defects-v2/native-ascending-aorta.npz');v,f=source.vertices,source.faces;b=recover_station_rings(v,f)[0];rings,replay=align_station_correspondence(v,f,b['rings']);centers=np.array([v[r].mean(0)for r in rings]);radius=np.array([np.linalg.norm(v[r]-c,axis=1).max()for r,c in zip(rings,centers)]);base_inner=v.copy();
for r,c in zip(rings,centers):
 d=v[r]-c;base_inner[r]=v[r]-.0009*d/np.linalg.norm(d,axis=1)[:,None]
moving=list(range(1,14));nvar=len(moving)*6;area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)*.0009
def fields(x):
 outer=v.copy();inner=base_inner.copy()
 for j,i in enumerate(moving):
  q=x[j*6:j*6+6];r=rings[i];rotation=Rotation.from_rotvec(q[3:]*.001/radius[i]).as_matrix();translation=q[:3]*.001;outer[r]=(v[r]-centers[i])@rotation.T+centers[i]+translation;inner[r]=(base_inner[r]-centers[i])@rotation.T+centers[i]+translation
 return outer,inner
def constraints(x):
 outer,inner=fields(x);aa,bb=outer[f],inner[f];delta=bb-aa;u,w=aa[:,1]-aa[:,0],aa[:,2]-aa[:,0];du,dw=delta[:,1]-delta[:,0],delta[:,2]-delta[:,0];c0=np.einsum('fi,fvi->fv',np.cross(u,w),delta);c1=np.einsum('fi,fvi->fv',np.cross(du,w)+np.cross(u,dw),delta);c2=np.einsum('fi,fvi->fv',np.cross(du,dw),delta);return np.array([(c0+c1*t+c2*t*t)/area[:,None]for t in [0.,.25,.5,.75,1.]]).ravel()
x0=np.zeros(nvar);value=constraints(x0);jac=np.empty((len(value),nvar));step=1e-4
for j in range(nvar):
 x=x0.copy();x[j]=step;jac[:,j]=(constraints(x)-value)/step
# Six-component L1 norm <=1 at each station guarantees every actual point
# moves <=1 mm under the exact proper rotation and translation together.
A=np.c_[jac,np.zeros_like(jac)];rhs=-1e-4-value;absolute=[];ab=[]
for j in range(nvar):
 row=np.zeros(nvar*2);row[j]=1;row[nvar+j]=-1;absolute.extend([row.copy(),np.r_[-row[:nvar],row[nvar:]]]);ab.extend([0.,0.])
for j in range(len(moving)):
 row=np.zeros(nvar*2);row[nvar+j*6:nvar+j*6+6]=1;absolute.append(row);ab.append(1.)
atomic_json(D/'construction-plan.json',dict(schema='male-aortic-minimal-rigid-stations/v1',source_geometry_sha256=source.geometry_sha,source_faces_and_uvs_literal=True,source_station_count=25,exact_source_phase_and_diagonal_replay=True,replayed_source_intervals=len(replay),replayed_source_faces=sum(len(q['sourceFaceIds'])for q in replay),moving_stations=moving,fixed_stations=[0,*range(14,25)],maximum_point_move_m=.001,positive_proper_rotations_only=True,no_radius_reduction=True,all_heart_and_other_vessels_fixed=True,whole_context_acceptance_pending=True,anatomical_acceptance=False))
solve=linprog(np.r_[np.zeros(nvar),np.ones(nvar)],A_ub=np.r_[A,np.array(absolute)],b_ub=np.r_[rhs,ab],bounds=[(-1,1)]*nvar+[(0,1)]*nvar,method='highs',options={'time_limit':10.})
report=dict(schema='male-aortic-minimal-rigid-station-result/v1',linear_solver_success=bool(solve.success),linear_solver_status=int(solve.status),linear_solver_message=str(solve.message),station_budget='Each station: combined translation/rotation L1 <=1 mm equivalent; exact world displacement checked separately',source_phase_diagonal_replay_passed=True,elapsed_seconds=time.monotonic()-T,installed=False,anatomical_acceptance=False)
if solve.success:
 x=solve.x[:nvar];outer,inner=fields(x);displacements=np.linalg.norm(outer-v,axis=1);assert displacements.max()<=.001+1e-12;attrs={k:q.copy()for k,q in source.attributes.items()};norm=attrs['normals'].astype(float)
 for j,i in enumerate(moving):norm[rings[i]]=norm[rings[i]]@Rotation.from_rotvec(x[j*6+3:j*6+6]*.001/radius[i]).as_matrix().T
 attrs['normals']=np.rint(norm).astype('<i2');candidate=Mesh(outer,f,source.frame,attrs);candidate.save(D/'candidate-outer.npz');atomic_npz(D/'wall-fields.npz',fields=np.array([outer,inner]),faces=f,frame=np.asarray(source.frame),rigid_station_parameters=x,rings=np.array(rings));cross,_=self_contact_summary(candidate,example_limit=1000);prism=finite_wall_cells(np.array([outer,inner]),f);atomic_json(D/'self-contact.json',cross);atomic_json(D/'prism-check.json',prism);report.update(maximum_move_m=float(displacements.max()),source_self_crossings=25,candidate_self_crossings=cross['crossing_face_pairs'],invalid_finite_wall_faces=prism['invalid_counts_by_interval'],passed=bool(cross['defect_pairs']==0 and prism['finite_wall_cells_pass']),candidate_geometry_sha256=candidate.geometry_sha)
else:report.update(passed=False,source_outer_defects_retained=25,reason='No linearized positive-wall candidate satisfies the conservative rigid-station 1 mm bound; this is not a proof against all nonlinear 1 mm motions.')
atomic_json(D/'result.json',report);print(json.dumps(report))
