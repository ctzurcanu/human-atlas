"""Only the explicitly permitted five-course blood-union diagnostic.

No wall, fit, course displacement, catalogue or installation is produced.
"""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,topology,winding_contains
from geometry_tools.algorithms.branch_junction_envelopes import close_branch
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,signed_volume
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
import manifold3d as md
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';C=P/'three-course-finite-educational-v1';I=P/'temporal-source-incidence-v2';OUT=I/'five-course-fluid-diagnostic-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);g=json.loads((C/'generation-complete.json').read_text());origin=np.array(g['origin_m']);existing={};inputs=[];ids={};solids=[];op=0
 def checkpoint(name,arrays=None):
  nonlocal op
  op+=1;stem=f'{op:02d}-{name}'
  if arrays is not None:atomic_npz(OUT/(stem+'-input.npz'),**arrays)
  atomic_json(OUT/(stem+'-pending.json'),dict(operation=name,elapsed_seconds=time.monotonic()-start,input_archive=stem+'-input.npz'if arrays is not None else None));return stem
 for row in g['classifiers']:
  if not row['role'].endswith('fluid'):continue
  inputrow=next(r for r in g['native_inputs']if r['role']==row['role']);q=np.load(C/inputrow['checkpoint_input']);world=q['positions_mm']/1000.+origin;existing[row['course']]=Mesh(world,q['faces'],FRAME);stem=checkpoint('make-'+row['role'],dict(positions_mm=q['positions_mm'],faces=q['faces']));solid,oid=make_solid(q['positions_mm'],q['faces']);solids.append(solid);ids[row['role']]=oid;inputs.append(dict(role=row['role'],original_id=oid,disk_faces=row['disk_faces'],junction_faces=row['junction_faces'],mouths=row['mouths']));atomic_json(OUT/(stem+'-done.json'),dict(status=str(solid.status()),original_id=oid))
 tz=np.load(I/'complete-source-charts.npz');v=tz['temporal_vertices'];f=tz['temporal_faces'];containment=[]
 for index,rootstation,parent in [(0,0,'M3'),(1,48,'temporal0')]:
  role='temporal'+str(index);rings=tz['component'+str(index)+'_rings'];faceids=tz['component'+str(index)+'_source_face_ids'];used=np.unique(f[faceids]);lookup=np.full(len(v),-1,int);lookup[used]=np.arange(len(used));outer=v[used];inner=outer.copy();faces=lookup[f[faceids]]
  for ring in rings:
   q=v[ring];center=q.mean(0);R=float(np.linalg.norm(q-center,axis=1).mean());thickness=min(.00012*np.sqrt(R/.001),.35*R);inner[lookup[ring]]=center+(q-center)*((R-thickness)/R)
  loops=boundary_loops(faces);target=set(lookup[rings[rootstation]].tolist());which=next(i for i,loop in enumerate(loops)if set(loop.tolist())==target);_,fluid=close_branch(outer,inner,faces,{which:v[rings[rootstation]].mean(0)},origin=origin)
  if signed_volume(fluid['positions'],fluid['faces'])<0:fluid['faces']=fluid['faces'][:,[0,2,1]]
  mesh=Mesh(fluid['positions'],fluid['faces'],FRAME);mesh.save(OUT/(role+'-diagnostic-fluid.npz'));quality=audit(mesh,require_closed=True,expected_components=1,expected_genus=[0]);atomic_json(OUT/(role+'-input-audit.json'),quality);assert quality['passed'];fan=fluid['positions'][np.unique(fluid['faces'][fluid['junction_faces']])];proof=winding_contains(existing[parent],fan,boundary_tolerance_m=1e-10);assert all(proof['inside'])and not any(proof['boundary_ambiguous']);containment.append(dict(parent=parent,child=role,proof=proof));existing[role]=mesh;vp=(fluid['positions']-origin)*1000.;stem=checkpoint('make-'+role+'-fluid',dict(positions_mm=vp,faces=fluid['faces']));solid,oid=make_solid(vp,fluid['faces']);solids.append(solid);ids[role+'-fluid']=oid;inputs.append(dict(role=role+'-fluid',original_id=oid,disk_faces=fluid['disk_faces'].tolist(),junction_faces=fluid['junction_faces'].tolist(),mouths=fluid['mouths']));atomic_json(OUT/(stem+'-done.json'),dict(status=str(solid.status()),original_id=oid))
 atomic_json(OUT/'strict-internal-fan-containment.json',containment);stem=checkpoint('union-five-fluid-courses');joined=md.Manifold.batch_boolean(solids,md.OpType.Add);assert joined.status()==md.Error.NoError;checkpoint('extract-diagnostic-fluid-union');mesh=joined.to_mesh64();origins,faceids=provenance(mesh);vertices=np.asarray(mesh.vert_properties)[:,:3]/1000.+origin;faces=np.asarray(mesh.tri_verts,dtype=np.int64);diagnostic=Mesh(vertices,faces,FRAME);diagnostic.save(OUT/'five-course-fluid-union.npz');atomic_npz(OUT/'diagnostic-provenance.npz',original_ids=origins,source_face_ids=faceids,frame=np.array(FRAME));t=topology(diagnostic);ports=[];internal=[]
 for row in inputs:
  own=origins==row['original_id'];internal.extend(np.flatnonzero(own&np.isin(faceids,row['junction_faces'])).tolist())
  for mouth in row['mouths']:
   if mouth['kind']!='external-fluid-classifier':continue
   count=int(np.sum(own&np.isin(faceids,mouth['newFaceIds'])));ports.append(dict(role=row['role'],mouth=mouth['mouth'],surviving_classifier_disk_faces=count))
 result=dict(schema='five-course-blood-union-topology-diagnostic/v1',installed=False,no_wall_constructed=True,no_fit_or_source_geometry_modified=True,scope='Diagnostic occupied-fluid classifier only, not an anatomical candidate or authorization to merge outside-root contacts.',topology=t,expected_tree_genus=0,expected_terminal_classifier_mouths=6,terminal_classifier_mouths=ports,internal_root_fan_faces_surviving=len(internal),tree_topology_test_passed=t['closed_oriented_manifold']and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==0 and len(ports)==6 and all(r['surviving_classifier_disk_faces']>0 for r in ports),outside_parent_contacts_remain_disqualifying='Read common-parent-root-zone-v1/containment-proof.json; this diagnostic cannot authorize arbitrary outside-parent arterial communication.',elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'topology-proof.json',result);print(json.dumps(result))

if __name__=='__main__':main()
