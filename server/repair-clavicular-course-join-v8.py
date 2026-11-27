"""One C1 join repair on the saved, obstacle-cleared route; no route search."""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames

start=time.monotonic();work=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracoacromial-clavicular-branch-v1'
source=work/'left-native-course-router-v7';out=work/'left-native-course-joined-v8';out.mkdir(exist_ok=False)
archive=source/'candidate-courses-and-source-witnesses.npz';z=dict(np.load(archive))
a,b=195,214;c=z['main_centres'];fields=z['main_fields'];original_fields=fields.copy();original_centres=c.copy()
plan=dict(schema='clavicular-c1-source-bound-join-repair/v8',source_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),changed_stations=[a+1,b-1],fixed_endpoint_stations=[a,b],all_other_centres_and_profiles_literal=True,fixed_radii=True,method='Cubic Hermite endpoint positions and actual source section normals; rotation-minimizing transport with continuous angular correction to the literal target profile.',reason='61 self-crossing pairs localize to the old spline reversal immediately before the preserved distal join. No foreign-context crossing was detected in the saved route.',no_parameter_search=True,installed=False)
(out/'plan-before-join-repair.json').write_text(json.dumps(plan,indent=2)+'\n')
def normal(k):
 r=fields[k]-c[k];n=np.cross(r[0],r[1]);n/=np.linalg.norm(n)
 if n@(c[min(k+1,len(c)-1)]-c[max(k-1,0)])<0:n=-n
 return n
length=np.linalg.norm(c[b]-c[a]);spline=CubicHermiteSpline([0,length],[c[a],c[b]],[normal(a),normal(b)])
t=np.linspace(0,length,b-a+1);curve=spline(t);tangent=spline(t,1);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
frames=rotation_minimizing_frames(curve,tangent,initial_reference=fields[a,0]-c[a]);target=fields[b,0]-c[b];target/=np.linalg.norm(target)
phase=np.arctan2(target@frames[-1,:,1],target@frames[-1,:,0]);angle=np.arange(16)[None]*np.pi/8+np.linspace(0,phase,len(curve))[:,None]
replacement=curve[:,None]+z['main_radii'][a:b+1,None,None]*(np.cos(angle)[:,:,None]*frames[:,None,:,0]+np.sin(angle)[:,:,None]*frames[:,None,:,1])
replacement[0]=fields[a];replacement[-1]=fields[b];curve[0]=c[a];curve[-1]=c[b]
c=c.copy();fields=fields.copy();c[a:b+1]=curve;fields[a:b+1]=replacement
assert np.array_equal(c[:a+1],original_centres[:a+1])and np.array_equal(c[b:],original_centres[b:])
assert np.array_equal(fields[:a+1],original_fields[:a+1])and np.array_equal(fields[b:],original_fields[b:])
z.update(main_centres=c,main_fields=fields,main_vertices=fields.reshape(-1,3),main_initial_centres=c.copy())
np.savez_compressed(out/'candidate-courses-and-source-witnesses.npz',**z)
for prefix in ['main','muscular']:Mesh(z[prefix+'_vertices'],z[prefix+'_faces'],'model:local-ta98:assembled-metres').save(out/(prefix+'-outer-course.npz'))
spec=json.loads((source/'course-spec.json').read_text());spec['course_archive']=str(out/'candidate-courses-and-source-witnesses.npz');spec['course_archive_sha256']=hashlib.sha256(Path(spec['course_archive']).read_bytes()).hexdigest();spec['join_repair_plan']=str(out/'plan-before-join-repair.json');spec['actual_render']=str(out/'actual-course-context-and-sections.png');spec['construction_authorized']=False
(out/'course-spec.json').write_text(json.dumps(spec,indent=2)+'\n')
check=audit(Mesh(z['main_vertices'],z['main_faces'],'model:local-ta98:assembled-metres'),require_closed=False,expected_components=1,max_pairs=600000,max_tests=2000000)
(out/'main-self-audit.json').write_text(json.dumps(check,indent=2)+'\n');assert check['passed'],check
proof=dict(passed=True,complete=True,only_prefix_join_changed=True,literal_v5_distal_and_twig=True,maximum_profile_motion_m=float(np.linalg.norm(fields-original_fields,axis=2).max()),elapsed_seconds=time.monotonic()-start)
(out/'join-repair-proof.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
