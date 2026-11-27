"""One bounded source-profile correction against the complete fixed femur."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,winding_contains
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.section_course_clearance import smooth_course_metric
from geometry_tools.algorithms.section_course_socp import solve_weighted_section_socp
from geometry_tools.algorithms.tube_collar_rounding import _rotation

BASE=(ROOT/sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1/course-correction-v1'

def tangent(points):
    t=np.gradient(points,axis=0);return t/np.linalg.norm(t,axis=1)[:,None]

def main():
    began=time.monotonic();plan_path=BASE/'plan.json';plan=json.loads(plan_path.read_bytes())
    for key in ['course','femur']:assert digest_file(ROOT/plan[key])==plan[key+'_sha256']
    data=np.load(ROOT/plan['course']);raw=np.load(ROOT/plan['femur']);frame='model:local-ta98:assembled-metres'
    bone=Mesh(raw['vertices'],raw['faces'],frame);assert audit(bone,require_closed=True,expected_components=1,expected_genus=[0])['passed']
    v=data['vertices'];f=data['faces'];rings=data['phase_aligned_ring_vertex_ids'];centers=data['centres'];old=data['profiles'];offset=old-centers[:,None]
    fixed=np.array(plan['fixed_source_root_and_ascending_stations']);station=data['station_of_vertex'];K=smooth_course_metric(len(rings),smoothing=plan['smoothing'],fixed=fixed)
    original_tangent=tangent(centers);index=SurfaceIndex(bone);d=np.zeros_like(centers);records=[];rows=[];face_rows=[]
    domain=np.nextafter(np.asarray(abs(v).max(0)+plan['maximum_motion_metres'],dtype='<f4'),np.float32(np.inf))
    storage_padding=float(np.linalg.norm(np.spacing(domain).astype(float)/2))+1e-10 if plan.get('account_for_float32_storage') else 0.
    report=dict(schema='native-lcfa-bone-clearance-course-fit/v1',plan_sha256=digest_file(plan_path),script_sha256=digest_file(Path(__file__)),iterations=records,solver_storage_padding_metres=storage_padding,passed=False,installed=False,anatomical_acceptance=False)
    def construct():
        nt=tangent(centers+d);rot=np.array([np.eye(3)if j in fixed else _rotation(a,b)for j,(a,b)in enumerate(zip(original_tangent,nt))]);q=np.einsum('sij,svj->svi',rot,offset)+centers[:,None]+d[:,None]
        q[fixed]=old[fixed];candidate=v.copy();candidate[rings]=q.astype('<f4').astype(float);return candidate,rot,nt
    try:
        for iteration in range(plan['maximum_linearizations']+1):
            candidate,rot,nt=construct();flat=candidate[rings].reshape(-1,3);membership=winding_contains(bone,flat,distance_backend='open3d');inside=np.array(membership['inside']);ambiguous=np.array(membership['boundary_ambiguous']);assert not ambiguous.any()
            nearest,_=index.closest_points(flat,certify=True);gap=nearest['distance_m'];normal=flat-nearest['points'];assert np.all(gap>1e-12);normal/=gap[:,None];normal[inside]*=-1;signed=np.where(inside,-gap,gap)
            mesh=Mesh(candidate,f,frame);hits,summary=pair_contacts(bone,mesh,0.);assert summary['complete_within_pair_budget'];bad=np.isin(hits['relation'],['crossing','coplanar-overlap','touching']);active=np.flatnonzero(signed<plan['clearance_metres'])
            state=dict(iteration=iteration,inside_controls=int(inside.sum()),minimum_signed_vertex_distance_metres=float(signed.min()),bone_contacts=summary,active_controls=len(active));records.append(state);atomic_json(BASE/'checkpoint.json',report,overwrite=(BASE/'checkpoint.json').exists())
            atomic_npz(BASE/f'candidate-checkpoint-{iteration:02}.npz',vertices=candidate,faces=f,frame=np.array(frame),source_vertices=v,rings=rings,source_profiles=old,source_centres=centers,station_displacements=d,station_rotations=rot,fixed_station_indices=fixed,source_to_exact_welded=data['source_to_exact_welded'])
            if not len(active) and not bad.any():report['passed']=True;break
            if iteration==plan['maximum_linearizations']:break
            for k in active:
                j,c=divmod(int(k),rings.shape[1]);rows.append((j,c,normal[k].copy(),nearest['points'][k].copy()))
            for bf,af in hits['pairs'][bad]:
                tri=bone.triangles[bf];n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n)
                for vertex in f[af]:
                    j=int(station[vertex]);c=int(np.flatnonzero(rings[j]==vertex)[0]);face_rows.append((j,c,n.copy(),tri[0].copy()))
            all_rows=rows+face_rows;W=np.zeros((len(all_rows),len(rings)));normals=[];bounds=[]
            for i,(j,c,n,q)in enumerate(all_rows):
                W[i,j]=1;normals.append(n);bounds.append(plan['clearance_metres']+storage_padding-n@(centers[j]+rot[j]@offset[j,c]-q))
            displacement,solver=solve_weighted_section_socp(K,W,normals,bounds,maximum_motion_m=plan['maximum_motion_metres'],time_limit_seconds=plan['solver_seconds_per_round'],maximum_iterations=100)
            state['solver']=solver
            if displacement is None:report['failure']='Bounded course halfspaces are infeasible or the solver did not complete';break
            d=displacement
        if report['passed']:
            quality=audit(mesh,require_closed=False,expected_components=1);report['whole_source_quality']=quality;report['passed'] &= quality['passed']
            protected=rings[fixed].ravel();assert np.array_equal(candidate[protected],v[protected]);before=old-old.mean(1)[:,None];after=candidate[rings]-candidate[rings].mean(1)[:,None]
            replay=np.einsum('sji,svj->svi',rot,after);error=np.linalg.norm(replay-before,axis=2);area=lambda p:np.cross(p,np.roll(p,-1,axis=1)).sum(1)/2
            old_area=abs(np.einsum('si,si->s',area(before),original_tangent));new_area=abs(np.einsum('si,si->s',area(after),nt));relative=abs(new_area-old_area)/old_area
            maximum_motion=float(np.linalg.norm(candidate-v,axis=1).max());report.update(protected_controls_literal=True,maximum_motion_metres=maximum_motion,maximum_rigid_profile_replay_error_metres=float(error.max()),maximum_projected_axial_area_change_relative=float(relative.max()),source_indices_literal=True,source_uvs_literal=True)
            report['passed'] &= maximum_motion<=plan['maximum_motion_metres'] and relative.max()<.001
            atomic_npz(BASE/'candidate.npz',vertices=candidate,faces=f,frame=np.array(frame),source_vertices=v,rings=rings,source_profiles=old,source_centres=centers,station_displacements=d,station_rotations=rot,fixed_station_indices=fixed,source_to_exact_welded=data['source_to_exact_welded'])
    except Exception as error:report['passed']=False;report['failure']=dict(type=type(error).__name__,message=str(error))
    report['elapsed_seconds']=time.monotonic()-began;atomic_json(BASE/'proof.json',report);print(json.dumps({k:report[k]for k in ['passed','elapsed_seconds','solver_storage_padding_metres']}|{'iterations':[dict(iteration=x['iteration'],inside=x['inside_controls'],crossings=x['bone_contacts']['crossing_face_pairs'],active=x['active_controls'],minimum_gap_m=x['minimum_signed_vertex_distance_metres'])for x in records],'maximum_motion_metres':report.get('maximum_motion_metres'),'maximum_area_change_relative':report.get('maximum_projected_axial_area_change_relative'),'failure':report.get('failure')}))

if __name__=='__main__':main()
