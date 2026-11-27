"""One bounded inferior-bursa LCFA course; profiles and original roots stay fixed."""
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
from geometry_tools.solver_diagnostics import json_safe_solver_diagnostics

BASE=(ROOT/sys.argv[1]).resolve()
def tangent(p):
    a=np.gradient(p,axis=0);return a/np.linalg.norm(a,axis=1)[:,None]

def main():
    began=time.monotonic();plan=json.loads((BASE/'plan.json').read_bytes());records=[]
    report=dict(schema='native-lcfa-coupled-course-fit/v1',plan_sha256=digest_file(BASE/'plan.json'),script_sha256=digest_file(Path(__file__)),source_catalogue_sha256=plan['source_catalogue_sha256'],input_sha256={key:plan[key+'_sha256'] for key in ['course','initial_candidate','femur','bursa','tendon']},passed=False,installed=False,anatomical_acceptance=False,iterations=records)
    try:
        assert digest_file(ROOT/plan['source_catalogue'])==plan['source_catalogue_sha256']
        for k in ['course','initial_candidate','femur','bursa','tendon']:assert digest_file(ROOT/plan[k])==plan[k+'_sha256']
        z=np.load(ROOT/plan['course']);initial=np.load(ROOT/plan['initial_candidate']);frame='model:local-ta98:assembled-metres'
        raw=np.load(ROOT/plan['femur']);bone=Mesh(raw['vertices'],raw['faces'],frame);bursa=Mesh.load(ROOT/plan['bursa']);tendon=Mesh.load(ROOT/plan['tendon']);obstacles=[('femur',bone,True),('bursa',bursa,True),('tendon',tendon,False)]
        for name,mesh,closed in obstacles:
            if closed:assert audit(mesh,require_closed=True,expected_components=1,expected_genus=[0])['passed']
        indices={name:SurfaceIndex(mesh)for name,mesh,_ in obstacles}
        v=z['vertices'];f=z['faces'];rings=z['phase_aligned_ring_vertex_ids'];centers=z['centres'];profiles=z['profiles'];offset=profiles-centers[:,None];fixed=np.array(plan['fixed_source_root_and_ascending_stations']);station=z['station_of_vertex'];old_tangent=tangent(centers);d=initial['station_displacements'].copy();K=smooth_course_metric(len(rings),smoothing=plan['smoothing'],fixed=fixed)
        normals={name:np.cross(mesh.triangles[:,1]-mesh.triangles[:,0],mesh.triangles[:,2]-mesh.triangles[:,0])for name,mesh,_ in obstacles}
        for name in normals:normals[name]/=np.linalg.norm(normals[name],axis=1)[:,None]
        # Positive geometric side is fixed by source course samples adjacent to
        # both ends of the actual tendon-contact interval, not per-point flips.
        anchors=np.array(plan['tendon_side_source_anchor_stations']);near,_=indices['tendon'].closest_points(centers[anchors],certify=True);signed=np.einsum('ij,ij->i',centers[anchors]-near['points'],normals['tendon'][near['face_ids']]);assert np.all(signed>plan['clearance_metres']);report['tendon_side_anchors']=[dict(station=int(j),signed_geometric_side_metres=float(s),distance_metres=float(a))for j,s,a in zip(anchors,signed,near['distance_m'])]
        padding=float(np.linalg.norm(np.spacing(np.nextafter(np.asarray(abs(v).max(0)+plan['maximum_motion_metres'],dtype='<f4'),np.float32(np.inf))).astype(float)/2))+1e-10
        report['solver_storage_padding_metres']=padding
        for iteration in range(plan['maximum_linearizations']+1):
            nt=tangent(centers+d);rot=np.array([np.eye(3)if j in fixed else _rotation(a,b)for j,(a,b)in enumerate(zip(old_tangent,nt))]);q=np.einsum('sij,svj->svi',rot,offset)+centers[:,None]+d[:,None];q[fixed]=profiles[fixed];candidate=v.copy();candidate[rings]=q.astype('<f4').astype(float);flat=candidate[rings].reshape(-1,3);mesh=Mesh(candidate,f,frame);rows=[];state=dict(iteration=iteration,obstacles=[]);records.append(state);all_clear=True
            for name,obstacle,closed in obstacles:
                near,_=indices[name].closest_points(flat,certify=True);gap=near['distance_m'];assert np.all(gap>1e-12);normal=flat-near['points'];normal/=gap[:,None]
                if closed:
                    membership=winding_contains(obstacle,flat,distance_backend='open3d');inside=np.array(membership['inside']);assert not np.array(membership['boundary_ambiguous']).any();normal[inside]*=-1;signed=np.where(inside,-gap,gap);eligible=np.ones(len(flat),bool)
                else:
                    radial=normal.copy();normal=normals[name][near['face_ids']].copy();interior=np.min(near['barycentric'],axis=1)>1e-9;signed=np.einsum('ij,ij->i',flat-near['points'],normal);normal[~interior]=radial[~interior];signed[~interior]=gap[~interior];inside=np.zeros(len(flat),bool);eligible=gap<=plan['tendon_local_patch_distance_metres']
                hits,contact=pair_contacts(obstacle,mesh,0.);assert contact['complete_within_pair_budget'];bad=np.isin(hits['relation'],['crossing','coplanar-overlap','touching']);active=np.flatnonzero(eligible&(signed<plan['clearance_metres']));all_clear=all_clear and not len(active)and not bad.any()
                state['obstacles'].append(dict(name=name,inside_controls=int(inside.sum()),active_controls=len(active),minimum_local_signed_distance_metres=float(signed[eligible].min())if eligible.any()else None,contacts=contact))
                for k in active:
                    j,c=divmod(int(k),rings.shape[1]);rows.append((name,j,c,normal[k].copy(),near['points'][k].copy(),int(near['face_ids'][k]),bool(closed or interior[k])))
                for bf,af in hits['pairs'][bad]:
                    for vertex in f[af]:
                        j=int(station[vertex]);c=int(np.flatnonzero(rings[j]==vertex)[0]);rows.append((name,j,c,normals[name][bf].copy(),obstacle.triangles[bf,0].copy(),int(bf),True))
            # Shading normal inconsistencies away from the active material patch
            # do not determine the geometric outside. Reject active ambiguity.
            active_tendon_faces=sorted({row[-2]for row in rows if row[0]=='tendon' and row[-1]})
            if active_tendon_faces:
                corners=tendon.attributes['normals'][tendon.faces[active_tendon_faces]].astype(float);corners/=np.linalg.norm(corners,axis=2)[:,:,None];cos=np.einsum('fci,fi->fc',corners,normals['tendon'][active_tendon_faces]);state['tendon_active_normal_cosine_min']=float(cos.min());assert cos.min()>0
            atomic_npz(BASE/f'candidate-checkpoint-{iteration:02}.npz',vertices=candidate,faces=f,frame=np.array(frame),source_vertices=v,rings=rings,source_profiles=profiles,source_centres=centers,station_displacements=d,station_rotations=rot,fixed_station_indices=fixed,source_to_exact_welded=z['source_to_exact_welded'])
            atomic_json(BASE/'checkpoint.json',report,overwrite=(BASE/'checkpoint.json').exists())
            corridor=plan['inferior_corridor'];core=np.asarray(corridor['stations']);caps=np.asarray(corridor['maximum_center_y_metres']);corridor_violation=float(np.max(centers[core,1]+d[core,1]-caps));state['inferior_corridor_max_violation_m']=corridor_violation
            all_clear=all_clear and corridor_violation<=1e-10
            if all_clear:report['passed']=True;break
            if iteration==plan['maximum_linearizations']:break
            W=np.zeros((len(rows),len(rings)));ns=[];bounds=[]
            for i,(name,j,c,n,point,face,oriented)in enumerate(rows):
                W[i,j]=1;ns.append(n);bounds.append(plan['clearance_metres']+padding-n@(centers[j]+rot[j]@offset[j,c]-point))
            CW=np.zeros((len(core),len(rings)));CW[np.arange(len(core)),core]=1;W=np.vstack([W,CW]);ns.extend([[0.,-1.,0.]]*len(core));bounds.extend((centers[core,1]-caps).tolist())
            displacement,solver=solve_weighted_section_socp(K,W,ns,bounds,maximum_motion_m=plan['maximum_motion_metres'],time_limit_seconds=plan['solver_seconds_per_round'],maximum_iterations=100);state['solver']=json_safe_solver_diagnostics(solver);state['constraint_rows']=len(rows)
            if displacement is None:report['failure']='Coupled local halfspaces infeasible or bounded solver did not complete';break
            d=displacement
        if report['passed']:
            quality=audit(mesh,require_closed=False,expected_components=1);assert quality['passed'];protected=rings[fixed].ravel();assert np.array_equal(candidate[protected],v[protected]);before=profiles-profiles.mean(1)[:,None];after=candidate[rings]-candidate[rings].mean(1)[:,None];replay=np.einsum('sji,svj->svi',rot,after);area=lambda a:np.cross(a,np.roll(a,-1,axis=1)).sum(1)/2;oldarea=abs(np.einsum('si,si->s',area(before),old_tangent));newarea=abs(np.einsum('si,si->s',area(after),nt));relative=abs(newarea-oldarea)/oldarea;motion=float(np.linalg.norm(candidate-v,axis=1).max());assert motion<=plan['maximum_motion_metres']and relative.max()<.001
            report.update(whole_source_quality=quality,protected_controls_literal=True,maximum_motion_metres=motion,maximum_rigid_profile_replay_error_metres=float(np.linalg.norm(replay-before,axis=2).max()),maximum_projected_axial_area_change_relative=float(relative.max()),source_indices_literal=True)
            atomic_npz(BASE/'candidate.npz',vertices=candidate,faces=f,frame=np.array(frame),source_vertices=v,rings=rings,source_profiles=profiles,source_centres=centers,station_displacements=d,station_rotations=rot,fixed_station_indices=fixed,source_to_exact_welded=z['source_to_exact_welded'])
    except Exception as error:
        import traceback
        report['passed']=False;report['failure']=dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc())
    report['passed']=bool(report['passed']);report['elapsed_seconds']=time.monotonic()-began;atomic_json(BASE/'proof.json',report);print(json.dumps(dict(passed=report['passed'],elapsed_seconds=report['elapsed_seconds'],failure=report.get('failure'),maximum_motion_metres=report.get('maximum_motion_metres'),iterations=[dict(iteration=s['iteration'],obstacles=[dict(name=o['name'],active=o['active_controls'],inside=o['inside_controls'],crossings=o['contacts']['crossing_face_pairs'],minimum_signed=o['minimum_local_signed_distance_metres'])for o in s['obstacles']],solver=s.get('solver'))for s in records])))

if __name__=='__main__':main()
