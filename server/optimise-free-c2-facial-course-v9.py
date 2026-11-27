"""Free intermediate C2 control optimisation in fixed native facial anatomy.

At most five declared sequential-convex iterations per route, within one
50-second proposal stage. No voxel intermediate is treated as a fixed point.
"""
from pathlib import Path
from copy import deepcopy
import sys,json,time,importlib.util,argparse
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.algorithms.tube_sections import arclength,rotation_minimizing_frames,bend_diagnostics
from geometry_tools.algorithms.tube_continuity import loft_faces,finite_wall_cells
from geometry_tools.algorithms.tube_c2_course import C2Course

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'
OUT=BASE/'male-left-facial-cervical-routes-free-c2-v9'


def main():
    global OUT,EXPECTED
    parser=argparse.ArgumentParser();parser.add_argument('--correct-open-sheets',action='store_true');parser.add_argument('--targeted-blockers',action='store_true');parser.add_argument('--reseeded-plan');parser.add_argument('--output-version');parser.add_argument('--expected-source');parser.add_argument('--registry-context-inventory');parser.add_argument('--only-route');parser.add_argument('--active-trust-mm',type=float);args=parser.parse_args()
    if args.active_trust_mm is not None:assert 0<args.active_trust_mm<=1. and args.registry_context_inventory and args.only_route
    if args.expected_source:
        assert len(args.expected_source)==64 and args.registry_context_inventory and args.reseeded_plan
        EXPECTED=args.expected_source
    current_stage=1 if args.targeted_blockers else(2 if args.correct_open_sheets else 1)
    if args.correct_open_sheets:OUT=BASE/'male-left-facial-cervical-routes-free-c2-v10'
    if args.targeted_blockers:OUT=BASE/'male-left-facial-cervical-routes-targeted-c2-v13';args.correct_open_sheets=True
    if args.reseeded_plan:OUT=BASE/'male-left-facial-cervical-routes-targeted-c2-v14';current_stage=2
    if args.output_version:
        assert '/'not in args.output_version and '..'not in args.output_version
        OUT=BASE/('male-left-facial-cervical-routes-'+args.output_version)
    began=time.monotonic();assert not OUT.exists();OUT.mkdir()
    oldfile=ROOT/args.reseeded_plan if args.reseeded_plan else BASE/('male-left-facial-cervical-routes-free-c2-v10/routes.json'if args.targeted_blockers else'male-left-facial-cervical-routes-native-domain-v7/routes.json');plan=deepcopy(json.loads(oldfile.read_bytes()));source=ROOT/plan['source_catalogue'];assert digest_file(source)==EXPECTED
    oldsha=plan['source_catalogue_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    spec=importlib.util.spec_from_file_location('facial_occupancy',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    inputs=deepcopy(json.loads((ROOT/plan['source_inputs']).read_bytes()));literal=[]
    for row in inputs['fixed_sources']:
        if row.get('missing_native_source'):continue
        assert row['source_record']==parts[row['part_id']],('Native facial source record changed',row['part_id'])
        mesh=reader.read_part(atlas,parts[row['part_id']],frame=plan['frame'],catalogue_path=source,catalogue_sha=EXPECTED)
        assert mesh.geometry_sha==row['geometry_sha256'],('Native facial source geometry changed',row['part_id'])
        literal.append(dict(part_id=row['part_id'],geometry_sha256=mesh.geometry_sha,source_record_literal=True))
    v7=json.loads((BASE/'male-left-facial-cervical-routes-native-domain-v7/stage1-evidence.json').read_bytes());fixed={}
    for route in v7['routes']:
        for row in route.get('fixed_boundaries',[]):
            pid=row['part_id']
            if pid in fixed:continue
            mesh=reader.read_part(atlas,parts[pid],frame=plan['frame'],catalogue_path=source,catalogue_sha=EXPECTED)
            assert mesh.geometry_sha==row['geometry_sha256'],('Native blocker changed',pid)
            classifier,mode=helper.occupancy_classifier(mesh,parts[pid].get('system')in{'arterial','venous'})
            fixed[pid]=(mesh,classifier,mode)
    for row in plan.get('extra_fixed_meshes',[]):
        file=ROOT/row['mesh'];assert digest_file(file)==row['mesh_sha256']
        mesh=Mesh.load(file);assert mesh.geometry_sha==row['geometry_sha256']
        classifier,mode=helper.occupancy_classifier(mesh,True);assert classifier is not None
        fixed[row['part_id']]=(mesh,classifier,mode)
    if args.targeted_blockers:
        full=json.loads((BASE/'male-left-facial-cervical-finite-network-free-c2-v10/feasibility-v2/evidence.json').read_bytes())
        if full['source_catalogue_sha256']!=EXPECTED:
            proof=ROOT/plan['literal_context_rebase_proof'];assert digest_file(proof)==plan['literal_context_rebase_proof_sha256'];rebase=json.loads(proof.read_bytes());assert rebase['previous_catalogue_sha256']==full['source_catalogue_sha256']and rebase['current_catalogue_sha256']==EXPECTED
        for pid in ['ZA:Hypoglossal nerve (XII).l','ZA:Platysma.l:Tendon']:
            row=next(r for r in full['context']if r['part_id']==pid)
            mesh=reader.read_part(atlas,parts[pid],frame=plan['frame'],catalogue_path=source,catalogue_sha=EXPECTED)
            assert mesh.geometry_sha==row['target_geometry_sha256']
            classifier,mode=helper.occupancy_classifier(mesh,parts[pid].get('system')in{'arterial','venous'})
            fixed[pid]=(mesh,classifier,mode)
    if args.registry_context_inventory:
        inventory=json.loads((ROOT/args.registry_context_inventory).read_bytes());assert inventory['source_catalogue_sha256']==EXPECTED and inventory['registry_revision']==plan['context_registry_revision']
        for row in inventory['all_registered_overlap_candidates']:
            if not row['source_category'].startswith('native-')or row['hidden']or row['suppressed']:continue
            pid=row['part_id']
            if pid=='ZA:Facial artery.l':continue # Declared parent junction, separately audited by the finite union.
            assert row['catalogue_record_matches_registry']and row['frame']==plan['frame']
            if pid in fixed:assert fixed[pid][0].geometry_sha==row['geometry_sha256'];continue
            mesh=reader.read_part(atlas,parts[pid],frame=plan['frame'],catalogue_path=source,catalogue_sha=EXPECTED);assert mesh.geometry_sha==row['geometry_sha256'];classifier,mode=helper.occupancy_classifier(mesh,parts[pid].get('system')in{'arterial','venous'});fixed[pid]=(mesh,classifier,mode)
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);parent_classifier,_=helper.occupancy_classifier(parent,True);origin=parent.bounds.mean(0)
    rings=np.array(plan['fixed_parent']['source_station_components'][0]['rings']);centers=parent.vertices[rings].mean(1)
    binding=dict(previous_source_sha256=oldsha,current_source_sha256=EXPECTED,fixed_source_records=literal,
                 fixed_blocker_geometry=[dict(part_id=k,geometry_sha256=v[0].geometry_sha)for k,v in fixed.items()],
                 all_fixed_sources_literal=True,no_native_geometry_moved=True)
    atomic_json(OUT/'literal-source-rebase.json',binding)
    inputs['source_catalogue_sha256']=EXPECTED;atomic_json(OUT/'inputs.json',inputs);plan.update(source_catalogue_sha256=EXPECTED,source_inputs=str((OUT/'inputs.json').relative_to(ROOT)),source_inputs_sha256=digest_file(OUT/'inputs.json'))
    contract=deepcopy(json.loads((ROOT/plan['construction_contract']).read_bytes()));contract.update(schema='free-c2-native-facial-construction/v9',source_catalogue_sha256=EXPECTED,
        literal_rebase_proof=str((OUT/'literal-source-rebase.json').relative_to(ROOT)),literal_rebase_proof_sha256=digest_file(OUT/'literal-source-rebase.json'),
        fixed_anatomical_constraints=['Native facial-parent coordinates/source profiles','Actual recipient endpoint and terminal calibre','Native skeleton, muscles, veins, mucosa and ducts'],
        intermediate_controls='ALL free cubic-spline controls; none of the voxel points, previous SG/SP point witnesses or fillet intermediates is pinned.',
        continuous_curvature_proof='C2 cubic derivative Bezier speed-projection lower bounds and acceleration endpoint second-order-cone norm bounds over whole intervals.',
        sequential_convex_iteration_limit=5,stage_internal_seconds=50,maximum_new_saved_stages=3,source_calibres_not_reduced=True,installed=False)
    contract['open_sheet_proposal_method']='Actual unsigned closest-triangle distance; an open sheet is not promoted to an infinite occupied halfspace. Final actual crossings and closed-source occupied volume remain mandatory.'if args.correct_open_sheets else 'Signed local source-sheet plane proposal'
    if args.targeted_blockers:contract.update(schema='targeted-native-facial-c2-contract/v13',corrected_routes_only=['tonsillar','glandular-1'],literal_passed_routes=['ascending-palatine','glandular-0'],new_explicit_fixed_blockers=['ZA:Hypoglossal nerve (XII).l','ZA:Platysma.l:Tendon'])
    if args.reseeded_plan:contract.update(schema='native-finite-tendon-detour-c2/v14',corrected_routes_only=['glandular-1'],literal_passed_routes=['ascending-palatine','tonsillar','glandular-0'],finite_tendon_detour_input=str(oldfile.relative_to(ROOT)),finite_tendon_detour_sha256=digest_file(oldfile))
    if args.only_route:contract.update(schema='current-registered-targeted-facial-c2/v27',corrected_routes_only=[args.only_route],literal_passed_routes=[r['name']for r in plan['routes']if r['name']!=args.only_route],registry_context_inventory=args.registry_context_inventory,registry_context_inventory_sha256=digest_file(ROOT/args.registry_context_inventory),all_registered_native_blockers_included=True,unregistered_frame_or_ownership_sources_remain_unknown=True)
    if args.active_trust_mm is not None:contract.update(schema='conservative-active-set-facial-c2/v28',local_control_step_radius_mm=args.active_trust_mm,omission_proof='At each sampled spline parameter, motion is bounded by the sum of absolute free-control spline weights times the enforced per-control step radius. Exact signed/unsigned distance is 1-Lipschitz; constraints farther than this bound plus full profile radius and clearance are safely inactive for that step. AABB lower bounds safely cull more distant source queries. All physical sources remain in each exhaustive sample check; full emitted triangle/volume checks are separate.',all_fixed_object_count=len(fixed))
    report_path='reports/ta98-male/facial-c2-construction-contract-'+args.output_version+'.json'if args.output_version else('reports/ta98-male/facial-targeted-c2-construction-contract-v14.json'if args.reseeded_plan else('reports/ta98-male/facial-targeted-c2-construction-contract-v13.json'if args.targeted_blockers else('reports/ta98-male/facial-free-c2-construction-contract-v10.json'if args.correct_open_sheets else 'reports/ta98-male/facial-free-c2-construction-contract-v9.json')))
    atomic_json(OUT/'construction-contract.json',contract);atomic_json(ROOT/report_path,contract)
    reports=[];passed=True
    for route in plan['routes']:
        if (args.only_route and route['name']!=args.only_route)or(not args.only_route and args.targeted_blockers and route['name']in({'ascending-palatine','tonsillar','glandular-0'}if args.reseeded_plan else{'ascending-palatine','glandular-0'})):
            assert digest_file(ROOT/route['guide'])==route['guide_sha256']
            reports.append(dict(name=route['name'],passing_guide_retained_literal=True,source_guide_sha256=route['guide_sha256'],iterations=[]));continue
        old=np.load(ROOT/route['guide']);oldpoints=old['points'];arc=arclength(oldpoints)*1000;length=arc[-1]
        count=max(8,min(50,int(np.ceil(length/2.5))+1));knots=np.linspace(0,length,count)
        reference=np.column_stack([np.interp(knots,arc,(oldpoints[:,d]-origin[d])*1000)for d in range(3)])
        delta=parent.vertices[rings[route['parent_station']]]-centers[route['parent_station']]
        _,_,vt=np.linalg.svd(delta,full_matrices=False);normal=vt[-1]
        if normal@(centers[route['parent_station']+1]-centers[route['parent_station']-1])<0:normal=-normal
        initial=delta[0]-normal*(delta[0]@normal);initial/=np.linalg.norm(initial);basis=np.c_[initial,np.cross(normal,initial),normal]
        profile=delta@basis;profile[:,2]=0;profile/=np.linalg.norm(profile,axis=1).mean()
        if 'source_profile_template'in route:
            profile=np.asarray(route['source_profile_template'],float);assert profile.shape==delta.shape and np.all(profile[:,2]==0)
        outer_max=route['illustrative_radius_m']*np.linalg.norm(profile,axis=1).max()
        minimum_radius=max(.00073336,outer_max/.75);course=C2Course(reference,length,count,minimum_radius*1000,speed_lower=.85,trust_radius_mm=8.)
        coefficient=reference.copy();parameters=np.linspace(0,length,max(40,int(np.ceil(length/.25))+1));oldcheck=np.column_stack([np.interp(parameters,arc,oldpoints[:,d])for d in range(3)])
        report=dict(name=route['name'],free_intermediate_controls=count-2,fixed_endpoint_controls=2,source_outer_profile_maximum_radius_m=float(outer_max),minimum_bend_radius_m=float(minimum_radius),iterations=[]);reports.append(report)
        displacement=(np.abs(course.basis(parameters)[:,1:-1]).sum(1)*args.active_trust_mm/1000)if args.active_trust_mm is not None else None
        def measure(q,active_step=False):
            movement=displacement if active_step and displacement is not None else np.zeros(len(q))
            pd=np.sqrt(igl.point_mesh_squared_distance(q,parent.vertices,parent.faces)[0]);pi=abs(igl.winding_number(parent_classifier.vertices,parent_classifier.faces,q))>.99
            parent_core=pi&(pd>outer_max+.00007+movement);constraints=[];bad=[];minimum=np.inf;activation=[]
            for pid,(mesh,classifier,mode)in fixed.items():
                intended={'ascending-palatine':{'ZA:Palatopharyngeus muscle.l','ZA:Soft palate'},'tonsillar':{'ZA:Palatine tonsil.l'},'glandular-0':{'ZA:Submandibular gland.l'},'glandular-1':{'ZA:Submandibular gland.l'}}
                if pid in intended[route['name']]:
                    activation.append(dict(part_id=pid,reason='Declared legitimate recipient tissue interface',samples=len(q)));continue
                if args.active_trust_mm is not None:
                    box_delta=np.maximum(mesh.bounds[0]-q,0)+np.maximum(q-mesh.bounds[1],0);box_distance=np.linalg.norm(box_delta,axis=1);possible=box_distance<=outer_max+.00007+movement;selected=np.flatnonzero(possible&~parent_core)
                    active_row=dict(part_id=pid,samples=len(q),safe_parent_core_samples=int(parent_core.sum()),safe_AABB_distance_omissions=int((~possible&~parent_core).sum()),reason='Actual source AABB distance lower bound exceeds full profile radius, clearance and rigorous local spline motion bound.',exact_distance_omissions=0,active_constraints=0);activation.append(active_row)
                else:selected=np.flatnonzero(np.all(q>=mesh.bounds[0]-.003,axis=1)&np.all(q<=mesh.bounds[1]+.003,axis=1)&~parent_core)
                if route['name']=='tonsillar'and pid=='ZA:Superior pharyngeal constrictor.l':
                    selected=selected[np.linalg.norm(q[selected]-np.array(route['recipient']['point_m']),axis=1)>.010]
                if not len(selected):continue
                d2,faces,near=igl.point_mesh_squared_distance(q[selected],mesh.vertices,mesh.faces);distance=np.sqrt(d2);vector=q[selected]-near
                n=vector/np.maximum(distance[:,None],1e-18)
                signed=distance.copy()
                if classifier is not None:
                    inside=abs(igl.winding_number(classifier.vertices,classifier.faces,q[selected]))>.99;signed[inside]*=-1;n[inside]*=-1
                elif not args.correct_open_sheets:
                    tri=mesh.triangles[faces];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
                    sign=np.sign(np.einsum('ij,ij->i',oldcheck[selected]-near,normal));sign[sign==0]=1
                    n=normal*sign[:,None];signed=np.einsum('ij,ij->i',vector,n)
                minimum=min(minimum,float(signed.min()-outer_max))
                reject=signed<outer_max+.00005
                if reject.any():bad.append(dict(part_id=pid,sample_count=int(reject.sum()),minimum_source_profile_clearance_m=float(signed.min()-outer_max)))
                active=signed<outer_max+.00007+movement[selected]if args.active_trust_mm is not None else np.ones(len(selected),bool)
                if args.active_trust_mm is not None:
                    active_row.update(exact_distance_omissions=int((~active).sum()),active_constraints=int(active.sum()),minimum_omitted_signed_distance_margin_m=float((signed[~active]-outer_max-.00007-movement[selected[~active]]).min())if(~active).any()else None,exact_distance_omission_reason='Exact signed/unsigned source distance exceeds full profile radius + clearance + complete local spline-step displacement bound; 1-Lipschitz distance certificate.')
                for j,i in enumerate(selected):
                    if not active[j]:continue
                    rhs=(outer_max+.00007-signed[j])*1000+n[j]@((q[i]-origin)*1000)
                    constraints.append((float(parameters[i]),n[j],float(rhs),dict(part_id=pid,sample=int(i),parameter_mm=float(parameters[i]),actual_query_point_m=q[i].tolist(),source_face=int(faces[j]),source_surface_point_m=near[j].tolist(),signed_or_unsigned_distance_m=float(signed[j]),closed_occupancy_mode=mode)))
            return constraints,bad,minimum,activation
        for iteration in range(5):
            if time.monotonic()-began>44:
                report['failure']='44-second optimisation-work budget reached';passed=False;break
            q=CubicSpline(knots,coefficient,axis=0,bc_type='natural')(parameters)/1000+origin
            constraints,bad,clearance,activation=measure(q,active_step=args.active_trust_mm is not None)
            if args.active_trust_mm is not None:
                atomic_npz(OUT/(route['name']+'-active-step-'+str(iteration)+'.npz'),control_center_mm=coefficient,query_points_m=q,parameter_mm=parameters,spline_motion_bound_m=displacement,frame=np.array(plan['frame']))
                atomic_json(OUT/(route['name']+'-active-step-'+str(iteration)+'.json'),dict(fixed_objects=len(fixed),local_control_step_radius_mm=args.active_trust_mm,source_activations=activation,active_native_constraints=len(constraints),previous_violations=bad,all_objects_retained=True,active_geometry_surrogates_are_local_only=True))
            proposed,solve=course.solve(constraints,time_limit=min(5.,max(.25,48-(time.monotonic()-began))),**dict(step_center=coefficient,step_radius_mm=args.active_trust_mm)if args.active_trust_mm is not None else{})
            report['iterations'].append(dict(iteration=iteration,linearised_native_constraints=len(constraints),previous_violations=bad,solver=solve))
            if proposed is None:
                report['failure']='Continuous C2 radius/speed/native-clearance cone constraints are infeasible or the fixed solve budget expired';passed=False;break
            if args.active_trust_mm is not None:
                actual_step=float(np.linalg.norm(proposed-coefficient,axis=1).max());assert actual_step<=args.active_trust_mm+1e-6
                atomic_npz(OUT/(route['name']+'-candidate-step-'+str(iteration)+'.npz'),controls_mm=proposed,frame=np.array(plan['frame']));report['iterations'][-1]['maximum_actual_control_step_mm']=actual_step
            coefficient=proposed;q=CubicSpline(knots,coefficient,axis=0,bc_type='natural')(parameters)/1000+origin
            _,bad,clearance,post_activation=measure(q);report['iterations'][-1].update(actual_sampled_violations=bad,minimum_sampled_profile_clearance_m=clearance,exhaustive_fixed_object_count=len(fixed))
            if not bad and solve['certificate']['continuous_curvature_bound_pass']:break
        else:
            passed=False;report['failure']='Five declared sequential-convex iterations did not clear fixed native surfaces'
        if not passed:break
        points_mm,tangents,t=course.sample(coefficient,chord_error_mm=.003,max_step_mm=.25);path=points_mm/1000+origin;path[[0,-1]]=oldpoints[[0,-1]]
        frames=rotation_minimizing_frames(path,tangents,initial_reference=initial);section=np.einsum('kij,nj->kni',frames,profile)
        field=np.array([path[:,None]+r*section for r in [route['illustrative_radius_m'],route['illustrative_radius_m']-route['illustrative_wall_m']]])
        faces=loft_faces(len(path),len(profile));prisms=finite_wall_cells(field.reshape(2,-1,3),faces);bends=bend_diagnostics(path,field,tangents,.85)
        certificate=course.certificate(coefficient);report.update(continuous_certificate=certificate,finite_prisms=prisms,maximum_curvature_support=float(bends['curvature_support'].max()),curve_stations=len(path),source_endpoints_literal=bool(np.array_equal(path[[0,-1]],oldpoints[[0,-1]])))
        atomic_npz(OUT/(route['name']+'-free-c2-controls.npz'),control_points_mm=coefficient,physical_knots_mm=knots,origin_m=origin,frame=np.array(plan['frame']))
        atomic_npz(OUT/(route['name']+'-free-c2-fields.npz'),fields=field.reshape(2,-1,3),faces=faces,guide_points=path,analytic_tangents=tangents,frame=np.array(plan['frame']))
        if not certificate['continuous_curvature_bound_pass']or not prisms['finite_wall_cells_pass']or len(bends['failed_stations']):
            passed=False;report['failure']='Whole C2 curvature certificate or actual finite profile Jacobian failed';break
        file=OUT/(route['name']+'-guide.npz');atomic_npz(file,points=path,analytic_tangents=tangents,outer_radii=np.full(len(path),route['illustrative_radius_m']),inner_radii=np.full(len(path),route['illustrative_radius_m']-route['illustrative_wall_m']),frame=np.array(plan['frame']))
        route.update(guide=str(file.relative_to(ROOT)),guide_sha256=digest_file(file),stations=len(path),length_m=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()));route.pop('guide_retained_from_passing_native_context',None)
    result=dict(schema='targeted-free-c2-native-facial-courses/v13'if args.targeted_blockers else('free-c2-native-facial-courses/v10'if args.correct_open_sheets else 'free-c2-native-facial-courses/v9'),stage=current_stage,maximum_stages=3,source_catalogue_sha256=EXPECTED,all_free_c2_curvature_and_profile_gates_pass=passed,
                routes=reports,no_intermediate_anatomical_point_pinned=True,calibres_not_reduced=True,installed=False,blender_queued=False,elapsed_seconds=time.monotonic()-began)
    if passed:
        plan.update(construction_contract=str((OUT/'construction-contract.json').relative_to(ROOT)),construction_contract_sha256=digest_file(OUT/'construction-contract.json'),predecessor_route_contract=str(oldfile.relative_to(ROOT)),predecessor_route_contract_sha256=digest_file(oldfile));atomic_json(OUT/'routes.json',plan)
    assert digest_file(source)==EXPECTED;atomic_json(OUT/('stage2-evidence.json'if args.reseeded_plan else('stage1-evidence.json'if args.targeted_blockers else('stage2-evidence.json'if args.correct_open_sheets else 'stage1-evidence.json'))),result)
    print(json.dumps(dict(passed=passed,fixed_source_records_literal=len(literal),fixed_blockers_literal=len(fixed),routes=[dict(name=r['name'],free_controls=r.get('free_intermediate_controls'),retained_literal=r.get('passing_guide_retained_literal',False),failure=r.get('failure'),iterations=len(r['iterations']),certificate=r.get('continuous_certificate'),invalid_prisms=r.get('finite_prisms',{}).get('invalid_counts_by_interval'))for r in reports],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
