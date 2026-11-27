"""Single analytic radius-constrained trial of the verified native domain.

No radius search, source-calibre reduction, anatomical motion or silent retry.
The circular fillet radius is certified analytically; finite wall prisms and
the entire resulting geometry still need their independent context gates.
"""
from pathlib import Path
from copy import deepcopy
import sys,json,time,runpy,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_sections import filleted_path,rotation_minimizing_frames,path_tangents,bend_diagnostics
from geometry_tools.algorithms.tube_continuity import loft_faces,finite_wall_cells


def simplify(points,tolerance):
    """Bound every removed source waypoint against its replacement chord."""
    if len(points)<=2:return points
    a,b=points[[0,-1]];d=b-a;fraction=np.clip((points-a)@d/max(d@d,1e-30),0,1)
    distance=np.linalg.norm(points-a-fraction[:,None]*d,axis=1);index=int(distance.argmax())
    if distance[index]<=tolerance:return points[[0,-1]]
    return np.r_[simplify(points[:index+1],tolerance)[:-1],simplify(points[index:],tolerance)]


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion'
    predecessor=base/'male-left-facial-cervical-routes-native-domain-v7/routes.json'
    plan=deepcopy(json.loads(predecessor.read_bytes()));source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256'];assert digest_file(source)==expected
    out=base/'male-left-facial-cervical-routes-curvature-v8';assert not out.exists();out.mkdir()
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);rings=np.array(plan['fixed_parent']['source_station_components'][0]['rings']);centers=parent.vertices[rings].mean(1)
    contract=deepcopy(json.loads((ROOT/plan['construction_contract']).read_bytes()));contract.update(schema='analytic-curvature-facial-contract/v8',predecessor_route_sha256=digest_file(predecessor),
        analytic_strategy='One deterministic 0.20 mm bounded polyline simplification followed by exact tangent-continuous circular fillets.',
        simplification_tolerance_m=.00020,analytic_curvature_support_limit=.75,endpoint_positions_literal=True,source_profile_isometry=True,
        calibres_scaled=False,source_anatomy_modified=False,full_context_acceptance_pending=True,maximum_new_saved_stages=2,installed=False)
    atomic_json(out/'construction-contract.json',contract);atomic_json(ROOT/'reports/ta98-male/facial-curvature-construction-contract-v8.json',contract)
    reports=[]
    try:
        for route in plan['routes']:
            if route.get('guide_retained_from_passing_native_context'):
                reports.append(dict(name=route['name'],passing_predecessor_guide_retained_literal=True));continue
            old=np.load(ROOT/route['guide']);points=old['points'];controls=simplify(points,.00020)
            delta=parent.vertices[rings[route['parent_station']]]-centers[route['parent_station']]
            _,_,vt=np.linalg.svd(delta,full_matrices=False);normal=vt[-1]
            if normal@(centers[route['parent_station']+1]-centers[route['parent_station']-1])<0:normal=-normal
            reference=delta[0]-normal*(delta[0]@normal);reference/=np.linalg.norm(reference)
            basis=np.c_[reference,np.cross(normal,reference),normal];profile=delta@basis;profile[:,2]=0;profile/=np.linalg.norm(profile,axis=1).mean()
            maximum_outer_radius=route['illustrative_radius_m']*np.linalg.norm(profile,axis=1).max()
            minimum_bend_radius=maximum_outer_radius/.75
            report=dict(name=route['name'],source_guide_sha256=route['guide_sha256'],source_stations=len(points),simplified_control_points=len(controls),
                        outer_profile_maximum_radius_m=float(maximum_outer_radius),minimum_analytic_bend_radius_m=float(minimum_bend_radius),controls_m=controls.tolist())
            reports.append(report)
            path,tangents=filleted_path(controls,minimum_bend_radius,max_step_m=.00025,max_turn_rad=np.deg2rad(10))
            assert np.array_equal(path[[0,-1]],points[[0,-1]])
            frames=rotation_minimizing_frames(path,tangents,initial_reference=reference);section=np.einsum('kij,nj->kni',frames,profile)
            fields=np.array([path[:,None]+r*section for r in [route['illustrative_radius_m'],route['illustrative_radius_m']-route['illustrative_wall_m']]])
            faces=loft_faces(len(path),len(profile));prism=finite_wall_cells(fields.reshape(2,-1,3),faces)
            bends=bend_diagnostics(path,fields,tangents,.85)
            # The constructor uses measured polyline tangents, so verify that
            # exact replay too rather than accepting only analytic tangents.
            replay_frames=rotation_minimizing_frames(path,path_tangents(path),initial_reference=reference)
            replay_section=np.einsum('kij,nj->kni',replay_frames,profile)
            replay=np.array([path[:,None]+r*replay_section for r in [route['illustrative_radius_m'],route['illustrative_radius_m']-route['illustrative_wall_m']]])
            replay_prism=finite_wall_cells(replay.reshape(2,-1,3),faces)
            report.update(stations=len(path),analytic_fillet_prism=prism,constructor_replay_prism=replay_prism,
                          maximum_sampled_curvature_support=float(bends['curvature_support'].max()),analytic_radius_m=float(minimum_bend_radius),endpoints_literal=True)
            atomic_npz(out/(route['name']+'-analytic-fields.npz'),fields=fields.reshape(2,-1,3),faces=faces,guide_points=path,analytic_tangents=tangents,frame=np.array(plan['frame']))
            assert prism['finite_wall_cells_pass']and replay_prism['finite_wall_cells_pass']and not len(bends['failed_stations']), 'Analytic fillet or exact constructor replay failed finite profile constraints'
            file=out/(route['name']+'-guide.npz');atomic_npz(file,points=path,analytic_tangents=tangents,outer_radii=np.full(len(path),route['illustrative_radius_m']),inner_radii=np.full(len(path),route['illustrative_radius_m']-route['illustrative_wall_m']),frame=np.array(plan['frame']))
            route.update(guide=str(file.relative_to(ROOT)),guide_sha256=digest_file(file),stations=len(path),length_m=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()))
        plan.update(construction_contract=str((out/'construction-contract.json').relative_to(ROOT)),construction_contract_sha256=digest_file(out/'construction-contract.json'),predecessor_route_contract=str(predecessor.relative_to(ROOT)),predecessor_route_contract_sha256=digest_file(predecessor))
        atomic_json(out/'routes.json',plan)
        assert time.monotonic()-began<32,'32-second preconstruction budget exhausted'
        sys.argv=['construct-facial-cervical-network.py','--route-version','curvature-v8','--output-version','curvature-v8']
        runpy.run_path(str(ROOT/'server/construct-facial-cervical-network.py'),run_name='__main__')
        result=dict(stage=1,maximum_stages=2,source_catalogue_sha256=expected,all_analytic_and_finite_field_checks_pass=True,routes=reports,
                    finite_network_constructed=True,installed=False,elapsed_seconds=time.monotonic()-began)
    except Exception as error:
        result=dict(stage=1,maximum_stages=2,source_catalogue_sha256=expected,all_analytic_and_finite_field_checks_pass=False,routes=reports,
                    exact_failure_type=type(error).__name__,exact_failure=str(error),stopped_without_radius_or_source_changes=True,
                    finite_network_constructed=False,blender_queued=False,installed=False,elapsed_seconds=time.monotonic()-began)
        atomic_json(out/'stage1-evidence.json',result)
        print(json.dumps(dict(passed=False,failure=result['exact_failure'],routes=[{k:r[k]for k in ['name','source_stations','simplified_control_points','minimum_analytic_bend_radius_m']if k in r}for r in reports],elapsed_seconds=result['elapsed_seconds'])))
        return
    assert digest_file(source)==expected;atomic_json(out/'stage1-evidence.json',result)
    print(json.dumps(dict(passed=True,finite_network_constructed=True,elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
