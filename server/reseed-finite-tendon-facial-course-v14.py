"""Stage 2: finite-sheet detour seed, then free C2 control optimisation."""
from pathlib import Path
from copy import deepcopy
import sys,json,time,runpy,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';out=base/'male-left-facial-cervical-routes-finite-tendon-seed-v14';assert not out.exists();out.mkdir()
    prior=base/'male-left-facial-cervical-routes-free-c2-v10/routes.json';plan=deepcopy(json.loads(prior.read_bytes()));source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256'];assert digest_file(source)==expected
    solved=base/'male-left-facial-cervical-routes-targeted-c2-v13';report=json.loads((solved/'stage1-evidence.json').read_bytes());assert not report['all_free_c2_curvature_and_profile_gates_pass']
    tonsil=next(r for r in report['routes']if r['name']=='tonsillar');assert tonsil['continuous_certificate']['continuous_curvature_bound_pass']and tonsil['finite_prisms']['finite_wall_cells_pass']
    for route in plan['routes']:
        if route['name']=='tonsillar':
            file=solved/'tonsillar-guide.npz';z=np.load(file);route.update(guide=str(file.relative_to(ROOT)),guide_sha256=digest_file(file),stations=len(z['points']),length_m=float(np.linalg.norm(np.diff(z['points'],axis=0),axis=1).sum()))
    route=next(r for r in plan['routes']if r['name']=='glandular-1');old=np.load(ROOT/route['guide']);baseline=old['points'];parent=Mesh.load(ROOT/plan['fixed_parent']['mesh'])
    helper=load('facial_occupancy_v14','server/check-facial-cervical-feasibility.py');search=load('facial_domain_v14','server/plan-native-reference-guided-facial-routes.py')
    pc,_=helper.occupancy_classifier(parent,True);radius=.000400021+.00010;step=.0006;margin=.008
    axes=[np.arange(baseline[:,k].min()-margin,baseline[:,k].max()+margin+step/2,step)for k in range(3)];shape=tuple(map(len,axes));q=np.stack(np.meshgrid(*axes,indexing='ij'),-1).reshape(-1,3)
    domain=cKDTree(baseline).query(q)[0]<=margin;points=q[domain];selected=np.flatnonzero(domain);safe=np.ones(len(points),bool)
    parent_distance=np.sqrt(igl.point_mesh_squared_distance(points,parent.vertices,parent.faces)[0]);parent_core=(abs(igl.winding_number(pc.vertices,pc.faces,points))>.99)&(parent_distance>radius)
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();binding=json.loads((solved/'literal-source-rebase.json').read_bytes());blockers=[]
    for row in binding['fixed_blocker_geometry']:
        if time.monotonic()-began>28:raise RuntimeError('28-second detour-query work budget reached')
        pid=row['part_id'];mesh=reader.read_part(atlas,parts[pid],frame=plan['frame'],catalogue_path=source,catalogue_sha=expected);assert mesh.geometry_sha==row['geometry_sha256']
        take=np.flatnonzero(safe&np.all(points>=mesh.bounds[0]-radius,axis=1)&np.all(points<=mesh.bounds[1]+radius,axis=1))
        if not len(take):continue
        p=points[take];distance=np.sqrt(igl.point_mesh_squared_distance(p,mesh.vertices,mesh.faces)[0]);bad=distance<=radius
        classifier,mode=helper.occupancy_classifier(mesh,parts[pid].get('system')in{'arterial','venous'})
        if classifier is not None:bad|=abs(igl.winding_number(classifier.vertices,classifier.faces,p))>.99
        bad&=~parent_core[take];safe[take[bad]]=False;blockers.append(dict(part_id=pid,geometry_sha256=mesh.geometry_sha,queried=len(take),rejected=int(bad.sum()),occupancy_mode=mode))
    mask=np.zeros(len(q),bool);mask[selected]=safe
    path,proof=search.lattice_path(mask,baseline[0],baseline[-1],shape,q,step,began)
    result=dict(stage=2,maximum_stages=3,source_catalogue_sha256=expected,source_calibre_unchanged=True,finite_sheet_used_as_finite_triangles=True,infinite_plane_surrogate=False,
                proposal_points=len(points),safe_points=int(safe.sum()),fixed_boundaries=blockers,search=proof,installed=False,blender_queued=False)
    if path is None:
        result.update(stopped_before_C2=True,elapsed_seconds=time.monotonic()-began);atomic_json(out/'stage2-evidence.json',result);print(json.dumps(result));return
    path=path[np.r_[True,np.linalg.norm(np.diff(path,axis=0),axis=1)>1e-10]]
    file=out/'glandular-1-detour-guide.npz';atomic_npz(file,points=path,outer_radii=np.full(len(path),route['illustrative_radius_m']),inner_radii=np.full(len(path),route['illustrative_radius_m']-route['illustrative_wall_m']),frame=np.array(plan['frame']))
    route.update(guide=str(file.relative_to(ROOT)),guide_sha256=digest_file(file),stations=len(path),length_m=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()))
    plan.update(source_inputs=str((solved/'inputs.json').relative_to(ROOT)),source_inputs_sha256=digest_file(solved/'inputs.json'),construction_contract=str((solved/'construction-contract.json').relative_to(ROOT)),construction_contract_sha256=digest_file(solved/'construction-contract.json'))
    atomic_json(out/'routes.json',plan);result.update(proposal_length_m=route['length_m'],no_intermediate_detour_points_pinned=True,detour_seconds=time.monotonic()-began);atomic_json(out/'stage2-evidence.json',result)
    assert time.monotonic()-began<28,'No new C2 solve starts after the declared detour budget'
    sys.argv=['optimise-free-c2-facial-course-v9.py','--targeted-blockers','--reseeded-plan',str((out/'routes.json').relative_to(ROOT))]
    runpy.run_path(str(ROOT/'server/optimise-free-c2-facial-course-v9.py'),run_name='__main__')
    result.update(stage2_total_seconds=time.monotonic()-began);atomic_json(out/'stage2-evidence.json',result)


if __name__=='__main__':main()
