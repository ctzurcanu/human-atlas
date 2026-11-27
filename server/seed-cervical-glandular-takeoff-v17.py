"""One native-source takeoff proposal; previous passing guides remain literal."""
from pathlib import Path
from copy import deepcopy
import sys,json,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';out=base/'male-left-facial-cervical-routes-cervical-seed-v17';assert not out.exists();out.mkdir()
    takeoff=base/'male-left-facial-cervical-takeoff-v16';witness=json.loads((takeoff/'takeoff-evidence.json').read_bytes());contract=json.loads((takeoff/'construction-contract.json').read_bytes())
    plan=deepcopy(json.loads((base/'male-left-facial-cervical-routes-free-c2-v10/routes.json').read_bytes()));source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256'];assert digest_file(source)==expected==witness['source_catalogue_sha256']
    solved=base/'male-left-facial-cervical-routes-targeted-c2-v13'
    tonsillar=next(r for r in plan['routes']if r['name']=='tonsillar');file=solved/'tonsillar-guide.npz';z=np.load(file);tonsillar.update(guide=str(file.relative_to(ROOT)),guide_sha256=digest_file(file),stations=len(z['points']),length_m=float(np.linalg.norm(np.diff(z['points'],axis=0),axis=1).sum()))
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);route=next(r for r in plan['routes']if r['name']=='glandular-1');start=np.array(witness['selected']['center_m']);stop=np.array(route['recipient']['point_m'])
    route.update(parent_station=witness['selected']['station'],source_profile_template=witness['source_profile_template'],source_profile_template_original_station=48)
    # This straight anatomical-domain seed is not a required course or a pinned intermediate.
    baseline=np.linspace(start,stop,40);helper=load('cervical_occupancy_v17','server/check-facial-cervical-feasibility.py');search=load('cervical_search_v17','server/plan-native-reference-guided-facial-routes.py')
    fixedbranch=base/'male-left-facial-cervical-finite-network-free-c2-v10/glandular-0-fields.npz';bank=np.load(fixedbranch);sheet=Mesh(bank['fields'][0],bank['faces'],plan['frame']);gland0,mode=helper.occupancy_classifier(sheet,True);assert gland0 is not None
    file=out/'literal-glandular-0-outer-classifier.npz';gland0.save(file);extra=dict(part_id='CANDIDATE:glandular-0:passed-fixed-envelope',mesh=str(file.relative_to(ROOT)),mesh_sha256=digest_file(file),geometry_sha256=gland0.geometry_sha,source_fields=str(fixedbranch.relative_to(ROOT)),source_fields_sha256=digest_file(fixedbranch),qualification='Temporary closed occupancy classifier of literal passing glandular-0 outer carrier; not rendered tissue.')
    pc,_=helper.occupancy_classifier(parent,True);radius=route['illustrative_radius_m']*np.linalg.norm(np.array(route['source_profile_template']),axis=1).max()+.00010;step=.0006;margin=.006
    axes=[np.arange(baseline[:,k].min()-margin,baseline[:,k].max()+margin+step/2,step)for k in range(3)];shape=tuple(map(len,axes));q=np.stack(np.meshgrid(*axes,indexing='ij'),-1).reshape(-1,3)
    domain=cKDTree(baseline).query(q)[0]<=margin;points=q[domain];selected=np.flatnonzero(domain);safe=np.ones(len(points),bool)
    parent_distance=np.sqrt(igl.point_mesh_squared_distance(points,parent.vertices,parent.faces)[0]);parent_core=(abs(igl.winding_number(pc.vertices,pc.faces,points))>.99)&(parent_distance>radius)
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();binding=json.loads((solved/'literal-source-rebase.json').read_bytes());blockers=[]
    meshes=[]
    for row in binding['fixed_blocker_geometry']:
        pid=row['part_id'];mesh=reader.read_part(atlas,parts[pid],frame=plan['frame'],catalogue_path=source,catalogue_sha=expected);assert mesh.geometry_sha==row['geometry_sha256'];classifier,mode=helper.occupancy_classifier(mesh,parts[pid].get('system')in{'arterial','venous'});meshes.append((pid,mesh,classifier,mode))
    meshes.append((extra['part_id'],gland0,gland0,'closed-passing-candidate-envelope'))
    for pid,mesh,classifier,mode in meshes:
        if time.monotonic()-began>40:raise RuntimeError('40-second native takeoff-domain budget reached')
        take=np.flatnonzero(safe&np.all(points>=mesh.bounds[0]-radius,axis=1)&np.all(points<=mesh.bounds[1]+radius,axis=1))
        if not len(take):continue
        p=points[take];distance=np.sqrt(igl.point_mesh_squared_distance(p,mesh.vertices,mesh.faces)[0]);bad=distance<=radius
        if classifier is not None:bad|=abs(igl.winding_number(classifier.vertices,classifier.faces,p))>.99
        bad&=~parent_core[take];safe[take[bad]]=False;blockers.append(dict(part_id=pid,geometry_sha256=mesh.geometry_sha,queried=len(take),rejected=int(bad.sum()),occupancy_mode=mode))
    mask=np.zeros(len(q),bool);mask[selected]=safe;path,proof=search.lattice_path(mask,start,stop,shape,q,step,began)
    result=dict(schema='cervical-native-takeoff-seed/v17',stage='course-domain',source_catalogue_sha256=expected,selected_native_parent_station=route['parent_station'],source_profile_template_original_station=48,source_calibre_unchanged=True,fixed_parent_and_recipient_literal=True,proposal_points=len(points),safe_points=int(safe.sum()),fixed_boundaries=blockers,search=proof,installed=False,blender_queued=False)
    if path is None:
        result.update(stopped_before_C2=True,elapsed_seconds=time.monotonic()-began);atomic_json(out/'stage-evidence.json',result);print(json.dumps(result));return
    path=path[np.r_[True,np.linalg.norm(np.diff(path,axis=0),axis=1)>1e-10]];file=out/'glandular-1-seed-guide.npz';atomic_npz(file,points=path,outer_radii=np.full(len(path),route['illustrative_radius_m']),inner_radii=np.full(len(path),route['illustrative_radius_m']-route['illustrative_wall_m']),frame=np.array(plan['frame']))
    route.update(guide=str(file.relative_to(ROOT)),guide_sha256=digest_file(file),stations=len(path),length_m=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()))
    contract.update(corrected_routes_only=['glandular-1'],parent_departure_policy='Native cervical interval 26–34; selected 34; original parent vertices unchanged. Original station-48 section profile transported without calibre reduction.',native_authored_controls={**contract['native_authored_controls'],'glandular-1':dict(parent_station=34,controls_m=path.tolist(),intermediates_are_initialisation_only=True)},extra_fixed_meshes=[extra],new_geometry_constructed=False)
    atomic_json(out/'construction-contract.json',contract);plan.update(extra_fixed_meshes=[extra],source_inputs=str((solved/'inputs.json').relative_to(ROOT)),source_inputs_sha256=digest_file(solved/'inputs.json'),construction_contract=str((out/'construction-contract.json').relative_to(ROOT)),construction_contract_sha256=digest_file(out/'construction-contract.json'));atomic_json(out/'routes.json',plan)
    result.update(seed_length_m=route['length_m'],seed_stations=len(path),no_intermediate_control_pinned=True,elapsed_seconds=time.monotonic()-began);atomic_json(out/'stage-evidence.json',result);assert digest_file(source)==expected;print(json.dumps({k:v for k,v in result.items()if k!='fixed_boundaries'}))


if __name__=='__main__':main()
