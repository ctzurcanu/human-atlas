"""Source-bound cervical takeoff inventory/contract before new course work."""
from pathlib import Path
from copy import deepcopy
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';out=base/'male-left-facial-cervical-takeoff-v16';assert not out.exists();out.mkdir()
    plan=deepcopy(json.loads((base/'male-left-facial-cervical-routes-free-c2-v10/routes.json').read_bytes()));source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256'];assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);rings=np.asarray(plan['fixed_parent']['source_station_components'][0]['rings']);centers=parent.vertices[rings].mean(1)
    ids=['ZA:Submandibular gland.l','ZA:Mandible:Bone','ZA:Stylohyoid muscle.l','ZA:Platysma.l:Tendon'];meshes={pid:reader.read_part(atlas,parts[pid],frame=plan['frame'],catalogue_path=source,catalogue_sha=expected)for pid in ids}
    route=next(r for r in plan['routes']if r['name']=='glandular-1');recipient=np.asarray(route['recipient']['point_m']);stations=np.arange(24,37);rows=[]
    for station in stations:
        row=dict(station=int(station),center_m=centers[station].tolist(),recipient_distance_m=float(np.linalg.norm(centers[station]-recipient)),source_parent_geometry_sha256=parent.geometry_sha,source_ring_indices=rings[station].tolist(),surface_relations=[])
        for pid,m in meshes.items():
            d2,face,near=igl.point_mesh_squared_distance(centers[station:station+1],m.vertices,m.faces)
            row['surface_relations'].append(dict(part_id=pid,source_geometry_sha256=m.geometry_sha,distance_m=float(np.sqrt(d2[0])),face=int(face[0]),surface_point_m=near[0].tolist()))
        gland=row['surface_relations'][0]['distance_m'];bone=row['surface_relations'][1]['distance_m'];stylo=row['surface_relations'][2]['distance_m'];tendon=row['surface_relations'][3]['distance_m']
        row['allowed_cervical_departure']=bool(26<=station<=34 and station!=32 and gland<=.005 and min(bone,stylo,tendon)>.00055 and np.linalg.norm(centers[station]-centers[32])>.0012)
        rows.append(row)
    allowed=[r for r in rows if r['allowed_cervical_departure']];assert allowed,'No actual gland-adjacent source cervical takeoff meets fixed context clearance'
    selected=min(allowed,key=lambda r:r['recipient_distance_m']+2*r['surface_relations'][0]['distance_m']);station=selected['station']
    old_delta=parent.vertices[rings[48]]-centers[48];_,_,vt=np.linalg.svd(old_delta,full_matrices=False);normal=vt[-1]
    if normal@(centers[49]-centers[47])<0:normal=-normal
    ref=old_delta[0]-normal*(old_delta[0]@normal);ref/=np.linalg.norm(ref);basis=np.c_[ref,np.cross(normal,ref),normal];template=old_delta@basis;template[:,2]=0;template/=np.linalg.norm(template,axis=1).mean()
    contract=deepcopy(json.loads((base/'male-left-facial-cervical-routes-targeted-c2-v13/construction-contract.json').read_bytes()));contract.update(schema='cervical-glandular-takeoff-contract/v16',
        allowed_native_parent_station_interval=[26,34],allowed_actual_stations=[r['station']for r in allowed],excluded_passed_glandular0_root_station=32,
        anatomical_basis='Native facial course beside the actual submandibular gland before the lower anterior mandibular turn; explicit gland/mandible/stylohyoid/tendon face witnesses define the permitted interval.',
        selected_takeoff_station=station,selected_takeoff_native_coordinates_m=centers[station].tolist(),source_profile_template_original_station=48,source_profile_template=template.tolist(),
        root_location_is_construction_choice=True,native_parent_surface_or_controls_moved=False,recipient_endpoint_or_calibre_changed=False,literal_passed_routes=['ascending-palatine','tonsillar','glandular-0'],
        new_geometry_constructed=False,installed=False)
    atomic_json(out/'construction-contract.json',contract);atomic_json(ROOT/'reports/ta98-male/facial-cervical-takeoff-contract-v16.json',contract)
    result=dict(schema='native-cervical-facial-takeoff/v16',source_catalogue_sha256=expected,allowed_stations=contract['allowed_actual_stations'],selected=selected,station_witnesses=rows,
                old_root_station=48,old_root_location_m=centers[48].tolist(),source_profile_template=template.tolist(),source_parent_vertices_unchanged=True,recipient_unchanged=True,installed=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(source)==expected;atomic_json(out/'takeoff-evidence.json',result);print(json.dumps(dict(allowed=result['allowed_stations'],selected_station=station,selected_native_point=centers[station].tolist(),surface_distances_mm={r['part_id']:r['distance_m']*1000 for r in selected['surface_relations']},elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
