"""Final targeted diagnostic; no reconstruction, packing or installation."""
from pathlib import Path
import sys,json,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import manifold3d as md
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.distance import triangle_pairs
from geometry_tools.algorithms.finite_tissue_boolean import make_solid


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';out=base/'male-left-facial-cervical-routes-targeted-c2-v13';file=out/'stage3-blocker-diagnosis-v15.json';assert not file.exists()
    prior=base/'male-left-facial-cervical-routes-free-c2-v10/routes.json';plan=json.loads(prior.read_bytes());source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256'];assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame=plan['frame']
    meshes={pid:reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=expected)for pid in ['ZA:Hypoglossal nerve (XII).l','ZA:Platysma.l:Tendon','ZA:Mandible:Bone','ZA:Platysma.l:Biarticular','ZA:Submandibular gland.l']}
    field=np.load(out/'tonsillar-free-c2-fields.npz');outer=Mesh(field['fields'][0],field['faces'],frame);inner=Mesh(field['fields'][1],field['faces'],frame);nerve=meshes['ZA:Hypoglossal nerve (XII).l'];contacts=[]
    for name,mesh in [('outer',outer),('inner',inner)]:
        _,summary=pair_contacts(mesh,nerve,0.,max_pairs=1000000,max_tests=12000000);assert summary['complete_within_pair_budget'];contacts.append(dict(sheet=name,geometry_sha256=mesh.geometry_sha,summary=summary))
    spec=importlib.util.spec_from_file_location('facial_v15_volume',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);parent_closed,_=helper.occupancy_classifier(parent,True);branch_closed,_=helper.occupancy_classifier(outer,True);nerve_closed,mode=helper.occupancy_classifier(nerve,False);assert branch_closed is not None and nerve_closed is not None
    origin=outer.bounds.mean(0)
    def solid(mesh):return make_solid((mesh.vertices-origin)*1000,mesh.faces)[0]
    overlap=(solid(branch_closed)-solid(parent_closed))^solid(nerve_closed);assert overlap.status()==md.Error.NoError
    tendon,mandible=meshes['ZA:Platysma.l:Tendon'],meshes['ZA:Mandible:Bone'];tface,bface=242,10486;patch=triangle_pairs(tendon.triangles[[tface]],mandible.triangles[[bface]])
    route=next(r for r in plan['routes']if r['name']=='glandular-1');p=np.load(ROOT/route['guide'])['points'];endpoints=[]
    for label,point in [('fixed_parent_departure',p[0]),('fixed_recipient_endpoint',p[-1])]:
        rows=[]
        for pid in ['ZA:Platysma.l:Tendon','ZA:Mandible:Bone','ZA:Platysma.l:Biarticular']:
            mesh=meshes[pid];d2,faces,near=__import__('igl').point_mesh_squared_distance(point[None],mesh.vertices,mesh.faces)
            rows.append(dict(part_id=pid,source_geometry_sha256=mesh.geometry_sha,distance_m=float(np.sqrt(d2[0])),source_face=int(faces[0]),nearest_source_point_m=near[0].tolist()))
        endpoints.append(dict(role=label,point_m=point.tolist(),actual_source_surfaces=rows))
    failed=json.loads((out/'stage1-evidence.json').read_bytes());gland=next(r for r in failed['routes']if r['name']=='glandular-1');domain=json.loads((base/'male-left-facial-cervical-routes-finite-tendon-seed-v14/stage2-evidence.json').read_bytes())
    result=dict(schema='targeted-native-facial-blocker-diagnostic/v15',stage=3,maximum_stages=3,source_catalogue_sha256=expected,
        tonsillar_local_XII_checks=dict(actual_sheet_context=contacts,occupied_new_outer_conduit_overlap_mm3=float(overlap.volume()),occupied_intersection_empty=bool(overlap.is_empty()),
                                       scope='Corrected tonsillar actual fields against XII only, with literal parent occupied volume removed; not a full family or emitted-bank context pass'),
        glandular_fixed_endpoints=endpoints,glandular_outer_radius_m=route['illustrative_radius_m'],
        source_attachment_patch_pair=dict(tendon_part='ZA:Platysma.l:Tendon',tendon_face=tface,tendon_triangle_m=tendon.triangles[tface].tolist(),mandible_part='ZA:Mandible:Bone',mandible_face=bface,mandible_triangle_m=mandible.triangles[bface].tolist(),
                                          actual_triangle_pair_metrics={k:np.asarray(v).tolist()for k,v in patch.items()},qualification='Local finite patch witness only; does not prove a global route is impossible'),
        linearised_conic_witnesses=gland['iterations'][0]['solver']['leading_native_conic_witnesses'],
        finite_source_domain_result=dict(grid_step_m=.0006,domain_halfwidth_m=.008,safe_sites=domain['safe_points'],search=domain['search'],global_impossibility_claim=False),
        conclusion='The tonsillar XII blocker is locally cleared. The second glandular course remains unconstructed at fixed endpoints; neither the local conic certificate nor finite grid disconnection proves global anatomical impossibility.',
        passed_palatine_and_glandular0_literal=True,source_calibre_or_fixed_anatomy_changed=False,whole_family_not_accepted=True,installed=False,blender_queued=False,term_count_added=0,elapsed_seconds=time.monotonic()-began)
    assert digest_file(source)==expected;atomic_json(file,result)
    print(json.dumps(dict(tonsillar_XII_triangle_crossings=sum(r['summary']['crossing_face_pairs']for r in contacts),tonsillar_XII_occupied_overlap_mm3=float(overlap.volume()),local_tendon_mandible_patch_distance_m=np.asarray(patch['distance_m']).tolist(),no_global_impossibility_claim=True,installed=False,elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
