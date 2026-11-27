"""Final bounded rejection/localization of an incomplete native field trial."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion'
    route_dir=base/'male-left-facial-cervical-routes-native-domain-v7';candidate=base/'male-left-facial-cervical-finite-network-native-domain-v7'
    plan=json.loads((route_dir/'routes.json').read_bytes());expected=plan['source_catalogue_sha256'];source=ROOT/plan['source_catalogue'];assert digest_file(source)==expected
    report_file=candidate/'stage3-rejection-evidence.json';assert not report_file.exists()
    fields_file=candidate/'ascending-palatine-fields.npz';field=np.load(fields_file);p=field['guide_points'];faces=field['faces'];outer=Mesh(field['fields'][0],faces,str(field['frame']))
    previous=json.loads((candidate/'ascending-palatine-paired-field.json').read_bytes());bad=np.array(previous['invalid_faces_by_interval'][0],int)
    count=field['fields'].shape[1]//len(p);stations=np.unique(faces[bad].min(1)//count)
    a,b=np.diff(p,axis=0)[:-1],np.diff(p,axis=0)[1:];angle=np.arccos(np.clip(np.einsum('ij,ij->i',a,b)/(np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)),-1,1))
    cross=np.linalg.norm(np.cross(a,b),axis=1);radius=np.full(len(a),np.inf);np.divide(np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)*np.linalg.norm(a+b,axis=1),2*cross,out=radius,where=cross>1e-25)
    segments=[]
    for group in np.split(stations,np.flatnonzero(np.diff(stations)>1)+1):
        if not len(group):continue
        lo,hi=int(group.min()),int(group.max());local=np.arange(max(0,lo-1),min(len(angle),hi+1))
        worst=int(local[np.argmin(radius[local])]);segments.append(dict(first_interval=lo,last_interval=hi,source_points_bounds_m=[p[lo:hi+2].min(0).tolist(),p[lo:hi+2].max(0).tolist()],
            worst_three_point_station=worst+1,worst_turn_angle_degrees=float(np.degrees(angle[worst])),minimum_three_point_curvature_radius_m=float(radius[worst]),
            actual_three_center_points_m=p[worst:worst+3].tolist()))
    atlas=json.loads(source.read_bytes());parts={r['id']:r for r in atlas['parts']};reader=AtlasReader()
    stage1=json.loads((route_dir/'stage1-evidence.json').read_bytes());blockers=stage1['routes'][0]['fixed_boundaries'];contacts=[];incomplete=False
    for blocker in blockers:
        if time.monotonic()-began>45:incomplete=True;break
        pid=blocker['part_id'];target=reader.read_part(atlas,parts[pid],frame=outer.frame,catalogue_path=source,catalogue_sha=expected)
        if not(np.all(target.bounds[1]>=outer.bounds[0])and np.all(target.bounds[0]<=outer.bounds[1])):continue
        _,summary=pair_contacts(outer,target,0.,max_pairs=100000,max_tests=1000000)
        contacts.append(dict(part_id=pid,source_geometry_sha256=target.geometry_sha,summary=summary))
    result=dict(schema='rejected-native-facial-field-localization/v7',source_catalogue_sha256=expected,stage=3,maximum_stages=3,
        rejected_at_stage=2,failure='676 triangular finite wall cells invert at source-calibre tight turns; no Boolean union or emitted owner bank was created.',
        preserved_candidate_field=str(fields_file.relative_to(ROOT)),preserved_candidate_file_sha256=digest_file(fields_file),actual_outer_sheet_geometry_sha256=outer.geometry_sha,
        invalid_wall_faces=len(bad),invalid_longitudinal_intervals=len(stations),bad_interval_groups=segments,
        actual_unjoined_outer_source_context=contacts,context_incomplete_due_internal_budget=incomplete,
        context_scope='Actual unjoined outer field triangles only. No claim of full emitted-bank or whole occupied-volume acceptance: finite wall cell construction failed first.',
        remaining_routes_not_constructed=['tonsillar','glandular-0','glandular-1'],source_calibres_not_reduced=True,skeleton_and_source_unchanged=True,
        further_series_stopped=True,blender_queued=False,installed=False,term_count_added=0,elapsed_seconds=time.monotonic()-began)
    assert digest_file(source)==expected;atomic_json(report_file,result)
    print(json.dumps(dict(invalid_wall_faces=len(bad),invalid_longitudinal_intervals=len(stations),bad_groups=len(segments),minimum_curvature_radius_mm=min(r['minimum_three_point_curvature_radius_m']for r in segments)*1000,
                         actual_context_crossings=[dict(part_id=r['part_id'],pairs=r['summary']['crossing_face_pairs'])for r in contacts if r['summary']['crossing_face_pairs']],
                         context_incomplete=incomplete,elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
