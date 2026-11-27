"""Final exact fit inequalities for the rejected analytic fillet attempt."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';out=base/'male-left-facial-cervical-routes-curvature-v8'
    previous=json.loads((out/'stage1-evidence.json').read_bytes());assert not previous['all_analytic_and_finite_field_checks_pass']
    file=out/'stage2-exact-constraints.json';assert not file.exists();row=previous['routes'][-1];p=np.array(row['controls_m']);radius=row['minimum_analytic_bend_radius_m']
    delta=np.diff(p,axis=0);length=np.linalg.norm(delta,axis=1);direction=delta/length[:,None];angle=np.arccos(np.clip(np.einsum('ij,ij->i',direction[:-1],direction[1:]),-1,1));coefficient=np.r_[0,np.tan(angle/2),0]
    required=radius*(coefficient[:-1]+coefficient[1:]);bad=np.flatnonzero(required>=length-1e-12)
    maximum=np.divide(length,coefficient[:-1]+coefficient[1:],out=np.full_like(length,np.inf),where=(coefficient[:-1]+coefficient[1:])>0)
    route_dir=base/'male-left-facial-cervical-routes-native-domain-v7';plan=json.loads((route_dir/'routes.json').read_bytes());source=ROOT/plan['source_catalogue'];expected=previous['source_catalogue_sha256'];assert digest_file(source)==expected
    catalogue=json.loads(source.read_bytes());parts={part['id']:part for part in catalogue['parts']};reader=AtlasReader();witnesses=p[np.unique(np.r_[bad,bad+1])]
    domain=json.loads((route_dir/'stage1-evidence.json').read_bytes());nearest=[]
    for blocker in domain['routes'][0]['fixed_boundaries']:
        pid=blocker['part_id'];mesh=reader.read_part(catalogue,parts[pid],frame=plan['frame'],catalogue_path=source,catalogue_sha=expected)
        d2,faces,points=igl.point_mesh_squared_distance(witnesses,mesh.vertices,mesh.faces)
        nearest.append(dict(part_id=pid,source_geometry_sha256=mesh.geometry_sha,minimum_bad_corner_surface_distance_m=float(np.sqrt(d2.min())),
                            nearest_bad_corner_m=witnesses[int(d2.argmin())].tolist(),source_face=int(faces[int(d2.argmin())]),source_surface_point_m=points[int(d2.argmin())].tolist()))
    constraints=[dict(segment=int(i),first_control_m=p[i].tolist(),last_control_m=p[i+1].tolist(),segment_length_m=float(length[i]),
                      required_adjacent_fillet_trim_length_m=float(required[i]),trim_overlap_m=float(required[i]-length[i]),
                      maximum_fillet_radius_on_this_fixed_segment_m=float(maximum[i]),required_profile_safe_fillet_radius_m=radius)for i in bad]
    result=dict(schema='analytic-facial-fillet-rejection/v8',stage=2,maximum_stages=2,source_catalogue_sha256=expected,
        source_outer_profile_maximum_radius_m=row['outer_profile_maximum_radius_m'],minimum_required_centerline_bend_radius_m=radius,
        failed_segment_constraints=constraints,maximum_radius_supported_by_all_fixed_controls_m=float(maximum.min()),actual_native_corner_surface_witnesses=sorted(nearest,key=lambda r:r['minimum_bad_corner_surface_distance_m']),
        single_strategy_result='No continuous nonoverlapping circular fillet can retain these controls and the required radius. Trial is stopped, not thinned or accepted.',
        previous_clearance_domain_qualification='The predecessor outer field cleared 17 fixed source triangle blockers but had inverted wall cells; no new valid wall or emitted bank exists.',
        no_radius_reduction=True,no_source_profile_or_anatomy_changed=True,blender_queued=False,installed=False,term_count_added=0,elapsed_seconds=time.monotonic()-began)
    assert digest_file(source)==expected;atomic_json(file,result)
    print(json.dumps(dict(failed_segments=len(constraints),required_radius_mm=radius*1000,maximum_radius_supported_mm=float(maximum.min()*1000),
                         maximum_trim_overlap_mm=max(r['trim_overlap_m']for r in constraints)*1000,nearest_fixed_sources=[dict(part_id=r['part_id'],distance_mm=r['minimum_bad_corner_surface_distance_m']*1000)for r in result['actual_native_corner_surface_witnesses'][:4]],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
