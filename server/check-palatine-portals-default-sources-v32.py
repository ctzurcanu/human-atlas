"""Coarse portal source measurements against runtime-retained tissue only."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader


def main():
    D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';P=D/'global-corridor-analysis-v31/evidence.json';e=json.loads(P.read_bytes());binding=json.loads((D/'source-binding-v19/evidence.json').read_bytes());source=ROOT/binding['snapshot'];expected=binding['snapshot_sha256'];assert digest_file(source)==expected;atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};runtime=json.loads((D/'nearby-unresolved-participation-frame-v25.json').read_bytes());reader=AtlasReader();points=np.array([r['point_m']for r in e['coarse_portal_landmarks']]);frame='model:local-ta98:assembled-metres';safe=np.ones(len(points),bool);queries=[]
    for pid in ['ZA:Gingiva','ZA:Soft palate','O3M:Trapezius muscle.l']:
        row=next(r for r in runtime['rows']if r['id']==pid);mesh=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=expected);assert mesh.geometry_sha==row['registered_geometry_binding']['geometry_sha256'];d2,faces,near=igl.point_mesh_squared_distance(points,mesh.vertices,mesh.faces);distance=np.sqrt(d2);inside=abs(igl.winding_number(mesh.vertices,mesh.faces,points))>.99;safe&=~inside&(distance>=.0006201);queries.append(dict(part_id=pid,geometry_sha256=mesh.geometry_sha,minimum_distance_m=float(distance.min()),inside_count=int(inside.sum()),passing_points=int((~inside&(distance>=.0006201)).sum())))
    candidates=[]
    for i,r in enumerate(e['coarse_portal_landmarks']):
        if not safe[i]:continue
        candidates.append(r)
    chosen=[]
    for y in [1.635,1.645,1.650,1.655,1.660]:
        if candidates:chosen.append(min(candidates,key=lambda x:abs(x['point_m'][1]-y)))
    result=dict(schema='native-global-palatine-portal-default-check/v32',source_catalogue_sha256=expected,portal_evidence_sha256=digest_file(P),previous_physical_fixed_objects=93,additional_default_runtime_tissues_tested=queries,portal_count=len(points),passing_portals=len(candidates),selected_coarse_anatomical_witnesses=chosen,source_origin_interval_measurements=e['native_origin_station_measurements'],not_a_guide_or_continuous_corridor_proof=True,all_open_physical_volume_unknowns_retained=True,hair_scope='The per-instance hair audit certifies held-v17 swept bounds only. A new global proposal must re-evaluate actual shaft bounds against its own domain or emitted geometry.',proposed_method='Retain native station8 initially. Build a new source-bound global clearance domain spanning proximal cervical parent, actual SG/SP opposing-surface portal set, source superior-constrictor upper-margin witnesses and unchanged PP recipient. Include all93 physical sources plus retained current-default tissues; exclude alternatives only by executed loader proof. Search free space with endpoint source-parent core and intended recipient-tissue masks, then fit proper C2 controls with all intermediates free. Do not use rejected v17 as a trust-region center. If station8 cannot enter, report that graph witness before considering another proximal native station0–14 with original station8 profile transported by proper rotations.',no_new_takeoff_selected=True,no_curve_or_mesh_constructed=True,no_pack_or_install=True)
    out=D/'global-corridor-analysis-v31/default-context-portal-check-v32.json';atomic_json(out,result);print(json.dumps(dict(queries=queries,passing_portals=len(candidates),selected=[dict(point=r['point_m'],surface_gap_mm=r['local_opposing_surface_separation_m']*1000,ligament_distances_mm={x['part_id']:x['distance_m']*1000 for x in r['source_witnesses']if'ligament'in x['part_id']})for r in chosen])))


if __name__=='__main__':main()
