"""Analytic segment/surface audit of saved v33 diagnostic graph gaps only."""
from pathlib import Path
import json,sys,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.distance import triangle_pairs

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'global-palatine-graph-v33';c=json.loads((O/'construction-contract.json').read_bytes());e=json.loads((O/'domain-evidence.json').read_bytes());post=json.loads((O/'postmortem-v34.json').read_bytes());source=ROOT/c['source_catalogue'];assert digest_file(source)==c['source_catalogue_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';pdir=ROOT/'.local-models/geometry-workbench/circulatory-completion';plan=json.loads((pdir/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());spec=importlib.util.spec_from_file_location('palatine_v35_occupancy',ROOT/'server/check-facial-cervical-feasibility.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);fixed={}
    for row in e['source_queries']:
        pid=row['part_id']
        if pid in parts:m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=c['source_catalogue_sha256']);allow=parts[pid].get('system')in{'arterial','venous'}
        elif pid=='CANDIDATE:glandular-0:passed-fixed-envelope':
            r=plan['extra_fixed_meshes'][0];assert digest_file(ROOT/r['mesh'])==r['mesh_sha256'];m=Mesh.load(ROOT/r['mesh']);allow=True
        else:
            name=pid.split(':')[1];a=np.load(D/(name+'-fields.npz'));m=Mesh(a['fields'][0],a['faces'],frame);allow=True
        assert m.geometry_sha==row['geometry_sha256'];cl,mode=h.occupancy_classifier(m,allow);fixed[pid]=(m,cl,mode)
    z=np.load(O/'clearance-domain.npz');q=z['points'];free=z['clearance_lower_bound_m']>=c['source_profile_radius_m']+c['minimum_clearance_m']+float(z['step'])/2;ids=np.flatnonzero(free);end=np.array(c['recipient']['point_m']);ds,ix=cKDTree(q[ids]).query(end,k=8)
    segments=[dict(name='recipient-nearest-free-cell-'+str(j),first=end,second=q[ids[k]],exemption='ZA:Palatopharyngeus muscle.l')for j,(d,k)in enumerate(zip(ds,ix))if d<=.002]
    segments+=[dict(name='source-to-portal-component-'+str(r['component']),first=np.array(r['source_component_point_m']),second=np.array(r['portal_component_point_m']),exemption=None)for r in post['portal_component_witnesses']]
    radius=c['source_profile_radius_m'];threshold=radius+c['minimum_clearance_m'];rows=[]
    for s in segments:
        a,b=s['first'],s['second'];segment=np.array([a,b,a]);low=np.minimum(a,b)-threshold;high=np.maximum(a,b)+threshold;checks=[]
        for pid,(m,cl,mode)in fixed.items():
            if pid==s['exemption']:
                checks.append(dict(part_id=pid,intended_PP_terminal_interface_exempt=True,complete_segment_inside_two_mm_recipient_mask=bool(max(np.linalg.norm(a-end),np.linalg.norm(b-end))<=.002)));continue
            if np.any(m.bounds[1]<low)or np.any(m.bounds[0]>high):continue
            triangles=m.triangles;selected=np.flatnonzero(np.all(triangles.max(1)>=low,axis=1)&np.all(triangles.min(1)<=high,axis=1))
            if len(selected):
                r=triangle_pairs(np.broadcast_to(segment,(len(selected),3,3)),triangles[selected]);j=int(np.argmin(r['distance_m']));distance=float(r['distance_m'][j]);lower=float(distance-r['numerical_roundoff_m'][j]);witness=dict(face=int(selected[j]),source_point_m=r['second'][j].tolist(),segment_point_m=r['first'][j].tolist(),distance_m=distance,distance_lower_bound_m=lower)
            else:lower=threshold;witness=dict(no_source_triangle_AABB_within_expanded_segment=True)
            # If the complete segment stays farther than the radius from the
            # source, a verified outside starting point establishes disjoint
            # closed-source occupancy throughout that segment.
            inside_start=bool(abs(igl.winding_number(cl.vertices,cl.faces,a[None])[0])>.99)if cl is not None else None
            checks.append(dict(part_id=pid,geometry_sha256=m.geometry_sha,classifier=mode,inside_valid_closed_source_at_start=inside_start,surface_clear_at_actual_radius=lower>=threshold,occupied_clear_if_closed=None if cl is None else bool(lower>=threshold and not inside_start),**witness))
        failures=[r for r in checks if r.get('surface_clear_at_actual_radius')is False or r.get('inside_valid_closed_source_at_start')is True]
        rows.append(dict(name=s['name'],first_m=a.tolist(),second_m=b.tolist(),length_m=float(np.linalg.norm(b-a)),surface_and_valid_closed_source_gate_passed=not failures,failures=failures,near_source_checks=checks,open_source_occupancy_unresolved=True))
    result=dict(schema='global-palatine-graph-analytic-segment-audit/v35',source_catalogue_sha256=c['source_catalogue_sha256'],domain_sha256=e['domain_sha256'],postmortem_sha256=digest_file(O/'postmortem-v34.json'),radius_m=radius,margin_m=c['minimum_clearance_m'],
       source_meshes=99,method='Full finite straight segment represented as a degenerate triangle; existing batched triangle-pair minimization includes segment/face intersections, edge/edge, endpoint/face and vertex/segment distances. Conservative AABB exclusion uses the complete expanded segment. No nearest-vertex approximation or shrunken radius.',segments=rows,
       no_new_route_or_geometry=True,no_anatomical_acceptance=True,no_pack_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'analytic-segment-audit-v35.json',result);print(json.dumps(dict(segments=[dict(name=r['name'],pass_=r['surface_and_valid_closed_source_gate_passed'],failures=[dict(part_id=f['part_id'],distance_m=f.get('distance_m'),inside=f.get('inside_valid_closed_source_at_start'))for f in r['failures']])for r in rows],elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
