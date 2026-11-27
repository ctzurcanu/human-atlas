"""Read-only postmortem of the completed v33 grid; no route reconstruction."""
from pathlib import Path
import sys,json,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.ndimage import label,binary_dilation,generate_binary_structure
from scipy.spatial import cKDTree
from collections import Counter
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'global-palatine-graph-v33';c=json.loads((O/'construction-contract.json').read_bytes());e=json.loads((O/'domain-evidence.json').read_bytes());g=json.loads((O/'graph-evidence.json').read_bytes());assert g['completed']and not g['route_found'];assert digest_file(O/'clearance-domain.npz')==e['domain_sha256'];source=ROOT/c['source_catalogue'];assert digest_file(source)==c['source_catalogue_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();z=np.load(O/'clearance-domain.npz');q=z['points'];clear=z['clearance_lower_bound_m'];shape=tuple(z['shape']);threshold=c['source_profile_radius_m']+c['minimum_clearance_m']+float(z['step'])/2;free=clear>=threshold;component,count=label(free.reshape(shape));labels=component.ravel();ids=np.flatnonzero(free);start=np.array(c['start_m']);end=np.array(c['recipient']['point_m']);band=np.array(c['SG_SP_band_m']);di,ix=cKDTree(q[ids]).query(start);source_component=int(labels[ids[ix]]);portal_ids=np.flatnonzero(free&(cKDTree(band).query(q)[0]<=.002));portal_components=sorted(set(int(x)for x in labels[portal_ids]));rows=[]
    for n in portal_components:
        mask=component==n;boundary=binary_dilation(mask,structure=generate_binary_structure(3,1))&~mask;boundary_ids=np.flatnonzero(boundary.ravel());blocked=boundary_ids[~free[boundary_ids]];lim=Counter(int(x)for x in z['limiting_source'][blocked]);closest_d,closest_i=cKDTree(q[labels==source_component]).query(q[mask.ravel()]);k=int(np.argmin(closest_d));portal_region_point=q[np.flatnonzero(mask.ravel())[k]];root_region_point=q[np.flatnonzero(labels==source_component)[closest_i[k]]]
        rows.append(dict(component=n,nodes=int(mask.sum()),portal_nodes=int((labels[portal_ids]==n).sum()),closest_source_component_free_node_gap_m=float(closest_d[k]),source_component_point_m=root_region_point.tolist(),portal_component_point_m=portal_region_point.tolist(),blocked_boundary_limiting_sources=[dict(part_id=e['source_queries'][i]['part_id'],grid_boundary_points=v)for i,v in lim.most_common(8)],not_a_proof_of_global_anatomical_impossibility=True))
    spec=importlib.util.spec_from_file_location('postmortem_occupancy',ROOT/'server/check-facial-cervical-feasibility.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);endpoint=[]
    for r in e['source_queries']:
        pid=r['part_id']
        if pid.startswith('CANDIDATE:'):
            name=pid.split(':')[1]
            if name=='glandular-0':continue
            a=np.load(D/(name+'-fields.npz'));m=Mesh(a['fields'][0],a['faces'],str(z['frame']));allow=True
        elif pid in parts:
            bb=np.array(parts[pid]['bounds']);box=float(np.linalg.norm(np.maximum(bb[0]-end,0)+np.maximum(end-bb[1],0)))
            if box>.005:continue
            m=reader.read_part(atlas,parts[pid],frame=str(z['frame']),catalogue_path=source,catalogue_sha=c['source_catalogue_sha256']);allow=parts[pid].get('system')in{'arterial','venous'}
        else:continue
        assert m.geometry_sha==r['geometry_sha256'];cl,mode=h.occupancy_classifier(m,allow);d2,f,p=igl.point_mesh_squared_distance(end[None],m.vertices,m.faces);inside=bool(abs(igl.winding_number(cl.vertices,cl.faces,end[None])[0])>.99)if cl is not None else None;distance=float(np.sqrt(d2[0]));endpoint.append(dict(part_id=pid,geometry_sha256=m.geometry_sha,surface_distance_m=distance,inside_valid_classifier=inside,signed_distance_m=-distance if inside else distance,classifier=mode,face=int(f[0]),surface_point_m=p[0].tolist(),intended_PP_terminal_tissue=pid=='ZA:Palatopharyngeus muscle.l'))
    nearest_ids=np.flatnonzero(free);ds,ks=cKDTree(q[nearest_ids]).query(end,k=8);nearest=[dict(point_m=q[nearest_ids[k]].tolist(),distance_m=float(d),component=int(labels[nearest_ids[k]]))for d,k in zip(ds,ks)]
    result=dict(schema='global-palatine-graph-postmortem/v34',source_catalogue_sha256=c['source_catalogue_sha256'],domain_sha256=e['domain_sha256'],analysis_only_no_proposal_change=True,
      threshold_m=threshold,free_grid_components=count,source_component=source_component,source_component_nodes=int((labels==source_component).sum()),portal_component_witnesses=rows,
      recipient_endpoint_m=end.tolist(),recipient_endpoint_source_witnesses=sorted(endpoint,key=lambda r:r['signed_distance_m']),nearest_free_grid_nodes_to_recipient=nearest,
      known_scope_limits=['The grid is conservative at 1 mm spacing: no route is not a global impossibility proof.','Parent/recipient mask transitions were not accepted as finite junctions.','117 additional AABB-retrieved source IDs were not classified by the executed runtime loader.','Open/nonmanifold source occupied volumes and complete oral fluid space remain unknown.'],
      no_C2_or_tube_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'postmortem-v34.json',result);print(json.dumps(dict(source_component=source_component,portal_components=rows,recipient_closest_sources=result['recipient_endpoint_source_witnesses'][:7],nearest_free_grid_nodes_to_recipient=nearest[:2],elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
