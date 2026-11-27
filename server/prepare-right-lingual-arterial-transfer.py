"""Mirror the accepted family course, then bind actual right native sources.

The centreline is independently constructed; no private reference model is
copied. Profile sections come from the native right ECA. Recipient and nearby
tissue checks must pass on the right rather than inheriting left acceptance.
"""
from pathlib import Path
import json,sys,time,copy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
T=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';L=B/'male-left-lingual-routes-v5';P=B/'male-lingual-source-plan-v1';D=B/'male-right-lingual-routes-v1';D.mkdir(exist_ok=False);inputs=json.loads((P/'inputs.json').read_text());left=json.loads((L/'routes.json').read_text());records={r['role']:r for r in inputs['fixed_sources']};current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';sha=digest_file(current);assert sha=='978f5a0a25d505f1fcfcca9d2b4cd9fd07023ae22e1363f15131dcd7d0e750b1';a=json.loads(current.read_bytes());parts={p['id']:p for p in a['parts']};reader=AtlasReader()
for role in ['external-carotid-right','greater-horn-right','hyoglossus-right','genioglossus-right','sublingual-gland-right','inferior-longitudinal-right','tongue-mucosa']:
    r=records[role];assert parts[r['part_id']]==r['source_record'];m=reader.read_part(a,parts[r['part_id']],frame=inputs['frame'],catalogue_path=current,catalogue_sha=sha);assert m.geometry_sha==r['geometry_sha256']and m.source['bank_sha256']==r['bank_sha256']
mirror=np.array([-1.,1.,1.]);right=copy.deepcopy(left);right.update(schema='right-lingual-source-bound-route-transfer/v1',source_catalogue_sha256=sha,source_inputs_sha256=digest_file(P/'inputs.json'),fixed_native_carotid=records['external-carotid-right'],left_course_contract_sha256=digest_file(L/'routes.json'),left_install_verification=str(B/'male-left-lingual-install-proposal-v1/independent-after-install.json'),anatomical_acceptance=False,installed=False)
root=np.array(records['external-carotid-right']['source_station_components'][0]['centres_m'][left['source_origin_station']]);oldroot=np.array(left['fixed_native_carotid']['source_station_components'][0]['centres_m'][left['source_origin_station']]);delta=root-oldroot*mirror
for route in right['routes']:
    z=np.load(route['guide']);points=z['points']*mirror
    if route['name']=='lingual':
        distance=np.r_[0,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))];s=np.clip(distance/.008,0,1);weight=1-s*s*(3-2*s);points+=weight[:,None]*delta;points[0]=root
    if route['recipient']is not None:route['recipient']['point_m']=(np.array(route['recipient']['point_m'])*mirror).tolist();route['recipient']['source_geometry_sha256']=Mesh.load(ROOT/records['sublingual-gland-right'if route['name']=='sublingual'else'hyoglossus-right']['mesh']).geometry_sha if route['name']in ['sublingual','suprahyoid']else route['recipient']['source_geometry_sha256']
    if route['parent']=='external-carotid-left':route['parent']='external-carotid-right'
    out=D/(route['name']+'-guide.npz');atomic_npz(out,points=points,outer_radii=z['outer_radii'],inner_radii=z['inner_radii'],frame=np.array(inputs['frame']));route.update(guide=str(out),guide_sha256=digest_file(out))
for name in ['whole-tongue-fixed-reference.npz','inner-tongue-fixed-reference.npz']:Mesh.load(L/name).save(D/name)
right['interfaces']=[dict(parent=r['parent'],child=r['name'],blood_opening_required=True,verified=False)for r in right['routes']];right['actual_right_parent_origin_adjustment_m']=delta.tolist();right['elapsed_seconds']=time.monotonic()-T;atomic_json(D/'routes.json',right);assert digest_file(current)==sha;print(json.dumps(dict(output=str(D),source_catalogue_sha256=sha,right_parent_origin_delta_m=delta.tolist(),anatomical_acceptance=False,elapsed_seconds=right['elapsed_seconds'])))
