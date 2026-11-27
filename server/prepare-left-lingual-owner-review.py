"""Expose the emitted owners and one local longitudinal communicating cut."""
from pathlib import Path
import json,sys,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
args=argparse.ArgumentParser();args.add_argument('--side',choices=['left','right'],default='left');args.add_argument('--version',default='v5');args.add_argument('--route-version');args.add_argument('--layout',default='owned-bank-v2');args=args.parse_args()
D=ROOT/('.local-models/geometry-workbench/circulatory-completion/male-'+args.side+'-lingual-finite-network-'+args.version);O=D/'owner-review-v1';O.mkdir(exist_ok=False);layout=json.loads((D/args.layout/'tissue.json').read_text());assert layout['passed'];plan=json.loads((D.parent/('male-'+args.side+'-lingual-routes-'+(args.route_version or args.version))/'routes.json').read_text());whole=[];colors=[[.75,.12,.10,1],[.94,.30,.13,1],[.93,.55,.10,1],[.94,.17,.34,1],[.78,.13,.48,1],[.69,.19,.51,1],[.96,.43,.42,1]]
for i,(name,pid)in enumerate(layout['part_ids'].items()):
    path=D/args.layout/(pid.replace(':','_')+'.npz');whole.append(dict(mesh=str(path),name=name,part_id=pid,tissue_draw_groups=layout['parts'][name]['tissueDrawGroups'],color=colors[i]))
atomic_json(O/'whole-spec.json',dict(frame=layout['frame'],resolution=1100,wall_seconds=60,assembled_tissues=True,meshes=whole,views=[dict(name='Actual seven owner cells, five left arterial terms',direction=[0,.2,1]),dict(name='Medial native-source recipients and new branch departure',direction=[1,.2,.2])]))
station=plan['source_origin_station'];parent=np.array(plan['fixed_native_carotid']['source_station_components'][0]['centres_m']);root=parent[station];axis=parent[station+1]-parent[station-1];guide=np.load(D/'lingual-fields.npz')['guide_points'];k=int(np.argmin(abs(np.r_[0,np.cumsum(np.linalg.norm(np.diff(guide,axis=0),axis=1))]-.005)));departure=guide[k]-root;normal=np.cross(axis,departure);normal/=np.linalg.norm(normal);local=[];ids={layout['part_ids'][name]for name in ['external-carotid-'+args.side,'lingual']}
for name in ['external-carotid-'+args.side,'lingual']:
    pid=layout['part_ids'][name];source=Mesh.load(D/args.layout/(pid.replace(':','_')+'.npz'));keep=[]
    for group in layout['parts'][name]['tissueDrawGroups']:
        if group.get('neighbor')not in ids:keep.extend(range(group['start']//3,(group['start']+group['count'])//3))
    faces=source.faces[np.array(keep,int)];mask=np.linalg.norm(source.vertices[faces].mean(1)-root,axis=1)<.014;faces=faces[mask];m=Mesh(source.vertices,faces,source.frame,source.attributes);path=O/(name+'-diagnostic-crop.npz');m.save(path);local.append(dict(mesh=str(path),name=name+' actual emitted wall; diagnostic 14 mm crop',color=colors[list(layout['part_ids']).index(name)]))
atomic_json(O/'focused-takeoff-spec.json',dict(frame=layout['frame'],resolution=1100,wall_seconds=60,meshes=local,views=[dict(name='Actual finite ECA and new lingual takeoff, local surface',direction=(normal+.2*axis/np.linalg.norm(axis)).tolist()),dict(name='Longitudinal communicating blood lumen through actual ECA branch mouth',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=root.tolist(),clear_outer=True))]))
atomic_json(O/'evidence.json',dict(actual_bank_sha256=layout['bank_sha256'],bank_geometry_unchanged=True,diagnostic_crop_radius_m=.014,root_source_station=station,section_contains_actual_parent_axis_and_lingual_departure=True,section_normal=normal.tolist(),lumen_cap_added=False,anatomical_acceptance=False));print(json.dumps(dict(whole_spec=str(O/'whole-spec.json'),focused_spec=str(O/'focused-takeoff-spec.json'))))
