"""Finite entry-track checks against real tissue; static valve closure excluded."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import igl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection
ap=argparse.ArgumentParser();ap.add_argument('--source');ap.add_argument('--bank');ap.add_argument('--stage',required=True);ap.add_argument('--fit',default='outflow-protected-fit-v287');args=ap.parse_args()
assert bool(args.source)!=bool(args.bank);start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
exclude={'ZA:Left coronary leaflet','ZA:Non-coronary leaflet','ZA:Right coronary leaflet','ZA:Anterior semilunar leaflet of pulmonary valve','ZA:Left semilunar leaflet of pulmonary valve','ZA:Right semilunar leaflet of pulmonary valve'}
if args.bank:
 e=json.loads((h/args.bank/'emission.json').read_text());assert digest_file(h/args.bank/'cells.bin')==e['bank_sha256'];reader=AtlasReader();mini=dict(chunks=[dict(url='/local-models/'+str((h/args.bank/'cells.bin').relative_to(ROOT/'.local-models')))])
 meshes=[reader.read_part(mini,dict(row,chunk=0,geometryOrigin=[0.,0.,0.]),frame='heart-native-local-metres') for pid,row in e['parts'].items() if pid not in exclude]
 v=np.concatenate([m.vertices for m in meshes]);offset=np.r_[0,np.cumsum([len(m.vertices) for m in meshes])];f=np.concatenate([m.faces+offset[i] for i,m in enumerate(meshes)]);v,inv=np.unique(v,axis=0,return_inverse=True);f=inv[f]
 source=h/e['source_stage'];bound=e['bank_sha256']
else:
 source=h/args.source;z=np.load(source/'arrangement.npz');v=z['positions'];m=json.loads((source/'ownership.json').read_text());s={r['id']:z[r['faces']] for r in m['surfaces']};f=np.concatenate([s[u['surface']] if u['side']==1 else s[u['surface']][:,[0,2,1]] for pid,uses in m['cells'].items() if pid not in exclude for u in uses]);bound=digest_file(source/'arrangement.npz')
meta=json.loads((source/'ownership.json').read_text());_,first,count=np.unique(np.sort(f,axis=1),axis=0,return_index=True,return_counts=True);assert count.max()==2;f=f[first[count==1]]
mesh=Mesh(v,f,'heart-native-local-metres');tree=BVH(mesh);fields=json.loads((h/meta['outflow_fluid_stage']/'fields.json').read_text());placement=json.loads((h/args.fit/'placement.json').read_text());origin=np.array(meta['origin']);rows=[];tests=0
fit_metadata=h/args.fit/'ownership.json'
root_stages=json.loads(fit_metadata.read_text())['root_stages'] if fit_metadata.exists() else ['aorta-rounded-course-v276','pulmonary-spatial-course-v278']
for i,native in enumerate(root_stages):
 z=np.load(h/native/'source-walls.npz');centers=z['centers'] if 'centers' in z else np.array([z['inner_positions'][z['station_indices']==j].mean(0) for j in np.unique(z['station_indices'])]);curve=centers+np.array(placement['transforms'][i])[:3,3]-origin
 if curve[0,1]>curve[-1,1]:curve=curve[::-1]
 curve=np.r_[np.array(fields['records'][i]['target_point_local_m'])[None],curve];frames=rotation_minimizing_frames(curve)
 for radius in [0.,.003]:
  for angle in ([0.] if radius==0 else np.arange(8)*np.pi/4):
   track=curve+radius*(np.cos(angle)*frames[:,:,0]+np.sin(angle)*frames[:,:,1]);inside=abs(igl.winding_number(v,f,track));hits=[]
   for k,(a,b) in enumerate(zip(track,track[1:])):
    near=tree.query_bounds(np.minimum(a,b),np.maximum(a,b),distance_m=0.);tests+=len(near)
    if tests>250000 or time.monotonic()-start>25:raise TimeoutError('Outflow finite segment budget')
    for face in near:
     if exact_segment_triangle_intersection(a,b,v[f[face]])['intersects']:hits.append(dict(segment=k,face=int(face)))
   rows.append(dict(vessel=native,radius_m=radius,angle_rad=float(angle),stations_outside_tissue=bool(np.all(inside<.5)),tissue_crossings=hits,curve_m=track.tolist()))
result=dict(passed=all(r['stations_outside_tissue'] and not r['tissue_crossings'] for r in rows),geometry_binding=bound,
 excluded_static_semilunar_cusp_ids=sorted(exclude),tracks=rows,complete_segment_tests=True,triangle_tests=tests,
 elapsed_seconds=time.monotonic()-start,full_valve_motion_and_distal_circulation_not_established=True)
atomic_json(out/'paths.json',result);print(json.dumps(dict(passed=result['passed'],tracks=len(rows),crossings=sum(len(r['tissue_crossings']) for r in rows),elapsed_seconds=time.monotonic()-start)),flush=True)
raise SystemExit(0 if result['passed'] else 2)
