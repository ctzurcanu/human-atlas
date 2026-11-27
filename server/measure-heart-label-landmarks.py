"""Reference statistics for cardiac labels; never export donor triangles/UVs."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from ta98_glb_geometry import load_roots
from geometry_tools.model import AtlasReader
from geometry_tools.fit import landmark_transform
from geometry_tools.common import atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();started=time.monotonic()
out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
directory=ROOT.parent/'3d/ext/meshes';metadata_path=directory/'964db2dd4f98052f03baa9ca5f2dbcae/assets.json';metadata=json.loads(metadata_path.read_text());names={r['cid']:r['name']for r in json.loads((directory/'en.json').read_text())['items']}
source=ROOT/'.local-models/models/Male_2024.glb'
ids={24297,19627,18198,20944,20999,18980,20592,18331,18598,22809,23777,20644,23809,20454,19988,18054,22621,21663,24354,23271,22800,18687,16459,18229,17776,22137,20044,18445,23412,23932,23251,20879,23056,17580,20690,19455,22898,18032,22227}
rows={r['sni'][0]['id']:r for r in metadata['structures']if r.get('cid')in ids};meshes=load_roots(source,set(rows),allow_zero_weight_targets=True)
def stats(v,f):
 t=v[f];area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1);q=area@t.mean(1)/area.sum();return q,float(area.sum())
by={}
for mesh in meshes:
 cid=rows[mesh['node']]['cid'];q,weight=stats(mesh['unique'],mesh['welded']);by.setdefault(cid,[]).append((mesh,q,weight))
source_centroid={cid:sum(q*w for m,q,w in ms)/sum(w for m,q,w in ms)for cid,ms in by.items()}
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};groups={g['key']:g for g in a['finiteCardiacIntegration']['termGroups']};reader=AtlasReader()
links={24297:'mitral-leaflet-18',19627:'mitral-leaflet-19',18198:'tricuspid-leaflet-14',20944:'tricuspid-leaflet-15',20999:'tricuspid-leaflet-16',18980:'papillary-28',20592:'papillary-30',18331:'papillary-31'}
x=[];y=[]
for cid,key in links.items():
 pid=groups[key]['partIds'][0];m=reader.read_part(a,parts[pid],frame='assembled-metres',catalogue_path=path);x.append(source_centroid[cid]);y.append(stats(m.vertices,m.faces)[0])
matrix,fit=landmark_transform(np.array(x),np.array(y),similarity=True,scale_bounds=(.001,.1));rows_out=[]
for cid,ms in sorted(by.items()):
 points=np.concatenate([m['unique']for m,q,w in ms])@matrix[:3,:3].T+matrix[:3,3];q=source_centroid[cid]@matrix[:3,:3].T+matrix[:3,3];_,_,axes=np.linalg.svd(points-q,full_matrices=False);projected=(points-q)@axes.T
 # A few aggregate extrema and a centroid describe an anatomical reference
 # location. Neither source vertex sequences nor donor faces are emitted.
 extent=np.quantile(projected,[.05,.95],axis=0);support={name:points[:,axis].max().item()for axis,name in enumerate(['maxX','maxY','maxZ'])};support.update({name:points[:,axis].min().item()for axis,name in enumerate(['minX','minY','minZ'])})
 rows_out.append(dict(cid=cid,name=names[cid],centroid_m=q.tolist(),principal_axes=axes.tolist(),principal_extent_90_m=extent.tolist(),bounds_m=[points.min(0).tolist(),points.max(0).tolist()],support=support,reference_only=True,donor_mesh_exported=False))
report=dict(schema='cardiac-label-landmark-statistics/v1',source_glb=str(source),source_sha256=digest_file(source),metadata_sha256=digest_file(metadata_path),target_catalogue_sha256=digest_file(path),matrix=matrix.tolist(),fit=fit,associations=links,landmarks=rows_out,donor_faces_exported=0,donor_texture_pixels_exported=0,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started);atomic_json(out/'report.json',report);print('Cardiac reference statistics',len(rows_out),'registration',fit,'elapsed',report['elapsed_seconds'],flush=True)
