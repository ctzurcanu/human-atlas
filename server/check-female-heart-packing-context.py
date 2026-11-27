"""Compare one unchanged thoracic context against the old and proposed heart."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--context',choices=['lungs','diaphragm','oesophagus','bones'],required=True);args=ap.parse_args()
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);start=time.monotonic()
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(path.read_text());cpath=h/args.candidate/'atlas.json';candidate=json.loads(cpath.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-female-ta98:assembled-metres'
def read(atlas,p):return reader.read_part(atlas,p,frame=frame,catalogue_path=path)
heart_id='TA98HEART:female:epicardium';old=read(a,parts[heart_id]);new=read(candidate,next(p for p in candidate['parts']if p['id']==heart_id))
if args.context=='lungs':ids=a['pulmonarySegmentIntegration']['partIds']
elif args.context=='diaphragm':ids=a['finiteDiaphragmIntegration']['partIds']
elif args.context=='oesophagus':ids=a['oesophagealWallGeometry']['partIds']
else:
 bounds=np.array([np.minimum(old.bounds[0],new.bounds[0]),np.maximum(old.bounds[1],new.bounds[1])])
 ids=[p['id']for p in a['parts']if p['id'].startswith('LOCAL:female:')and p['system']=='skeletal'and not p.get('suppressed')and np.all(np.array(p['bounds'])[1]>=bounds[0])and np.all(np.array(p['bounds'])[0]<=bounds[1])]
assert ids
meshes=[read(a,parts[pid])for pid in ids]
for pid in ids:assert parts[pid]==next(p for p in candidate['parts']if p['id']==pid)
offsets=np.r_[0,np.cumsum([len(m.vertices)for m in meshes])]
v,iv=np.unique(np.concatenate([m.vertices for m in meshes]),axis=0,return_inverse=True)
f=np.concatenate([iv[m.faces+offsets[i]]for i,m in enumerate(meshes)])
_,first,count=np.unique(np.sort(f,axis=1),axis=0,return_index=True,return_counts=True)
assert count.max()<=2
f=f[first[count==1]];fixed=Mesh(v,f,frame);fixed.save(out/'fixed-context.npz')
rows=[]
for name,heart in [('before',old),('candidate',new)]:
 _,proof=pair_contacts(heart,fixed,0.,max_pairs=250000,max_tests=2000000)
 assert proof['complete_within_pair_budget']
 rows.append(dict(name=name,**proof));atomic_json(out/(name+'.json'),proof)
 print(name,args.context,'crossing pairs',proof['crossing_face_pairs'],flush=True)
atomic_json(out/'context.json',dict(context=args.context,source_catalogue_sha256=digest_file(path),candidate_catalogue_sha256=digest_file(cpath),fixed_part_ids=ids,fixed_geometry_sha256=fixed.geometry_sha,contacts=rows,
 candidate_no_crossings=not rows[-1]['crossing_face_pairs'],anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
