"""A finite five-size comparison after the saved translation-only rejection."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import manifold3d as md
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
ap=argparse.ArgumentParser();ap.add_argument('--joint',action='store_true');args=ap.parse_args()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/('final-context-v217'if args.joint else'final-context-v212');out.mkdir(exist_ok=False)
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(cat.read_text());reader=AtlasReader();parts={p['id']:p for p in a['parts']};frame='model:local-female-ta98:assembled-metres'
def solid(m):
 v,inv=np.unique(m.vertices,axis=0,return_inverse=True);s=md.Manifold(md.Mesh64(np.ascontiguousarray(v),np.ascontiguousarray(inv[m.faces],dtype=np.uint64)));assert s.status()==md.Error.NoError;return s
ids=list(a['finiteDiaphragmIntegration']['partIds'])+[p['id']for p in a['parts']if p['id'].startswith('GASTRIC:female:019:')]
fixed=None
for pid in ids:
 s=solid(reader.read_part(a,parts[pid],frame=frame,catalogue_path=cat));fixed=s if fixed is None else fixed+s
fit=a['liverTransferIntegration'];source=Mesh.load(b/'vessel-containment-v157/liver-envelope.npz');v=source.vertices*fit['scale']+fit['translation_m'];center=v.mean(0);liver=solid(Mesh(v-center,source.faces,frame))
previous=json.loads((b/'final-context-v210/fit.json').read_text());best=min(previous['trials'],key=lambda r:r['overlap_m3']);delta=np.asarray(best['translation_m']);trials=[];chosen=None;start=time.monotonic()
class Found(Exception):pass
def objective(x):
 global chosen
 factor=float(x[0]);delta=np.asarray(x[1:])/1000
 q=liver.scale([factor]*3).translate(center+delta)^fixed;assert q.status()==md.Error.NoError
 row=dict(additional_scale=factor,translation_m=delta.tolist(),center_m=center.tolist(),overlap_m3=float(q.volume()),empty=q.is_empty());trials.append(row)
 if row['empty']:chosen=row;raise Found()
 return max(0,row['overlap_m3'])*1e10
try:
 if args.joint:
  seed=np.r_[.985,delta*1000];simplex=np.r_[seed[None],seed[None]+np.diag([-.004,.3,.3,.3])]
  minimize(objective,seed,method='Nelder-Mead',bounds=[(.94,1),(-10,5),(-3,10),(-6,6)],options=dict(maxiter=160,maxfev=300,initial_simplex=simplex,xatol=1e-8,fatol=1e-16))
 else:
  for factor in [.9975,.995,.99,.985,.98]:objective(np.r_[factor,delta*1000])
except Found:pass
record=dict(passed=chosen is not None,chosen=chosen,trials=trials,source_sha256=digest_file(cat),fixed_parts=ids,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'fit.json',record)
if chosen is None:raise ValueError('Five bounded size candidates rejected; no installation')
