"""Find a bounded translation of the complete liver, never the skeleton."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
import manifold3d as md
from scipy.optimize import differential_evolution, minimize
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh

b=ROOT/'.local-models/geometry-workbench/liver-completion'
ap=argparse.ArgumentParser();ap.add_argument('--refine-saved',action='store_true');args=ap.parse_args()
out=b/('final-context-v210' if args.refine_saved else 'final-context-v209');out.mkdir(exist_ok=False)
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
a=json.loads(cat.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-female-ta98:assembled-metres'
def solid(m):
 v,inv=np.unique(m.vertices,axis=0,return_inverse=True)
 s=md.Manifold(md.Mesh64(np.ascontiguousarray(v),np.ascontiguousarray(inv[m.faces],dtype=np.uint64)))
 if s.status()!=md.Error.NoError:raise ValueError(('Invalid closed fixed tissue',m.source,s.status()))
 return s
fit=a['liverTransferIntegration'];source=Mesh.load(b/'vessel-containment-v157/liver-envelope.npz')
v=source.vertices*fit['scale']+fit['translation_m'];liver=solid(Mesh(v,source.faces,frame))
ids=list(a['finiteDiaphragmIntegration']['partIds'])
ids += [p['id'] for p in a['parts'] if p['id'].startswith('GASTRIC:female:019:')]
obstacles=[]
for pid in ids:
 m=reader.read_part(a,parts[pid],frame=frame,catalogue_path=cat)
 obstacles.append((pid,solid(m)))
fixed=obstacles[0][1]
for _,s in obstacles[1:]:fixed=fixed+s
assert fixed.status()==md.Error.NoError
trials=[];accepted=None;start=time.monotonic()
class Found(Exception):pass
def objective(x):
 global accepted
 q=liver.translate(x)^fixed;assert q.status()==md.Error.NoError
 row=dict(translation_m=list(map(float,x)),overlap_m3=float(q.volume()),empty=q.is_empty());trials.append(row)
 if q.is_empty():accepted=row;raise Found()
 return max(0,row['overlap_m3'])*1e9 + np.sum((np.asarray(x)/.01)**2)*.001
try:
 if args.refine_saved:
  prior=json.loads((b/'final-context-v209/fit.json').read_text())
  best=min(prior['trials'],key=lambda r:r['overlap_m3']);x=np.asarray(best['translation_m'])*1000
  simplex=np.r_[x[None],x[None]+np.eye(3)*.3]
  minimize(lambda y:objective(np.asarray(y)/1000),x,method='Nelder-Mead',options=dict(maxiter=100,maxfev=200,xatol=.0001,fatol=1e-12,initial_simplex=simplex))
 else:
  for x in [[-.008,0,0],[-.012,0,.004],[-.016,0,.008],[-.02,-.002,.008]]:objective(x)
  differential_evolution(objective,[(-.025,.005),(-.008,.008),(-.012,.018)],seed=209,popsize=4,maxiter=15,polish=False,workers=1,tol=1e-5)
except Found:pass
result=dict(passed=accepted is not None,chosen=accepted,trials=trials,source_sha256=digest_file(cat),fixed_parts=ids,rotation_degrees=0,additional_scale=1,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'fit.json',result)
print(json.dumps({k:v for k,v in result.items() if k!='trials'}),flush=True)
if not result['passed']:raise ValueError('Bounded joint stomach/diaphragm translation did not pass; candidate not installed')
