"""Bounded native-volume fitting with one uniform transform for the closure."""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
from scipy.optimize import differential_evolution
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'female-transfer-v202';a=json.loads((out/'atlas.json').read_text());fit=json.loads((out/'fit.json').read_text());reader=AtlasReader();frame='model:local-female-ta98:assembled-metres'
def solid(mesh):
 v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);s=md.Manifold(md.Mesh64(np.ascontiguousarray(v),np.ascontiguousarray(inv[mesh.faces],dtype=np.uint64)));assert s.status()==md.Error.NoError;return s
fixed=[solid(reader.read_part(a,next(p for p in a['parts']if p['id']==pid),frame=frame))for pid in a['finiteDiaphragmIntegration']['partIds']]
diaphragm=fixed[0]
for s in fixed[1:]:diaphragm=diaphragm+s
source=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');pivot=source.bounds.mean(0);initial_center=pivot*fit['scale']+fit['translation_m'];shape=solid(Mesh(source.vertices-pivot,source.faces,frame))
trials=[];accepted=None;started=time.monotonic()
class Found(Exception):pass
def objective(x):
 global accepted
 s=float(x[0]);delta=np.array(x[1:]);overlap=shape.scale([s,s,s]).translate(initial_center+delta)^diaphragm;assert overlap.status()==md.Error.NoError
 row=dict(scale=s,deltaMetres=delta.tolist(),intersectionVolumeMetresCubed=float(overlap.volume()),empty=overlap.is_empty());trials.append(row)
 if row['empty']:
  accepted=row;raise Found()
 return row['intersectionVolumeMetresCubed']*1e9+np.sum((delta/.01)**2)*.01
try:
 objective([.64,0,-.004,-.001])
 differential_evolution(objective,[(.62,fit['scale']),(-.008,.008),(-.020,.002),(-.008,.004)],seed=204,popsize=5,maxiter=12,polish=False,workers=1,tol=1e-5)
except Found:pass
record=dict(passed=accepted is not None,chosen=accepted,trialCount=len(trials),trials=trials,rotationDegrees=0,fixedSkeleton=True,fixedDiaphragm=True,scaleBounds=[.62,fit['scale']],translationBoundsMetres=[[-.008,.008],[-.020,.002],[-.008,.004]],elapsedSeconds=time.monotonic()-started)
atomic_json(out/'clearance-solver.json',record)
if accepted is None:raise ValueError('Bounded scale/translation fit did not clear the fixed diaphragm')
new=dict(fit);new['scale']=accepted['scale'];new['translation_m']=(initial_center+accepted['deltaMetres']-pivot*accepted['scale']).tolist();new['matrix']=(np.eye(4)).tolist()
for i in range(3):new['matrix'][i][i]=new['scale'];new['matrix'][i][3]=new['translation_m'][i]
new['diaphragmClearance']=accepted;new['refinedFromFit']=str(out/'fit.json');new['anatomicalAcceptance']=False
atomic_json(out/'fit-cleared-v204.json',new)
print(json.dumps({k:new[k]for k in ['scale','translation_m','diaphragmClearance']}),flush=True)
