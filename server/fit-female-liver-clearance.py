"""Native finite-volume search for a small anterior/caudal placement."""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'female-transfer-v202';a=json.loads((out/'atlas.json').read_text());fit=json.loads((out/'fit.json').read_text());reader=AtlasReader();frame='model:local-female-ta98:assembled-metres'
def solid(mesh):
 v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);s=md.Manifold(md.Mesh64(np.ascontiguousarray(v),np.ascontiguousarray(inv[mesh.faces],dtype=np.uint64)))
 assert s.status()==md.Error.NoError,'Closed source tissue required';return s
fixed=[]
for pid in a['finiteDiaphragmIntegration']['partIds']:
 p=next(p for p in a['parts']if p['id']==pid);fixed.append(solid(reader.read_part(a,p,frame=frame)))
diaphragm=fixed[0]
for s in fixed[1:]:diaphragm=diaphragm+s
source=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');v=source.vertices*fit['scale']+fit['translation_m'];envelope=solid(Mesh(v,source.faces,frame))
trials=[];chosen=None;started=time.monotonic()
offsets=[(0,dz)for dz in [-.002,-.004,-.006,-.008,-.010,-.012]]
for dy,dz in offsets:
 delta=[0.,dy,dz];overlap=envelope.translate(delta)^diaphragm;assert overlap.status()==md.Error.NoError
 r=dict(deltaMetres=delta,intersectionVolumeMetresCubed=float(overlap.volume()),triangles=int(overlap.num_tri()),empty=overlap.is_empty());trials.append(r)
 print(json.dumps(r),flush=True)
 if r['empty']:chosen=r;break
record=dict(passed=chosen is not None,chosen=chosen,trials=trials,fixedSkeleton=True,fixedDiaphragm=True,scope='Actual closed liver envelope against the union of all fixed female diaphragm tissue cells',elapsedSeconds=time.monotonic()-started)
atomic_json(out/'finite-volume-clearance-posterior.json',record)
if chosen is None:raise ValueError('No bounded tested translation clears the female diaphragm')
