"""Source-landmark functional partitions with a single reciprocal cut bank.

Functional portal fissures are planes, not sculpted grooves. This candidate
uses the fixed gallbladder/caval and native segment interfaces as controls.
"""
from pathlib import Path
import json,sys,time
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.algorithms.shared_tissue_axis_split import split_shared_cells_at_axis
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
from geometry_tools.tissue import pack

start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'functional-planes-v149';out.mkdir(exist_ok=False)
src=base/'closed-liver-cells-v102';z=np.load(src/'arrangement.npz');meta=json.loads((src/'assembly.json').read_text())
v=z['positions'];surfaces={sid:z[key]for sid,key in meta['surfaces'].items()};gross={'right','left','quadrate'};cells=defaultdict(list)
for name,refs in meta['cells'].items():
    target='body'if name in gross else'I'if name=='caudate'else name
    for r in refs:
        peer=r.get('neighbor')
        if name in gross and peer in gross:continue
        cells[target].append(dict(surface=r['surface'],side=r['side'],neighbor='body'if peer in gross else'I'if peer=='caudate'else peer))
used={r['surface']for refs in cells.values()for r in refs};surfaces={sid:f for sid,f in surfaces.items()if sid in used}
atlas=json.loads((base/'liver-surface-controls-v143/with-landmark.json').read_text());reader=AtlasReader();parts={p['id']:p for p in atlas['parts']}
capsule=reader.read_part(atlas,parts['TA98PART:ZA:Liver:capsule'],frame='model:local-ta98:assembled-metres')
groups={g['term']:g for g in atlas['liverFeatureIntegration']['termGroups']}
def surface_center(n):return capsule.triangles[groups[f'A05.8.01.{n:03d}']['surface']['faces']].mean((0,1))
gall=surface_center(13);caval=surface_center(9);direction=caval-gall;n=np.cross([0,1,0],direction);n/=np.linalg.norm(n)
if n[0]<0:n=-n
main=(n,float(n@gall))
ven=surface_center(10);teres=surface_center(14);d=ven-teres;n=np.cross([0,1,0],d);n/=np.linalg.norm(n)
if n[0]<0:n=-n
umbilical=(n,float(n@teres))
native={label:Mesh.load(base/f'common-field-v13/{number:03d}.npz')for label,number in [('V',1),('VI',0),('VIII',8),('II',4),('III',2)]}
def shared(a,b):
    x=native[a].vertices;y=native[b].vertices;lookup=set(map(tuple,y));return np.array([q for q in x if tuple(q)in lookup])
q=shared('V','VI');center=q.mean(0);n=np.linalg.svd(q-center,full_matrices=False)[2][-1];n[1]=0;n/=np.linalg.norm(n)
if n[2]<0:n=-n
right=(n,float(n@center))
q=shared('V','VIII');center=q.mean(0);n=np.linalg.svd(q-center,full_matrices=False)[2][-1]
if n[1]<0:n=-n
portal=(n,float(n@center))
left_level=float((native['II'].vertices[np.unique(native['II'].faces)].mean(0)[1]+native['III'].vertices[np.unique(native['III'].faces)].mean(0)[1])/2)
planes=[('main',{'body':['R','L']},main),('umbilical',{'L':['IV','lateral-left']},umbilical),
 ('right',{'R':['lateral-right','medial-right']},right),('left-level',{'lateral-left':['III','II']},(np.array([0.,1.,0.]),left_level)),
 ('right-portal-level',{'medial-right':['V','VIII'],'lateral-right':['VI','VII']},portal)]
evidence=[]
for name,split,(n,offset)in planes:
    result=split_shared_cells_at_axis(v,surfaces,cells,split,normal=n,coordinate=offset)
    v,surfaces,cells=result['positions'],result['surfaces'],result['cells']
    evidence.append(dict(name=name,normal=n.tolist(),offset=offset,evidence=result['evidence']))
    print(name,'split',flush=True)
assert set(cells)=={'I','II','III','IV','V','VI','VII','VIII','capsule'}
# Optional functional presentation owns the parenchyma; the ordinary capsule
# retains its existing bank and suppresses its adjacent functional boundary.
used={r['surface']for refs in cells.values()for r in refs};surfaces={sid:f for sid,f in surfaces.items()if sid in used}
mapped,_,conform=conform_quantized_surfaces(v,surfaces,[0,1.25,0]);surfaces=mapped
ids={label:'TA98LIVER:male:segment:'+label if label!='capsule'else'TA98PART:ZA:Liver:capsule'for label in cells}
result=pack(v,surfaces,cells,ids,'model:local-ta98:assembled-metres',out/'packed',origin=[0,1.25,0],max_pairs=2000000,max_tests=16000000,checkpoint_directory=out/'checkpoints')
atomic_json(out/'evidence.json',dict(passed=result['passed'],planes=evidence,conform=conform,sourceCatalogueSha256=digest_file(base/'liver-surface-controls-v143/with-landmark.json'),
    sourcePositionsLiteral=True,sourceBoundary='Repaired capsule-facing liver skin; existing caudate territory retained',
    functionalBoundariesProvisional=True,controls='Actual gallbladder fossa/caval groove, venosum/teres fissures, native V/VI and V/VIII interfaces; II/III source level',
    elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],bytes=result['total_bytes'],elapsedSeconds=time.monotonic()-start)),flush=True)
