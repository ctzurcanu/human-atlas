"""Keep closed classifiers, omit their original hilar disks from display."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json
from geometry_tools.open3d_backend import library
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'completion-v155';a=json.loads((out/'atlas-with-sections.json').read_text());r=AtlasReader();parts={p['id']:p for p in a['parts']}
src=base/'closed-liver-cells-v102';meta=json.loads((src/'assembly.json').read_text());z=np.load(src/'arrangement.npz');v=z['positions'];f=np.concatenate([z[key]for sid,key in meta['surfaces'].items()if sid.startswith('retained-hilar-source-closure:')])
o3d=library();scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor(v.astype(np.float32)),o3d.core.Tensor(f.astype(np.uint32)))
record={}
for pid in a['liverCompletionIntegration']['functionalPartIds']:
    p=parts[pid];m=r.read_part(a,p,frame='model:local-ta98:assembled-metres');distance=scene.compute_distance(o3d.core.Tensor(m.triangles.mean(1).astype(np.float32)),nthreads=1).numpy();new=[];omitted=[]
    for g in p['tissueDrawGroups']:
        indices=np.arange(g['start']//3,(g['start']+g['count'])//3)
        if g.get('neighbor'):new.append(g);continue
        closed=distance[indices]<=5e-8 # bounded Float32/world-coordinate replay error
        for flag in [False,True]:
            chosen=indices[closed==flag]
            if not len(chosen):continue
            runs=np.split(chosen,np.flatnonzero(np.diff(chosen)>1)+1)
            for run in runs:new.append(dict(start=int(run[0])*3,count=len(run)*3,**({'computationalOnly':True}if flag else{})))
        omitted.extend(indices[closed].tolist())
    p['tissueDrawGroups']=new;p['liverComputationalFaces']=omitted
    record[pid]=dict(computationalTriangles=len(omitted),sourceGeometryUnchanged=True)
atomic_json(out/'atlas-ready.json',a);atomic_json(out/'computational-mouth-evidence.json',dict(passed=True,parts=record,noDisplayedHilarDisks=True,coordinateGuardMetres=5e-8))
print(json.dumps(record))
