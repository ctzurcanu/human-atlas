from pathlib import Path
import json,runpy,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from ta98_glb_geometry import load_roots
mapping=runpy.run_path(str(ROOT/'server/heart-reference-types.py'))['SOURCE_CIDS']
metadata=json.loads((ROOT.parent/'3d/ext/meshes/964db2dd4f98052f03baa9ca5f2dbcae/assets.json').read_text())
rows={r['sni'][0]['id']:r for r in metadata['structures']if r.get('cid')in set(mapping.values())}
sources=load_roots(ROOT/'.local-models/models/Male_2024.glb',set(rows),allow_zero_weight_targets=True)
matrix=np.asarray(json.loads((ROOT/'.local-models/geometry-workbench/heart/cardiac-label-landmarks-v604/report.json').read_text())['matrix'])
result=[];started=time.monotonic()
for m in sources:
    v=m['unique']@matrix[:3,:3].T+matrix[:3,3];f=m['welded'];cases=[]
    for eps in [0.,1e-8,1e-7,1e-6]:
        _,remap=np.unique(v if eps==0 else np.rint(v/eps),axis=0,return_inverse=True)
        ff=remap[f];ff=ff[(ff[:,0]!=ff[:,1])&(ff[:,1]!=ff[:,2])&(ff[:,0]!=ff[:,2])]
        _,unique=np.unique(np.sort(ff,axis=1),axis=0,return_index=True);ff=ff[unique]
        edges=np.sort(np.concatenate([ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]]),axis=1);_,counts=np.unique(edges,axis=0,return_counts=True)
        cases.append(dict(weld_metres=eps,triangles=len(ff),open_edges=int((counts==1).sum()),nonmanifold_edges=int((counts>2).sum())))
    result.append(dict(cid=rows[m['node']]['cid'],node=m['node'],cases=cases))
path=ROOT/'.local-models/geometry-workbench/heart/typed-source-topology-v645.json';assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
for r in result:print(r['cid'],r['cases'],flush=True)
print('elapsed',round(time.monotonic()-started,2),'s',flush=True)
