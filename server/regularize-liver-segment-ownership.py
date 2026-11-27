"""Transfer an isolated calculation sliver to its largest tissue neighbour.

No position changes, source faces disappear only in reciprocal internal pairs.
The adjusted territory boundary remains explicitly provisional.
"""
from pathlib import Path
from collections import defaultdict
import sys,json,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.tissue import pack
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';src=base/'functional-planes-v149/packed';out=base/'functional-regions-v152';out.mkdir(exist_ok=False)
meta=json.loads((src/'tissue.json').read_text());z=np.load(src/'arrangement.npz');v=z['positions'];source={sid:z[key]for sid,key in meta['surfaces'].items()}
q=(v-[0,1.25,0]).astype('<f4');_,first,inv=np.unique(q,axis=0,return_index=True,return_inverse=True);canon=first[inv]
cellfaces={name:canon[np.concatenate([source[r['surface']]if r['side']==1 else source[r['surface']][:,[0,2,1]]for r in refs])]for name,refs in meta['cells'].items()}
owners=defaultdict(list)
for name,ff in cellfaces.items():
    for f in ff:owners[tuple(sorted(map(int,f)))].append((name,f))
changes=[]
for name,ff in cellfaces.items():
    if name=='capsule':continue
    edgefaces=defaultdict(list)
    for i,f in enumerate(ff):
        for a,b in zip(f,np.roll(f,-1)):edgefaces[tuple(sorted((int(a),int(b))))].append(i)
    pairs=np.array([ids for ids in edgefaces.values()if len(ids)==2],int).reshape(-1,2)
    graph=coo_matrix((np.ones(len(pairs)*2),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(ff),len(ff))).tocsr()
    n,components=connected_components(graph,directed=False)
    if n==1:continue
    counts=np.bincount(components);main=counts.argmax()
    for k in range(n):
        if k==main:continue
        ids=np.flatnonzero(components==k);assert len(ids)<=16,'No automatic transfer of a substantive native region'
        areas=defaultdict(float)
        for f in ff[ids]:
            for peer,g in owners[tuple(sorted(map(int,f)))]:
                if peer not in [name,'capsule']:
                    tri=v[f];areas[peer]+=float(np.linalg.norm(np.cross(tri[1]-tri[0],tri[2]-tri[0]))/2)
        assert areas;target=max(areas,key=areas.get)
        for f in ff[ids]:
            key=tuple(sorted(map(int,f)));owners[key]=[(target if a==name else a,g)for a,g in owners[key]]
        tri=v[ff[ids]];volume=float(abs(np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6))
        changes.append(dict(source=name,target=target,triangles=len(ids),volume_m3=volume,positionsMoved=False))
patches=defaultdict(list)
for uses in owners.values():
    if len(uses)==2 and uses[0][0]==uses[1][0]:continue
    if len(uses)==1:a,f=uses[0];patches[(a,None)].append(f)
    else:(a,f),(b,g)=uses;patches[(a,b)].append(f)
surfaces={};cells={name:[]for name in cellfaces}
for i,((a,b),ff)in enumerate(patches.items()):
    sid='territory:'+str(i);surfaces[sid]=np.array(ff,int);cells[a].append(dict(surface=sid,side=1,neighbor=b))
    if b:cells[b].append(dict(surface=sid,side=-1,neighbor=a))
result=pack(v,surfaces,cells,meta['part_ids'],meta['frame'],out/'packed',origin=[0,1.25,0],max_pairs=2000000,max_tests=16000000,checkpoint_directory=out/'checkpoints')
assert result['passed'],result['checks']['failures']
assert all(len(c['topology']['components'])==1 for c in result['checks']['cells'].values())
atomic_json(out/'evidence.json',dict(passed=True,oneConnectedComponentPerTerritory=True,sourcePositionsLiteral=True,changes=changes,
    reciprocalTissueInterfaces=True,allOccupiedVolumeChecksPassed=True,exteriorUnchanged=True,elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(passed=True,changes=changes,bytes=result['total_bytes'],elapsedSeconds=time.monotonic()-start)),flush=True)
