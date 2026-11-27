"""Retriangulate one source-parent chart without changing its boundary controls."""
from pathlib import Path
import argparse,json,sys,time
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);source=h/'pv-root-conditioned-complex-v174'
meta=json.loads((source/'arrangement.json').read_text());z=np.load(source/'arrangement.npz');p=z['positions'];o=np.array(meta['origin']);row_by_id={r['id']:r for r in meta['surfaces']};names=sorted(row_by_id);counts=[len(z[row_by_id[s]['array']])for s in names];offsets=np.r_[0,np.cumsum(counts)];f=np.concatenate([z[row_by_id[s]['array']]for s in names]);parents=np.concatenate([z[row_by_id[s]['parents']]for s in names]);owners=np.repeat(np.arange(len(names)),counts)
e=json.loads((h/'regional-bank-v102/emission.json').read_text());old=np.load(h/'regional-bank-v102/arrangement.npz');op=(old['positions']-o).astype('<f4').astype(float)+o;rv,rf=z['root_positions'],z['root_faces'];normals=[]
for sid in names:
    row=row_by_id[sid];pa=z[row['parents']];src=rv[rf[pa]]if row.get('parentSurface')=='occupied-outer'else op[old[e['surfaces'][sid]][pa]];n=np.cross(src[:,1]-src[:,0],src[:,2]-src[:,0])*row.get('parentNormalSign',1);n/=np.linalg.norm(n,axis=1)[:,None];normals.append(n)
n=np.concatenate(normals);diag=json.loads((source/'diagnostics.json').read_text());assert diag['geometry_sha256']==Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha
guard=LocalCrossingGuard(p,f,diag['remaining_global_defect_pairs'],maximum_displacement_m=.005,maximum_pair_tests=6000000)
changes=[]
def boundary(ff):
    c=Counter((int(a),int(b))for t in ff for a,b in zip(t,np.roll(t,-1)))
    return Counter({k:v-c[k[::-1]]for k,v in c.items()if v>c[k[::-1]]})
for iteration in range(24):
    if time.monotonic()-start>20.:break
    q=p[f];geo=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);bad=np.flatnonzero(np.einsum('ij,ij->i',geo,n)<=0)
    if not len(bad):break
    accepted=False
    for i in bad:
        sid=int(owners[i]);lo,hi=offsets[sid:sid+2];local=f[lo:hi]
        for a,b,c in [f[i],np.roll(f[i],-1),np.roll(f[i],-2)]:
            peers=np.flatnonzero(np.isin(local,[a,b]).sum(1)==2)+lo
            for j in peers:
                if i==j or parents[i]!=parents[j]:continue
                other=f[j];d=next((int(v)for v in other if v not in [a,b]),None)
                if d is None or d==c:continue
                if np.any((f==c).any(1)&(f==d).any(1)):continue
                proposal=np.array([[c,a,d],[c,d,b]],int)
                if boundary(f[[i,j]])!=boundary(proposal):continue
                pp=p[proposal];nn=np.cross(pp[:,1]-pp[:,0],pp[:,2]-pp[:,0])
                if np.any(np.einsum('ij,ij->i',nn,n[[i,j]])<=0):continue
                # Both replacements stay in the same declared original
                # source triangle chart. Its four boundary controls stay
                # literal; the broad BVH bound is for changed diagonals only.
                guarded=guard.check_faces([i,j],proposal,np.ones(2,bool))
                if guarded is None:continue
                guard.commit(guarded);f[[i,j]]=proposal;changes.append(dict(faces=[int(i),int(j)],source_chart=names[sid],source_parent=int(parents[i])));accepted=True;break
            if accepted:break
        if accepted:break
    if not accepted:break
q=p[f];remaining_bad=np.flatnonzero(np.einsum('ij,ij->i',np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),n)<=0)
arrays={k:z[k]for k in z.files}
for i,sid in enumerate(names):arrays[row_by_id[sid]['array']]=f[offsets[i]:offsets[i+1]]
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',meta)
proof=dict(passed=not len(remaining_bad)and not guard.remaining,reversed_source_facets=remaining_bad.tolist(),remaining_global_crossings=len(guard.remaining),diagonals=changes,all_positions_literal=True,elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'construction.json',proof);atomic_json(out/'diagnostics.json',dict(geometry_sha256=Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha,remaining_global_defect_pairs=sorted(guard.remaining),complete_guarded_worklist=True));print(json.dumps(proof),flush=True)
