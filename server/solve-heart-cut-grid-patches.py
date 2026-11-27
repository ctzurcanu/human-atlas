"""Small bounded integer-grid solves for the remaining coupled cut patches."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.bounded_float32 import untangle
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--previous',default='pv-root-conditioned-complex-v164');args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
bank=h/'pv-root-bank-v159';source=h/'pv-root-encoded-v158';previous=h/args.previous
e=json.loads((bank/'emission.json').read_text());z=np.load(bank/'arrangement.npz');anchor=z['positions'];o=np.array(e['origin']);f=np.concatenate([z[k]for k in e['surfaces'].values()])
p=np.load(previous/'arrangement.npz')['positions'];pr=json.loads((previous/('grouped-rounding.json'if(previous/'grouped-rounding.json').exists()else'construction.json')).read_text());pairs=np.array(pr.get('remaining_declared_defect_pairs',pr.get('global_guard',{}).get('remaining_global_defect_pairs')),int)
raw=np.load(h/'pv-root-interfaces-v157/arrangement.npz');word,first=np.unique((raw['positions']-o).astype('<f4'),axis=0,return_index=True)
assert np.array_equal(word.astype(float)+o,anchor);reference=raw['positions'][first]
base=np.load(h/'regional-bank-v102/arrangement.npz')['positions'];roots=[Mesh.load(h/f'pv-root-conditioned-source-v156/root-{i}.npz').vertices for i in range(4)]
dtype=np.dtype((np.void,12));fixed=np.flatnonzero(np.isin((anchor-o).astype('<f4').view(dtype).ravel(),(np.r_[base,*roots]-o).astype('<f4').view(dtype).ravel()));locked=set(map(int,fixed))
guard=LocalCrossingGuard(p,f,pairs,maximum_displacement_m=4e-8,maximum_pair_tests=3000000)
graph={}
for pair in pairs:
    vertices=set(map(int,f[pair].ravel()))-locked
    for v in vertices:graph.setdefault(v,set()).update(vertices-{v})
seen=set();components=[]
for v in graph:
    if v in seen:continue
    group={v};queue=[v]
    for a in queue:
        for b in graph[a]-group:group.add(b);queue.append(b)
    seen.update(group);components.append(sorted(group))
rows=[];bound=2e-8
print('Remaining cut patches',len(components),'; global crossings',len(pairs),flush=True)
for index,group in enumerate(components):
    if time.monotonic()-start>50.:break
    if len(group)>32:
        rows.append(dict(patch=index,passed=False,error='More than 32 mutable controls'));continue
    remaining=np.array(sorted(guard.remaining),int).reshape(-1,2)
    affected_pairs=remaining[np.isin(f[remaining],group).any(axis=(1,2))]
    if not len(affected_pairs):continue
    patch_deadline=min(start+50.,time.monotonic()+6.);extra_support=[];sat_pairs=affected_pairs.copy()
    try:
        solved=False
        for refinement in range(24):
            if time.monotonic()>=patch_deadline:raise TimeoutError('Finite patch refinement deadline')
            neighborhood=guard.pairs_for(group)
            incident=np.unique(np.concatenate([guard.incident[guard.starts[v]:guard.starts[v+1]]for v in group]))
            face_ids=np.unique(np.r_[neighborhood.ravel(),incident,sat_pairs.ravel()]);global_faces=f[face_ids]
            used,local=np.unique(global_faces,return_inverse=True);lf=local.reshape(-1,3)
            lookup=np.full(len(p),-1,int);lookup[used]=np.arange(len(used));local_group=lookup[group]
            target=reference[used];current=guard.positions[used]
            tri=target[lf];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];words=np.repeat(normal[:,None],3,axis=1)
            support=[]
            for v,g in zip(group,local_group):
                for axis in range(3):
                    for sign in [-1,1]:
                        support.append(dict(control=int(g),normal=(np.eye(3)[axis]*sign).tolist(),point=anchor[v].tolist(),minimum=-bound))
            support.extend(dict(control=int(lookup[v]),normal=n,point=anchor[v].tolist(),minimum=-bound)for v,n in extra_support)
            candidate,evidence=untangle(target,current,lf,lookup[f[sat_pairs]],words,
                origin=o,mutable_controls=local_group,new_steps=32,maximum_displacement_m=bound,
                time_limit=min(2.5,patch_deadline-time.monotonic()),feasibility_only=False,maximum_rounds=6,point_halfspaces=support)
            values=candidate[local_group];delta=values-anchor[group];distance=np.linalg.norm(delta,axis=1)
            if np.any(distance>bound):
                extra_support.extend((int(v),(-d/length).tolist())for v,d,length in zip(group,delta,distance)if length>bound)
                continue
            trial=guard.positions.copy();trial[group]=values
            assert len(np.unique(trial,axis=0))==len(trial), 'New literal coordinate alias'
            proposal=guard.check_move(group,values)
            if proposal is None:
                new=np.array(guard.last_new_pairs,int).reshape(-1,2)
                assert len(new), 'Global crossing counterexample required'
                sat_pairs=np.unique(np.sort(np.r_[sat_pairs,new],axis=1),axis=0)
                expanded=sorted(set(group)|(set(map(int,f[new].ravel()))-locked))
                if len(expanded)>32:raise ValueError('Counterexample requires more than 32 mutable controls')
                group=expanded;continue
            solved=True;break
        if not solved:raise ValueError('Finite patch constraint refinement budget')
        before=len(guard.remaining);guard.commit(proposal)
        row=dict(patch=index,passed=True,controls=len(group),crossings_before=before,crossings_after=len(guard.remaining),refinements=refinement,solve=evidence)
    except (ValueError,TimeoutError,AssertionError)as error:
        row=dict(patch=index,passed=False,controls=len(group),error=str(error))
    rows.append(row);print(json.dumps({k:v for k,v in row.items()if k!='solve'}),flush=True)
    atomic_json(out/'progress.json',dict(patches=rows,remaining_global_crossings=len(guard.remaining),elapsed_seconds=time.monotonic()-start),overwrite=True)
orig=np.load(source/'arrangement.npz');arrays={k:orig[k]for k in orig.files};arrays['positions']=guard.positions
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',json.loads((source/'arrangement.json').read_text()))
proof=dict(passed=not guard.remaining,patches=rows,global_guard=guard.evidence(),
    source_controls_literal=np.array_equal(guard.positions[fixed],anchor[fixed]),
    maximum_total_displacement_m=float(np.linalg.norm(guard.positions-anchor,axis=1).max()),
    elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'construction.json',proof);print(json.dumps({k:v for k,v in proof.items()if k not in ['patches','global_guard']}),flush=True)
