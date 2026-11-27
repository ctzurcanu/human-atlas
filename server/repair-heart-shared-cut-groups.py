"""Resolve coupled stored cut points under the same immutable source bound."""
from pathlib import Path
from itertools import combinations,product
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.grouped_float32_repair import repair_grouped_rounding
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--previous',default='pv-root-conditioned-complex-v162');args=ap.parse_args()
start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
bank=h/'pv-root-bank-v159';source=h/'pv-root-encoded-v158';previous=h/args.previous
e=json.loads((bank/'emission.json').read_text());z=np.load(bank/'arrangement.npz');anchor=z['positions'];o=np.array(e['origin'])
f=np.concatenate([z[key]for key in e['surfaces'].values()]);p=np.load(previous/'arrangement.npz')['positions']
report_path=previous/('grouped-rounding.json'if(previous/'grouped-rounding.json').exists()else'rounding.json')
report=json.loads(report_path.read_text());pairs=np.array(report.get('remaining_global_defect_pairs',report.get('global_crossing_guard',{}).get('remaining_global_defect_pairs')),int)
base=np.load(h/'regional-bank-v102/arrangement.npz')['positions'];roots=[Mesh.load(h/f'pv-root-conditioned-source-v156/root-{i}.npz').vertices for i in range(4)]
dtype=np.dtype((np.void,12));literal=np.r_[base,*roots];fixed=np.flatnonzero(np.isin((anchor-o).astype('<f4').view(dtype).ravel(),(literal-o).astype('<f4').view(dtype).ravel()));locked=set(map(int,fixed))
groups=set();graph={}
for a,b in pairs:
    vertices=set(map(int,f[[a,b]].ravel()))-locked
    for subset in [vertices,*[set(map(int,f[i]))-locked for i in [a,b]]]:
        if 1<=len(subset)<=16:groups.add(tuple(sorted(subset)))
    for v in vertices:graph.setdefault(v,set()).update(vertices-{v})
    for pair in combinations(sorted(vertices),2):groups.add(pair)
seen=set()
for v in graph:
    if v in seen:continue
    component={v};queue=[v]
    for a in queue:
        for b in graph[a]-component:component.add(b);queue.append(b)
    seen.update(component)
    if 1<=len(component)<=16:groups.add(tuple(sorted(component)))
groups=sorted(groups,key=lambda g:(len(g),g));assert len(groups)<=2000
q=anchor[f];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
words=np.repeat(normal[:,None],3,axis=1)
quantum=float(np.max(abs(np.spacing(np.max(abs((anchor-o).astype('<f4')),axis=0)))))
translations=np.array([v for v in product(range(-2,3),repeat=3)if any(v)],float)*quantum
guard=LocalCrossingGuard(p,f,pairs,maximum_displacement_m=2e-8,maximum_pair_tests=6000000)
print('Coupled cut groups',len(groups),'; remaining global crossings',len(pairs),flush=True)
r=repair_grouped_rounding(anchor,f,words,groups,origin=o,stored_positions=p,
    maximum_displacement_m=2e-8,locked=fixed,defect_face_pairs=pairs,
    translations_m=translations,global_guard=guard,maximum_group_size=16,
    max_evaluations=40000,maximum_rounds=128,wall_seconds=40.)
atomic_json(out/'grouped-rounding.json',r['evidence'])
orig=np.load(source/'arrangement.npz');arrays={k:orig[k]for k in orig.files};arrays['positions']=r['positions']
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',json.loads((source/'arrangement.json').read_text()))
proof=dict(passed=r['evidence']['passed'],remaining_global_crossings=len(guard.remaining),
    source_controls_literal=np.array_equal(r['positions'][fixed],anchor[fixed]),
    maximum_total_displacement_m=float(np.linalg.norm(r['positions']-anchor,axis=1).max()),
    previous_rounding_sha256=digest_file(report_path),
    elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'construction.json',proof);print(json.dumps(proof),flush=True)
