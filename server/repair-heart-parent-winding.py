"""Repair storage against source-parent orientation, never the damaged embedding."""
from pathlib import Path
from itertools import product
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.grouped_float32_repair import repair_grouped_rounding
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--source',default='pv-root-conditioned-complex-v171');ap.add_argument('--bank',default='pv-root-bank-v159');ap.add_argument('--local-source-adjustments',action='store_true');ap.add_argument('--cut-bound-nm',type=float,default=20.);ap.add_argument('--local-source-bound-nm',type=float,default=10.);args=ap.parse_args();start=time.monotonic()
assert 0<args.cut_bound_nm<=1000 and 0<args.local_source_bound_nm<=100
cut_bound=args.cut_bound_nm*1e-9;source_bound=args.local_source_bound_nm*1e-9
assert source_bound<=cut_bound
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);source=h/args.source;bank=h/args.bank
meta=json.loads((source/'arrangement.json').read_text());z=np.load(source/'arrangement.npz');p=z['positions'];o=np.array(meta['origin']);anchor=np.load(bank/'arrangement.npz')['positions'];rows={r['id']:r for r in meta['surfaces']}
old=json.loads((h/'regional-bank-v102/emission.json').read_text());oz=np.load(h/'regional-bank-v102/arrangement.npz');op=(oz['positions']-o).astype('<f4').astype(float)+o;rv,rf=z['root_positions'],z['root_faces']
faces=[];fields=[]
for sid in sorted(rows):
    r=rows[sid];f=z[r['array']];pa=z[r['parents']]
    src=rv[rf[pa]]if r.get('parentSurface')=='occupied-outer'else op[oz[old['surfaces'][sid]][pa]]
    n=np.cross(src[:,1]-src[:,0],src[:,2]-src[:,0])*r.get('parentNormalSign',1);n/=np.linalg.norm(n,axis=1)[:,None]
    faces.append(f);fields.append(np.repeat(n[:,None],3,axis=1))
f=np.concatenate(faces);words=np.concatenate(fields);q=p[f];actual=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);bad=np.flatnonzero(np.einsum('ij,ij->i',words[:,0],actual)<=0)
diag=json.loads((source/'diagnostics.json').read_text());assert diag['geometry_sha256']==Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha;pairs=np.array(diag['remaining_global_defect_pairs'],int).reshape(-1,2)
roots=[Mesh.load(h/meta.get('root_source_stage','pv-root-conditioned-source-v156')/f'root-{i}.npz').vertices for i in range(4)];typ=np.dtype((np.void,12));fixed=np.flatnonzero(np.isin((anchor-o).astype('<f4').view(typ).ravel(),(np.r_[oz['positions'],*roots]-o).astype('<f4').view(typ).ravel()));locked=set(map(int,fixed))
all_source_fixed=fixed.copy();control_bounds=np.full(len(p),cut_bound)
groups=set();selection=np.unique(np.r_[bad,pairs.ravel()]);graph={}
if args.local_source_adjustments:
    targeted=np.unique(f[selection]);neighbors=np.flatnonzero(np.isin(f,targeted).any(1));allowed=np.intersect1d(fixed,np.unique(f[neighbors]));assert len(allowed)<=512
    control_bounds[allowed]=source_bound;fixed=np.setdiff1d(fixed,allowed);locked=set(map(int,fixed))
for tri in f[selection]:
    mutable=set(map(int,tri))-locked
    for v in mutable:groups.add((v,));graph.setdefault(v,set()).update(mutable-{v})
    if 1<=len(mutable)<=16:groups.add(tuple(sorted(mutable)))
for pair in pairs:
    mutable=set(map(int,f[pair].ravel()))-locked
    if 1<=len(mutable)<=16:groups.add(tuple(sorted(mutable)))
seen=set()
for v in graph:
    if v in seen:continue
    component={v};queue=[v]
    for a in queue:
        for b in graph[a]-component:component.add(b);queue.append(b)
    seen.update(component)
    if 1<=len(component)<=16:groups.add(tuple(sorted(component)))
groups=sorted(groups,key=lambda x:(len(x),x));assert len(groups)<=2000
quantum=float(np.max(abs(np.spacing(np.max(abs((anchor-o).astype('<f4')),axis=0)))))
translations=np.array([v for v in product([-1,-.5,-.25,0,.25,.5,1],repeat=3)if any(v)],float)*quantum
translations=np.r_[translations,np.array([v for v in product([-2,0,2],repeat=3)if any(v)],float)*quantum]
if cut_bound>2e-8:
    for radius in [4,8,16,32]:
        translations=np.r_[translations,np.array([v for v in product([-radius,0,radius],repeat=3)if any(v)],float)*quantum]
assert len(translations)<=500
guard=LocalCrossingGuard(p,f,pairs,maximum_displacement_m=2*cut_bound,maximum_pair_tests=8000000)
print('Source-parent reversals',len(bad),'; crossings',len(pairs),'; candidate groups',len(groups),flush=True)
r=repair_grouped_rounding(anchor,f,words,groups,origin=o,stored_positions=p,
    maximum_displacement_m=cut_bound,locked=fixed,defect_face_pairs=pairs,
    translations_m=translations,global_guard=guard,maximum_group_size=16,
    max_evaluations=90000,maximum_rounds=256,wall_seconds=40.,maximum_displacement_by_vertex=control_bounds)
atomic_json(out/'source-winding-repair.json',r['evidence']);arrays={k:z[k]for k in z.files};arrays['positions']=r['positions'];atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',meta)
source_delta=np.linalg.norm(r['positions'][all_source_fixed]-anchor[all_source_fixed],axis=1)
report=dict(passed=r['evidence']['passed'],reversed_facets=r['evidence']['remaining_bad_faces'],
    remaining_global_crossings=len(guard.remaining),source_controls_literal=not np.any(source_delta),source_controls_changed=int(np.count_nonzero(source_delta)),maximum_source_adjustment_m=float(source_delta.max(initial=0)),
    maximum_adjustment_m=float(np.linalg.norm(r['positions']-anchor,axis=1).max()),
    cut_control_bound_m=cut_bound,source_control_bound_m=source_bound if args.local_source_adjustments else 0.,
    stop=r['evidence']['budget_or_stop_reason'],elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'construction.json',report);atomic_json(out/'diagnostics.json',dict(geometry_sha256=Mesh(r['positions'],f,'model:local-ta98:assembled-metres').geometry_sha,remaining_global_defect_pairs=sorted(guard.remaining),complete_guarded_worklist=True));print(json.dumps(report),flush=True)
