"""Bounded new-cut controls only; original source controls remain literal."""
from pathlib import Path
from itertools import combinations,product
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.grouped_float32_repair import repair_grouped_rounding
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--source-stage',default='coupled-attachments-v249');ap.add_argument('--reference-stage');ap.add_argument('--target-faces',nargs=2,type=int);ap.add_argument('--round-budget',type=int,default=48);ap.add_argument('--evaluation-budget',type=int,default=20000);args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.input;out=h/args.stage;out.mkdir(exist_ok=False);meta=json.loads((src/'ownership.json').read_text());z=np.load(src/'arrangement.npz');current=z['positions'];reference=h/args.reference_stage/'arrangement.npz' if args.reference_stage else src/'arrangement.npz';p=np.load(reference)['positions'];assert p.shape==current.shape;q=current.astype('<f4').astype(float);f=np.concatenate([z[r['faces']]for r in meta['surfaces']]);w=np.concatenate([z[r['normals']]for r in meta['surfaces']]);initial=np.load(h/args.source_stage/'source-arrangement.npz')['positions'];assert np.array_equal(p[:len(initial)],initial);locked=np.arange(len(initial));b=q[f];n=np.cross(b[:,1]-b[:,0],b[:,2]-b[:,0]);bad=np.flatnonzero(np.any(np.einsum('fki,fi->fk',w.astype(float),n)<=0,axis=1));groups=set()
selected=bad[slice(*args.target_faces)] if args.target_faces else bad
for i in selected:
 ids=[int(v)for v in f[i]if v>=len(initial)]
 for count in range(1,len(ids)+1):groups.update(combinations(ids,count))
print('Stored normal defects',len(bad),'targeted',len(selected),'groups',len(groups),flush=True)
quantum=float(np.abs(np.spacing(np.max(np.abs(q.astype('<f4')),axis=0))).max());deltas=np.array([v for v in product(range(-2,3),repeat=3)if any(v)],float)*quantum
result=repair_grouped_rounding(p,f,w,sorted(groups,key=lambda v:(len(v),v)),origin=np.zeros(3),maximum_displacement_m=2e-8,locked=locked,translations_m=deltas,maximum_group_size=3,max_evaluations=args.evaluation_budget,maximum_rounds=args.round_budget,wall_seconds=20.,stored_positions=q)
arrays={k:z[k]for k in z.files};arrays['positions']=result['positions'];atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'ownership.json',meta);proof=result['evidence'];proof['immutable_reference']=str(reference);proof['target_faces']=args.target_faces;proof['original_source_controls_literal']=np.array_equal(result['positions'][locked],initial.astype('<f4').astype(float));proof['full_global_geometry_checks_required']=True;atomic_json(out/'rounding.json',proof);print(json.dumps({k:proof[k]for k in ['passed','initial_bad_faces','remaining_bad_faces','maximum_displacement_m','budget_or_stop_reason','elapsed_seconds']}),flush=True)
