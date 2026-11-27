"""Repair diagnosed encoded cut crossings with an immutable full-chart BVH."""
from pathlib import Path
import argparse,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.bounded_float32_rounding import round_diagnosed_crossings
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--bank',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--source-stage',default='coupled-attachments-v249');ap.add_argument('--reference-stage',default='coupled-attachments-v249');args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.input;bank=h/args.bank;out=h/args.stage;out.mkdir(exist_ok=False);meta=json.loads((src/'ownership.json').read_text());z=np.load(src/'arrangement.npz');p=z['positions'];rows=meta['surfaces'];f=np.concatenate([z[r['faces']]for r in rows]);offsets=dict(zip([r['id']for r in rows],np.r_[0,np.cumsum([len(z[r['faces']])for r in rows])][:-1]));e=json.loads((bank/'emission.json').read_text());assert digest_file(bank/'cells.bin')==e['bank_sha256'];reader=AtlasReader();mini=dict(chunks=[dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')))]);cached={}
for path in(bank/'checkpoints').glob('*.json'):
 r=json.loads(path.read_text())
 if r['kind']=='cell-audit':cached[r['geometry'][0]]=r
pairs=set();failed=[]
for pid,layout in e['parts'].items():
 m=reader.read_part(mini,dict(layout,chunk=0,geometryOrigin=[0.,0.,0.]),frame='heart-native-local-metres');record=cached.get(m.geometry_sha)
 if not record:raise ValueError('Finish the current cell audit first: '+pid)
 proof=record['result'];si=proof['self_intersections'];assert si['complete_within_pair_budget']and len(si['examples'])==si['defect_pairs'];ids=np.concatenate([offsets[u['surface']]+np.arange(len(z[next(r['faces']for r in rows if r['id']==u['surface'])]))for u in meta['cells'][pid]])
 if not proof['passed']:
  assert proof['failures']==['self_intersection_check'];failed.append(pid)
 for a,b in si['examples']:pairs.add(tuple(sorted([int(ids[a]),int(ids[b])])))
anchor=np.load(h/args.reference_stage/'arrangement.npz')['positions'];fixed=np.arange(len(np.load(h/args.source_stage/'source-arrangement.npz')['positions']));print('Diagnosed chart crossings',len(pairs),'in',failed,flush=True)
r=round_diagnosed_crossings(p,f,sorted(pairs),origin=np.zeros(3),protected=fixed,maximum_displacement_m=2e-8,anchor_positions=anchor,seconds=20.,maximum_moves=32,maximum_pair_tests=6000000);arrays={k:z[k]for k in z.files};arrays['positions']=r['positions'];atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'ownership.json',meta);atomic_json(out/'crossing-rounding.json',r['evidence']);print(json.dumps({k:r['evidence'][k]for k in ['passed','stop_reason','maximum_total_anchor_displacement_m','remaining_global_defect_pair_count','elapsed_seconds']}),flush=True)
