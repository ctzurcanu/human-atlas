"""Condition only the new reciprocal cut controls, under a complete global guard."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.algorithms.bounded_float32_rounding import round_diagnosed_crossings
ap=argparse.ArgumentParser();ap.add_argument('--bank',required=True);ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--resume')
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';bank=h/args.bank;source=h/args.source;out=h/args.stage;out.mkdir(exist_ok=False)
e=json.loads((bank/'emission.json').read_text());z=np.load(bank/'arrangement.npz');p=z['positions'];f=np.concatenate([z[key]for key in e['surfaces'].values()]);m=Mesh(p,f,e['frame']);plan=json.loads((bank/'physical-scan-plan.json').read_text());assert plan['geometry_sha256']==m.geometry_sha
pairs=[]
for i in range(3):
    r=json.loads((bank/f'physical-scan-{i}.json').read_text())
    assert r['geometry_sha256']==m.geometry_sha and r['root_pair']==plan['pairs'][i]
    assert r['complete_within_pair_budget'] and len(r['examples'])==r['defect_pairs']
    pairs.extend(r['examples'])
o=np.array(e['origin']);dtype=np.dtype((np.void,12));base=np.load(h/'regional-bank-v102/arrangement.npz')['positions'];root=[Mesh.load(h/f'pv-root-conditioned-source-v156/root-{i}.npz').vertices for i in range(4)]
literal=np.r_[base,*root];protected=np.flatnonzero(np.isin((p-o).astype('<f4').view(dtype).ravel(),(literal-o).astype('<f4').view(dtype).ravel()))
anchor=p
if args.resume:
    previous=h/args.resume;pr=json.loads((previous/'rounding.json').read_text());saved=np.load(previous/'arrangement.npz')['positions']
    replay=p.copy()
    for move in pr['moves']:
        v=move['vertex'];assert np.array_equal(replay[v],move['from_m']);replay[v]=move['to_m']
    assert np.array_equal(replay,saved), 'Saved rounding coordinates must replay the recorded moves exactly'
    p=saved;pairs=pr['remaining_global_defect_pairs']
print('Complete global defect list',len(pairs),'; fixed source controls',len(protected),flush=True)
r=round_diagnosed_crossings(p,f,pairs,origin=o,protected=protected,seconds=40.,maximum_moves=256,
    maximum_pair_tests=12000000,maximum_displacement_m=2e-8,anchor_positions=anchor)
atomic_json(out/'rounding.json',r['evidence'])
original=np.load(source/'arrangement.npz');assert np.array_equal(original['positions'],anchor)
arrays={k:original[k]for k in original.files};arrays['positions']=r['positions'];atomic_npz(out/'arrangement.npz',**arrays)
atomic_json(out/'arrangement.json',json.loads((source/'arrangement.json').read_text()))
atomic_json(out/'construction.json',dict(passed=r['evidence']['passed'],
    literal_source_controls_fixed=r['evidence']['protected_controls_literal'],
    new_literal_aliases=len(np.unique(r['positions'],axis=0))!=len(p),
    elapsed_seconds=time.monotonic()-start,installed=False))
print(json.dumps({k:r['evidence'][k]for k in ['passed','stop_reason','remaining_global_defect_pair_count','maximum_total_anchor_displacement_m','accepted_moves','elapsed_seconds']}),flush=True)
