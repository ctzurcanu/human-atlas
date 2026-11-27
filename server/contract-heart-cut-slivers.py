"""Contract only diagnosed sub-storage cut edges, testing every material link."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.shared_edge_collapse import collapse_short_edges
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--previous',default='pv-root-conditioned-complex-v167');ap.add_argument('--rings',type=int,default=1,choices=[0,1,2]);ap.add_argument('--source-parent-fields',action='store_true');args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);previous=h/args.previous;bank=h/'pv-root-bank-v159'
meta=json.loads((previous/'arrangement.json').read_text());z=np.load(previous/'arrangement.npz');anchor=np.load(bank/'arrangement.npz')['positions'];p=z['positions'];o=np.array(meta['origin'])
rows=meta['surfaces'];row_by_id={r['id']:r for r in rows};emission=json.loads((bank/'emission.json').read_text());bz=np.load(bank/'arrangement.npz')
surfaces={sid:z[row_by_id[sid]['array']]for sid in sorted(row_by_id)};f=np.concatenate(list(surfaces.values()))
if(previous/'diagnostics.json').exists():
    proof=json.loads((previous/'diagnostics.json').read_text());assert proof['geometry_sha256']==Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha;pairs=proof['remaining_global_defect_pairs']
else:
    assert np.array_equal(f,np.concatenate([bz[key]for key in emission['surfaces'].values()])), 'Diagnostic face order must match the literal bank'
    proof=json.loads((previous/'construction.json').read_text());pairs=proof['global_guard']['remaining_global_defect_pairs']
base=np.load(h/'regional-bank-v102/arrangement.npz')['positions'];roots=[Mesh.load(h/f'pv-root-conditioned-source-v156/root-{i}.npz').vertices for i in range(4)]
dtype=np.dtype((np.void,12));fixed=np.flatnonzero(np.isin((anchor-o).astype('<f4').view(dtype).ravel(),(np.r_[base,*roots]-o).astype('<f4').view(dtype).ravel()))
fields={}
if args.source_parent_fields:
    old=json.loads((h/'regional-bank-v102/emission.json').read_text());oz=np.load(h/'regional-bank-v102/arrangement.npz');op=(oz['positions']-o).astype('<f4').astype(float)+o;rv,rf=z['root_positions'],z['root_faces']
for sid,faces in surfaces.items():
    if args.source_parent_fields:
        row=row_by_id[sid];pa=z[row['parents']];q=rv[rf[pa]]if row.get('parentSurface')=='occupied-outer'else op[oz[old['surfaces'][sid]][pa]];sign=row.get('parentNormalSign',1)
    else:q=p[faces];sign=1
    n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])*sign;n/=np.linalg.norm(n,axis=1)[:,None];fields[sid]=np.repeat(n[:,None],3,axis=1)
guard=LocalCrossingGuard(p,f,pairs,maximum_displacement_m=4e-8,maximum_pair_tests=4000000)
r=collapse_short_edges(p,surfaces,meta['cells'],fields,desired_positions=anchor,
    maximum_displacement_m=2e-8,maximum_steps=48,wall_seconds=20.,candidate_face_rings=args.rings,
    defect_face_pairs=pairs,preserve_material_membership=False,locked=fixed,global_guard=guard)
atomic_json(out/'contractions.json',r['evidence']);arrays={k:z[k]for k in z.files};arrays['positions']=p
new_rows=[]
for row in rows:
    faces=r['surfaces'][row['id']]
    if not len(faces):continue
    arrays[row['array']]=faces;arrays[row['parents']]=z[row['parents']][r['source_faces'][row['id']]]
    new_rows.append(row)
alive={r['id']for r in new_rows};cells={name:[u for u in uses if u['surface']in alive]for name,uses in meta['cells'].items()};assert all(cells.values())
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',dict(meta,surfaces=new_rows,cells=cells))
newf=np.concatenate([r['surfaces'][sid]for sid in surfaces if len(r['surfaces'][sid])]);assert np.array_equal(newf,guard.faces[guard.active])
lookup=np.full(len(f),-1,int);lookup[guard.active]=np.arange(guard.active.sum());remaining=np.array(sorted(guard.remaining),int).reshape(-1,2);mapped=lookup[remaining];assert np.all(mapped>=0)
atomic_json(out/'diagnostics.json',dict(geometry_sha256=Mesh(p,newf,'model:local-ta98:assembled-metres').geometry_sha,remaining_global_defect_pairs=mapped.tolist(),complete_guarded_worklist=True))
report=dict(passed=r['evidence']['passed']and not guard.remaining,
    remaining_global_crossings=len(guard.remaining),source_controls_literal=True,
    removed_physical_faces=r['evidence']['removed_physical_faces'],
    contractions=len(r['evidence']['contractions']),stop=r['evidence'].get('budget_or_stop_reason'),
    elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
