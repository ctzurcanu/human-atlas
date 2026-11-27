"""Collapse diagnosed submicrometre split stars onto literal carrier controls.

This is restricted to the complete remaining numerical defect worklist. It
does not smooth anatomy or weld surfaces by proximity. All carrier controls
stay literal, every affected shared cell is updated and closure is rechecked.
"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
from geometry_tools.audit import audit
from geometry_tools.model import Mesh
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'regional-conditioned-v99';out=b/'regional-clean-v101';out.mkdir(exist_ok=False);z=np.load(source/'cells.npz');p=z['positions'];e=json.loads((b/'regional-bank-v98/emission.json').read_text());layout=np.load(b/'regional-bank-v98/arrangement.npz');f=np.concatenate([layout[k]for k in e['surfaces'].values()]);r=json.loads((source/'rounding.json').read_text());affected=np.unique(f[np.asarray(r['remaining_global_defect_pairs']).ravel()]);base=np.load(b/'chart-coats-v82/cells.npz')['positions'];dist,nearest=cKDTree(base).query(p[affected]);bound=.00000025;lookup={tuple(v):i for i,v in enumerate(p)};remap=np.arange(len(p));moves=[]
for vertex,distance,control in zip(affected,dist,nearest):
 if not 0<distance<=bound:continue
 target=lookup[tuple(base[control])];remap[vertex]=target;moves.append(dict(vertex=int(vertex),carrier_control=int(control),target_vertex=int(target),maximum_displacement_m=float(distance)))
assert moves;metadata=json.loads((source/'construction.json').read_text());cells={str(i):remap[z['cell_'+str(i)]]for i in range(len(metadata['names']))};normalized=normalize_shared_boundary(p,cells,[0,1.375,.03125]);arrays={'positions':normalized['positions']};rows=[]
for i,name in enumerate(metadata['names']):
 faces=normalized['cells'][str(i)];arrays['cell_'+str(i)]=faces;used,inv=np.unique(faces,return_inverse=True);m=Mesh(arrays['positions'][used],inv.reshape(-1,3),e['frame']);m.save(out/('cell-'+str(i)+'.npz'));proof=audit(m,self_intersections=False,require_closed=True);atomic_json(out/('topology-'+str(i)+'.json'),proof);rows.append(dict(name=name,passed=proof['passed'],failures=proof['failures'],triangles=len(faces)));print(name,proof['passed'],proof['failures'],flush=True)
atomic_npz(out/'cells.npz',**arrays);atomic_json(out/'encoding.json',normalized['evidence']);atomic_json(out/'construction.json',dict(passed=all(q['passed']for q in rows),names=metadata['names'],cells=rows,source_sha256=digest_file(source/'cells.npz'),carrier_controls_unchanged=True,maximum_calculation_star_collapse_m=bound,collapses=moves,full_emitted_audits_pending=True,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False))
