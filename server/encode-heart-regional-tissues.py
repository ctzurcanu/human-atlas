"""Normalize and inventory the actual common Float32 cardiac boundary cells."""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
ap=argparse.ArgumentParser();ap.add_argument('--source',default='regional-tissues-v87');ap.add_argument('--stage',default='regional-encoded-v88');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/args.source;out=b/args.stage;out.mkdir(exist_ok=False);z=np.load(source/'cells.npz');metadata=json.loads((source/'construction.json').read_text());names=metadata['names'];o=np.array([0.,1.375,.03125]);raw={str(i):z['cell_'+str(i)]for i in range(len(names))if len(z['cell_'+str(i)])};r=normalize_shared_boundary(z['positions'],raw,o);atomic_json(out/'encoding.json',r['evidence']);arrays={'positions':r['positions']};rows=[]
for i,name in enumerate(names):
 if str(i)not in r['cells']:continue
 faces=r['cells'][str(i)];arrays['cell_'+str(i)]=faces;used,inv=np.unique(faces,return_inverse=True);m=Mesh(r['positions'][used],inv.reshape(-1,3),'model:local-ta98:assembled-metres');m.save(out/('cell-'+str(i)+'.npz'));proof=audit(m,self_intersections=False,require_closed=True);atomic_json(out/('topology-'+str(i)+'.json'),proof);rows.append(dict(index=i,name=name,triangles=len(faces),geometry_sha256=m.geometry_sha,passed=proof['passed'],failures=proof['failures'],volume_mL=proof['topology']['oriented_triangle_integral_m3']*1e6));print(i,name,rows[-1]['passed'],rows[-1]['failures'],flush=True);atomic_json(out/'progress.json',dict(cells=rows,elapsed_seconds=time.monotonic()-start),overwrite=True)
atomic_npz(out/'cells.npz',**arrays);atomic_json(out/'construction.json',dict(passed=all(row['passed']for row in rows),source_sha256=digest_file(source/'cells.npz'),names=names,cells=rows,positions=len(r['positions']),triangles=sum(row['triangles']for row in rows),elapsed_seconds=time.monotonic()-start,full_self_checks_pending=True,installed=False,anatomical_acceptance=False));print('Encoded cardiac topology stage complete',round(time.monotonic()-start,2),'s',flush=True)
