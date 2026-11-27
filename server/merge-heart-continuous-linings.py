"""Keep anatomical endocardium/epicardium continuous, not calculation charts."""
from pathlib import Path
import sys,time,json,argparse
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/args.source;out=b/args.stage;out.mkdir(exist_ok=False);z=np.load(source/'cells.npz');metadata=json.loads((source/'construction.json').read_text());names=metadata['names'][:8]+['endocardium','epicardium'];cells={str(i):z['cell_'+str(i)]for i in range(8)};cells['8']=np.concatenate([z['cell_'+str(i)]for i in range(8,16)if 'cell_'+str(i)in z]);cells['9']=np.concatenate([z['cell_'+str(i)]for i in range(16,24)if 'cell_'+str(i)in z]);r=normalize_shared_boundary(z['positions'],cells,[0,1.375,.03125]);arrays={'positions':r['positions']};rows=[]
for i,name in enumerate(names):
 f=r['cells'][str(i)];arrays['cell_'+str(i)]=f;ids,inv=np.unique(f,return_inverse=True);m=Mesh(r['positions'][ids],inv.reshape(-1,3),'model:local-ta98:assembled-metres');m.save(out/('cell-'+str(i)+'.npz'));proof=audit(m,self_intersections=False,require_closed=True);atomic_json(out/('topology-'+str(i)+'.json'),proof);row=dict(index=i,name=name,triangles=len(f),geometry_sha256=m.geometry_sha,passed=proof['passed'],failures=proof['failures'],volume_mL=proof['topology']['oriented_triangle_integral_m3']*1e6);rows.append(row);print(row,flush=True)
atomic_npz(out/'cells.npz',**arrays);atomic_json(out/'encoding.json',r['evidence']);atomic_json(out/'construction.json',dict(passed=all(row['passed']for row in rows),source_sha256=digest_file(source/'cells.npz'),names=names,cells=rows,physical_tissues=len(names),new_images=0,new_uvs=0,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False))
