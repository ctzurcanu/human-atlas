"""Normalize one shared Float32 bank and measure physical lining locations."""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='tissues-encoded-v55');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);source=b/'tissues-v53';z=np.load(source/'cells.npz');metadata=json.loads((source/'construction.json').read_text());names=metadata['cells'];o=np.array([0.,1.375,.03125]);cells={str(i):z['cell_'+str(i)]for i in range(len(names))};r=normalize_shared_boundary(z['positions'],cells,o);atomic_json(out/'encoding.json',r['evidence']);p=r['positions'];arrays={'positions':p};counts={}
for i,name in enumerate(names):
 f=r['cells'][str(i)];arrays['cell_'+str(i)]=f;used,inv=np.unique(f,return_inverse=True);Mesh(p[used],inv.reshape(-1,3),'model:local-ta98:assembled-metres').save(out/('cell-'+str(i)+'.npz'));counts[name]=len(f)
atomic_npz(out/'cells.npz',**arrays);print('Shared encoded cells saved',len(p),'positions',round(time.monotonic()-start,2),'s',flush=True)
m=Mesh.load(b/'whole-wall-v25/whole-wall.npz');inner=np.load(b/'wall-regions-v27/charts.npz')['endocardial_faces'];proof=[]
for kind,condition in [('endocardium',inner),('epicardium',~inner)]:
 scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces[condition].astype('uint32')))
 for i,name in enumerate(names):
  if not name.startswith(kind+':'):continue
  f=r['cells'][str(i)];used=np.unique(f);points=np.r_[p[used],p[f].mean(1)];d=scene.compute_distance(o3d.core.Tensor((points-o).astype('float32')),nthreads=1).numpy();row=dict(name=name,samples=len(points),distance_quantiles_m=np.quantile(d,[0,.5,.9,.99,1]).tolist(),samples_over_80um=int(np.sum(d>.00008)),maximum_m=float(d.max()));proof.append(row);print(name,'max distance mm',round(float(d.max()*1000),3),'over80um',row['samples_over_80um'],flush=True)
report=dict(constructed=True,source_sha256=digest_file(source/'cells.npz'),source_geometry_sha256=m.geometry_sha,names=names,positions=len(p),cell_triangles=counts,coating_location_checks=proof,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print('Encoded lining measurement finished',round(time.monotonic()-start,2),'s',flush=True)
