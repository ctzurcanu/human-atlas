"""Actual myocardial exclusion along each native root axis and source ring."""
from pathlib import Path
import sys,time,json
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh,AtlasReader
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'root-axis-v54';out.mkdir(exist_ok=False);m=Mesh.load(b/'whole-wall-v25/whole-wall.npz');o=np.array([0.,1.375,.03125]);scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces.astype('uint32')));gates=json.loads((b/'blood-space-v12/compartments.json').read_text())['gates'];rows=[]
for g in gates:
 p=np.array(g['center']);n=np.array(g['normal']);n/=np.linalg.norm(n);offset=np.arange(-.004,.020001,.0005);points=p+offset[:,None]*n;hit=scene.compute_signed_distance(o3d.core.Tensor((points-o).astype('float32')),nthreads=1,nsamples=3).numpy();row=dict(part_id=g['part_id'],center_m=p.tolist(),normal=n.tolist(),offsets_m=offset.tolist(),signed_tissue_distances_m=hit.tolist(),axis_tissue_intervals_mm=[float(v*1000)for v in offset[hit<-.00001]],minimum_tissue_clearance_mm=float(hit.min()*1000));rows.append(row);print(g['part_id'],'tissue intervals mm',row['axis_tissue_intervals_mm'],'min clearance mm',round(row['minimum_tissue_clearance_mm'],3),flush=True)
atomic_json(out/'evidence.json',dict(source_geometry_sha256=m.geometry_sha,roots=rows,scope='Myocardial axis exclusion, not a complete root-wall interface or valve certificate',elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False))
