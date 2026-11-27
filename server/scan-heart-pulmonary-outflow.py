"""Locate actual uncapped pulmonary outflow profiles; no vessel motion."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.simple_compartment_section import polygon_contains
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'pulmonary-sections-v47';out.mkdir(exist_ok=False)
m=Mesh.load(b/'whole-wall-v25/whole-wall.npz');g=next(g for g in json.loads((b/'blood-space-v12/compartments.json').read_text())['gates']if g['part_id']=='ZA:Pulmonary trunk');p=np.array(g['center']);n=np.array(g['normal']);n/=np.linalg.norm(n);axis=np.eye(3)[np.argmin(abs(n))];u=np.cross(n,axis);u/=np.linalg.norm(u);v=np.cross(n,u);basis=np.c_[u,v];rows=[]
for offset in np.linspace(0,.025,11):
 center=p+n*offset;scalar=(m.vertices-center)@n;scalar[abs(scalar)<2e-13]=0;cut=clip_scalar_surface(m.vertices,m.faces,scalar,0.);loops=boundary_loops(cut['negative_faces']);contours=[]
 for ids in loops:
  q=cut['positions'][ids];xy=(q-center)@basis;cross=xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1)
  contours.append(dict(points_m=q.tolist(),area_mm2=float(abs(cross.sum())*5e5),centroid_m=q.mean(0).tolist(),axis_inside=polygon_contains(xy,np.zeros(2)),distance_to_axis_mm=float(np.linalg.norm(xy.mean(0))*1000)))
 row=dict(offset_mm=float(offset*1000),plane_center_m=center.tolist(),normal=n.tolist(),contours=contours);rows.append(row);atomic_json(out/(str(round(offset*10000))+'.json'),row);print('depth',round(offset*1000,1),'enclosing',sum(c['axis_inside']for c in contours),'near',[(round(c['area_mm2']),round(c['distance_to_axis_mm'],1))for c in contours if c['distance_to_axis_mm']<25],flush=True)
atomic_json(out/'sections.json',dict(source_geometry_sha256=m.geometry_sha,sections=rows,elapsed_seconds=time.monotonic()-start,installed=False))
