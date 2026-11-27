"""Actual uncapped wall sections at the eight native vessel-root planes."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.simple_compartment_section import polygon_contains,polygon_clearance
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='root-sections-v33');args=ap.parse_args()
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';out=base/args.stage;out.mkdir(exist_ok=False);m=Mesh.load(base/'whole-wall-v25/whole-wall.npz');gates=json.loads((base/'blood-space-v12/compartments.json').read_text())['gates'];rows=[]
for g in gates:
 p=np.array(g['center']);n=np.array(g['normal']);n/=np.linalg.norm(n);axis=np.eye(3)[np.argmin(abs(n))];u=np.cross(n,axis);u/=np.linalg.norm(u);v=np.cross(n,u);basis=np.c_[u,v];scalar=(m.vertices-p)@n;scalar[abs(scalar)<=2e-13]=0;cut=clip_scalar_surface(m.vertices,m.faces,scalar,0.);rings=boundary_loops(cut['negative_faces']);sections=[]
 for ring in rings:
  world=cut['positions'][ring];xy=(world-p)@basis;cross=xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1);area=float(abs(cross.sum()/2));sections.append(dict(points_m=world.tolist(),area_m2=area,contains_native_root_centre=polygon_contains(xy,np.zeros(2)),clearance_m=polygon_clearance(xy,np.zeros(2)),centroid_distance_m=float(np.linalg.norm(world.mean(0)-p)),maximum_root_radius_m=float(np.linalg.norm(xy,axis=1).max())))
 enclosing=sorted([s for s in sections if s['contains_native_root_centre']],key=lambda q:q['area_m2']);row=dict(part_id=g['part_id'],plane_point_m=p.tolist(),plane_normal=n.tolist(),enclosing_contours=len(enclosing),candidate_inner_lumen=enclosing[0]if len(enclosing)>=2 else None,contours=sections);rows.append(row);atomic_json(out/(g['part_id'].replace(':','_').replace(' ','_')+'.json'),row);print(g['part_id'],'enclosing',len(enclosing),'clearance mm',None if len(enclosing)<2 else round(enclosing[0]['clearance_m']*1000,3),flush=True)
atomic_json(out/'sections.json',dict(source_geometry_sha256=m.geometry_sha,roots=rows,caps_created=0,source_changed=False,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False))
