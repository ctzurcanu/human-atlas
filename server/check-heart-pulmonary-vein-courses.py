"""Measure actual topological native PV rings, rather than PCA-axis guesses."""
from pathlib import Path
import sys,time,json
from collections import deque
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.tube_continuity import boundary_loops
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'pv-course-v57';out.mkdir(exist_ok=False);a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());m=Mesh.load(b/'whole-wall-v25/whole-wall.npz');o=np.array([0.,1.375,.03125]);scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces.astype('uint32')));gates=json.loads((b/'blood-space-v12/compartments.json').read_text())['gates'];reader=AtlasReader();rows=[]
for gate in gates:
 if 'pulmonary vein'not in gate['part_id']:continue
 pid=gate['part_id'];part=next(p for p in a['parts']if p['id']==pid);vessel=reader.read_part(a,part,frame=m.frame);v,inv=np.unique(vessel.vertices,axis=0,return_inverse=True);faces=inv[vessel.faces];loops=boundary_loops(faces);root=np.array(gate['center']);adj=[set()for _ in v]
 for f in faces:
  for k in range(3):adj[f[k]].add(int(f[(k+1)%3]));adj[f[(k+1)%3]].add(int(f[k]))
 for j,loop in enumerate(loops):
  if np.linalg.norm(v[loop].mean(0)-root)>.0001:continue
  distance=np.full(len(v),-1);distance[loop]=0;queue=deque(map(int,loop))
  while queue:
   k=queue.popleft()
   for q in adj[k]:
    if distance[q]<0:distance[q]=distance[k]+1;queue.append(q)
  levels=[np.flatnonzero(distance==d)for d in range(distance.max()+1)];assert all(len(q)==len(loop)for q in levels),'Unexpected source ring topology';centers=np.array([v[q].mean(0)for q in levels]);profiles=[v[q]for q in levels];samples=[]
  for k in range(len(centers)-1):
   count=max(2,int(np.ceil(np.linalg.norm(centers[k+1]-centers[k])/.0004))+1)
   for t in np.linspace(0,1,count):samples.append(centers[k]*(1-t)+centers[k+1]*t)
  points=np.array(samples);d=scene.compute_signed_distance(o3d.core.Tensor((points-o).astype('float32')),nthreads=1,nsamples=3).numpy();profile_points=np.concatenate([c+.8*(q-c)for c,q in zip(centers,profiles)]);radial=scene.compute_signed_distance(o3d.core.Tensor((profile_points-o).astype('float32')),nthreads=1,nsamples=3).numpy();row=dict(part_id=pid,branch=j,centers_m=centers.tolist(),source_ring_size=len(loop),samples=len(points),axis_minimum_distance_m=float(d.min()),axis_inside_tissue_count=int(np.sum(d<-.00001)),profile_inset_fraction=.8,profile_minimum_distance_m=float(radial.min()),profile_inside_tissue_count=int(np.sum(radial<-.00001)));rows.append(row);print(pid,j,'axis mm',round(float(d.min()*1000),3),'axis intrusions',row['axis_inside_tissue_count'],'80% profile intrusions',row['profile_inside_tissue_count'],flush=True);atomic_npz(out/(pid.replace(':','_').replace(' ','_')+'-'+str(j)+'.npz'),centers=centers,profiles=np.array(profiles),source_frame=np.array(m.frame))
atomic_json(out/'evidence.json',dict(source_geometry_sha256=m.geometry_sha,branches=rows,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False))
