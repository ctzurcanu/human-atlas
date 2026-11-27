"""Measure the native whole cardiac wall before any reconstruction."""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.tube_continuity import boundary_loops
base=ROOT/'.local-models/geometry-workbench/heart';out=base/'wall-seams-v2';out.mkdir(exist_ok=False);source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(source.read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();names=['ZA:Right atrium','ZA:Right ventricle','ZA:Left atrium','ZA:Left ventricle'];meshes=[r.read_part(a,parts[pid],frame='model:local-ta98:assembled-metres')for pid in names];vertices=np.concatenate([m.vertices for m in meshes]);offset=np.cumsum([0]+[len(m.vertices)for m in meshes]);faces=np.concatenate([m.faces+o for m,o in zip(meshes,offset)]);v,inv=np.unique(vertices,axis=0,return_inverse=True);f=inv[faces];m=Mesh(v,f,meshes[0].frame);edges,count=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0,return_counts=True);boundary=edges[count==1];bids=np.unique(boundary);print('Original whole wall:',len(v),'vertices',len(f),'faces',len(boundary),'open edges',flush=True)
rows=[];boundary_points=[]
for name,mesh in zip(names,meshes):
 vv,ii=np.unique(mesh.vertices,axis=0,return_inverse=True);ff=ii[mesh.faces];ee,cc=np.unique(np.sort(np.concatenate([ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]]),axis=1),axis=0,return_counts=True);points=vv[np.unique(ee[cc==1])];boundary_points.append(points)
for i in range(4):
 for j in range(i+1,4):
  d,ix=cKDTree(boundary_points[j]).query(boundary_points[i]);take=d<.0005;rows.append(dict(first=names[i],second=names[j],first_boundary_controls=len(d),within_half_mm=int(take.sum()),within_one_micron=int((d<1e-6).sum()),minimum_distance_m=float(d.min()),close_quantiles_m=np.quantile(d[take],[0,.5,1]).tolist()if take.any()else[]));print(rows[-1],flush=True)
atomic_npz(out/'wall.npz',vertices=v,faces=f,frame=np.array(m.frame),face_owner=np.repeat(np.arange(4),[len(m.faces)for m in meshes]))
try:loops=boundary_loops(f);loop_rows=[dict(indices=q.tolist(),centre=v[q].mean(0).tolist(),vertices=len(q),bounds=np.stack([v[q].min(0),v[q].max(0)]).tolist())for q in loops]
except ValueError as exc:loop_rows=[];print('Boundary graph:',exc,flush=True)
atomic_json(out/'seams.json',dict(source_catalogue_sha256=digest_file(source),part_ids=names,unmodified_whole_wall_topology=topology(m),boundary_loops=loop_rows,measured_candidate_interfaces=rows,geometry_modified=False,anatomical_acceptance=False))
