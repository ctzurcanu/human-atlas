"""Measure complete source port curves against the actual installed wall.

Each curve edge is represented as a degenerate triangle only for the distance
kernel, which explicitly retains segment distance. It is never tissue/solid
geometry and is never installed. No point-sample minimum is called a true gap.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import igl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_continuity import boundary_loops
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());pp={p['id']:p for p in a['parts']};reader=AtlasReader();wall=Mesh.load(b/'blood-spaces-v126/installed-boundary.npz');wb=BVH(wall);prior=json.loads((b/'blood-spaces-v126/blood-spaces.json').read_text());rows=[]
for g in prior['gates']:
 pid=g['ring_source'];m=reader.read_part(a,pp[pid],frame=wall.frame,catalogue_path=path);v,iv=np.unique(m.vertices,axis=0,return_inverse=True);f=iv[m.faces];rings=boundary_loops(f);center=np.asarray(g['center_m']);ring=min(rings,key=lambda q:np.linalg.norm(v[q].mean(0)-center));assert np.linalg.norm(v[ring].mean(0)-center)<1e-14
 points=v[ring];N=len(points);edges=np.c_[np.arange(N),np.roll(np.arange(N),-1),np.roll(np.arange(N),-1)];curve=Mesh(points,edges,wall.frame);nearest=BVH(curve).closest_mesh(wb);d2,face_ids,closest=igl.point_mesh_squared_distance(points,wall.vertices,wall.faces);w=igl.winding_number(wall.vertices,wall.faces,points);row=dict(part_id=g['part_id'],ring_source=pid,source_geometry_sha256=m.geometry_sha,whole_wall_geometry_sha256=wall.geometry_sha,source_control_ids=ring.tolist(),source_points_m=points.tolist(),minimum_complete_curve_triangle_distance_m=float(nearest['distance_m']),minimum_witness_curve=nearest['first'].tolist(),minimum_witness_wall=nearest['second'].tolist(),wall_face=int(nearest['second_face_id']),source_edge=int(nearest['first_face_id']),max_vertex_to_wall_distance_m=float(np.sqrt(d2.max())),minimum_vertex_to_wall_distance_m=float(np.sqrt(d2.min())),source_vertex_winding=w.tolist(),vertex_closest_points_m=closest.tolist(),wall_faces=face_ids.tolist(),tested_triangle_pairs=nearest['tested_triangle_pairs'],sample_distances_do_not_prove_whole_curve_contact=True,attachment_not_inferred=True);atomic_json(out/f'root-{len(rows)}.json',row);rows.append(row);print(pid,'minimum whole-curve gap mm',round(row['minimum_complete_curve_triangle_distance_m']*1000,6),'maximum vertex distance mm',round(row['max_vertex_to_wall_distance_m']*1000,3),flush=True)
atomic_json(out/'mouths.json',dict(computed=True,catalogue_sha256=digest_file(path),source_diagnostic_sha256=digest_file(b/'blood-spaces-v126/blood-spaces.json'),roots=rows,elapsed_seconds=time.monotonic()-start,geometry_changed=False,anatomical_acceptance=False))
