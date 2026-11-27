"""Identify literal native mouth correspondence, without calling proximity a join."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_continuity import boundary_loops
ap=argparse.ArgumentParser();ap.add_argument('--sex',required=True,choices=['male','female']);args=ap.parse_args();sex=args.sex
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/('native-mouths-v219-'+sex);out.mkdir(exist_ok=False)
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json';a=json.loads(cat.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres'
pairs=[('artery','ZA:Coeliac trunk'if sex=='male'else'LOCAL:female:common_hepatic_artery:1722'),('portal','TA98HEPATIC:'+sex+':portal-distal'),('hepatic_vein','ZA:Inferior vena cava (abdominal part)'if sex=='male'else'LOCAL:female:vena_cava:1797'),('biliary','TA98PART:FJ3079:A05.8.01.061'if sex=='male'else'TA98HEPATIC:female:common-hepatic-parent')]
rows=[]
for role,pid in pairs:
 first=reader.read_part(a,parts['TA98HEPATIC:'+sex+':'+role],frame=frame,catalogue_path=cat);second=reader.read_part(a,parts[pid],frame=frame,catalogue_path=cat)
 v1=np.unique(first.vertices[np.unique(first.faces)],axis=0);v2,inv=np.unique(second.vertices[np.unique(second.faces)],axis=0,return_inverse=True);used=np.unique(second.faces);mapping=np.full(len(second.vertices),-1);mapping[used]=inv;sf=mapping[second.faces];distance,which=cKDTree(v1).query(v2)
 try:loops=boundary_loops(sf)
 except ValueError:loops=[]
 summary=[]
 for loop in loops:
  summary.append(dict(points=len(loop),center_m=v2[loop].mean(0).tolist(),coordinates=v2[loop].tolist(),max_nearest_root_vertex_m=float(distance[loop].max()),all_vertices_literal=bool(np.all(distance[loop]==0))))
 row=dict(role=role,parent=pid,literal_shared_vertices=int(np.count_nonzero(distance==0)),parent_boundary_loops=summary,minimum_vertex_to_vertex_m=float(distance.min()),nearest_actual_surfaces=BVH(first).closest_mesh(BVH(second)))
 rows.append(row);print(json.dumps({k:v for k,v in row.items()if k not in ['nearest_actual_surfaces','parent_boundary_loops']}),flush=True)
atomic_json(out/'mouths.json',dict(sex=sex,interfaces=rows,not_anatomical_acceptance=True))
