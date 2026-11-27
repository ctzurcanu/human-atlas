"""Separately capped full-cell checks or complete centerline patency stage."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,winding_contains
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';C=P/'three-course-finite-educational-v1';A=C/'assembly-v1'

def cell(name):
 start=time.monotonic();out=C/(name+'-intrinsic-v1');out.mkdir(exist_ok=False);mesh=Mesh.load(A/(name+'.npz'));result=audit(mesh,require_closed=True,expected_components=1,expected_genus=[3 if name=='wall'else 0]);result['elapsed_seconds']=time.monotonic()-start;atomic_json(out/'audit.json',result);print(json.dumps(dict(cell=name,passed=result['passed'],failures=result['failures'],triangles=len(mesh.faces),defect_pairs=result['self_intersections'].get('defect_pairs'),elapsed_seconds=result['elapsed_seconds'])));assert result['passed']

def patency():
 start=time.monotonic();out=C/'patency-v1';out.mkdir(exist_ok=False);z=np.load(F/'three-course-source-charts.npz');wall=Mesh.load(A/'wall.npz');fluid=Mesh.load(A/'fluid.npz');paths=[z['M3_vertices'][z['M3_rings']].mean(1),z['branch_vertices'][z['main_rings']].mean(1),z['branch_vertices'][z['twig_rings']].mean(1)];connectors=[np.array([paths[1][-1],paths[0][36]]),np.array([paths[2][-1],paths[1][120]])];points=np.r_[*paths,*connectors];fluidproof=winding_contains(fluid,points,boundary_tolerance_m=1e-11);wallproof=winding_contains(wall,points,boundary_tolerance_m=1e-11);assert all(fluidproof['inside'])and not any(fluidproof['boundary_ambiguous']);assert not any(wallproof['inside'])and not any(wallproof['boundary_ambiguous']);lo=wall.triangles.min(1);hi=wall.triangles.max(1);tests=0;hits=[];segments=0
 for path_index,path in enumerate(paths+connectors):
  for segment,(first,last)in enumerate(zip(path[:-1],path[1:])):
   segments+=1;candidates=np.flatnonzero(np.all(hi>=np.minimum(first,last),axis=1)&np.all(lo<=np.maximum(first,last),axis=1))
   for fid in candidates:
    tests+=1
    if exact_segment_triangle_intersection(first,last,wall.triangles[fid])['intersects']:hits.append(dict(path=path_index,segment=segment,wall_face=int(fid)))
 assert not hits
 # Four true native terminal outer circles survive exactly. No internal
 # prototype annular rims are accepted as substitute terminal openings.
 keys={tuple(v)for v in wall.vertices[np.unique(wall.faces)]};ports=[]
 for role,vkey,rkey,station in [('M3','M3_vertices','M3_rings',0),('M3','M3_vertices','M3_rings',156),('main','branch_vertices','main_rings',0),('twig','branch_vertices','twig_rings',0)]:
  q=z[vkey][z[rkey][station]];assert all(tuple(p)in keys for p in q);ports.append(dict(course=role,station=station,literal_outer_controls=len(q),positive_inner_profile=True))
 assembly=json.loads((A/'assembly-proof.json').read_text());assert assembly['root_fan_or_fluid_disk_faces_in_wall']==0 and len(assembly['native_terminal_annuli'])==4
 result=dict(schema='full-native-three-course-finite-patency/v1',passed=True,installed=False,complete_source_centerline_segments=segments,exact_segment_triangle_tests=tests,wall_barrier_hits=hits,source_station_centers_inside_blood=len(points),fluid_winding=fluidproof,wall_winding=wallproof,connectors_length_m=[float(np.linalg.norm(p[1]-p[0]))for p in connectors],native_terminal_ports=ports,positive_inner_profiles='Authorized educational radial-thickness construction, not measured source anatomy.',scope='Continuous source-centerline clearance witness plus a single connected finite blood classifier, not a clinical flow resistance guarantee.',elapsed_seconds=time.monotonic()-start);atomic_json(out/'patency-proof.json',result);print(json.dumps({k:result[k]for k in ['passed','complete_source_centerline_segments','exact_segment_triangle_tests','source_station_centers_inside_blood','connectors_length_m','elapsed_seconds']}))

if __name__=='__main__':
 if len(sys.argv)!=2 or sys.argv[1]not in ['wall','outer','fluid','patency']:raise SystemExit('Choose one separately capped stage: wall, outer, fluid or patency')
 patency()if sys.argv[1]=='patency'else cell(sys.argv[1])
