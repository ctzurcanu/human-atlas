"""Explain a preserved failed prerequisite; do not retry any native operation."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.audit import winding_contains
from geometry_tools.spatial import BVH
I=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3/m3-postcentral-port-feasibility-v1/temporal-source-incidence-v2';OUT=I/'five-course-fluid-diagnostic-v1'

def main():
 start=time.monotonic();z=np.load(I/'complete-source-charts.npz');center=z['temporal_vertices'][z['component1_rings'][48]].mean(0);child=Mesh.load(OUT/'temporal1-diagnostic-fluid.npz');parent=Mesh.load(OUT/'temporal0-diagnostic-fluid.npz');anchor=np.flatnonzero(np.all(child.vertices==center,axis=1));assert len(anchor)==1;faces=np.flatnonzero(np.any(child.faces==anchor[0],axis=1));assert len(faces)==12;ids=np.unique(child.faces[faces]);points=child.vertices[ids];proof=winding_contains(parent,points,boundary_tolerance_m=1e-10);near,_=BVH(parent).closest_points(points);outside=~np.array(proof['inside'])&~np.array(proof['boundary_ambiguous']);rows=[dict(child_vertex=int(ids[i]),point_m=points[i].tolist(),actual_parent_surface_distance_m=float(near['distance_m'][i]),nearest_parent_point_m=near['points'][i].tolist(),parent_face=int(near['face_ids'][i]))for i in np.flatnonzero(outside)];result=dict(schema='failed-diagnostic-root-prerequisite-explanation/v1',failed_course='temporal1 (49-ring short course)',actual_parent='temporal0 station144',strict_internal_root_fan_inside_vertices=int(np.sum(proof['inside'])),outside_vertices=len(rows),ambiguous_vertices=int(np.sum(proof['boundary_ambiguous'])),outside_witnesses=rows,root_center_inside=bool(proof['inside'][int(np.flatnonzero(ids==anchor[0])[0])]),temporal0_M3_root_passed=True,five_course_union_not_run=True,diagnostic_not_retried=True,reason='The strict shortcut using a root-plane closing fan is not wholly in parent fluid. A true parent-centered root does not establish containment of its entire finite root profile; no topology conclusion may be drawn from an unrun union.',next_requirement='Review this root profile and permitted local sleeve construction separately before any complete network union.',no_fit_or_geometry_change=True,elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'failure-analysis.json',result);print(json.dumps(result))

if __name__=='__main__':main()
