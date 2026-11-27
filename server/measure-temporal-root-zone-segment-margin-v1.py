"""Whole-segment parent-surface clearance and stable blood exit; no fitting."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.distance import triangle_pairs
P=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3/m3-postcentral-port-feasibility-v1/temporal-source-incidence-v2/common-parent-root-zone-v1'

def main():
 start=time.monotonic();a=json.loads((P/'containment-proof.json').read_text());meshes={k:Mesh.load(P/('unchanged-M3-'+k+'-classifier.npz'))for k in ['outer','fluid']};rows=[]
 for index,s in enumerate(a['segments']):
  if not s['family'].startswith('new-'):continue
  endpoints=np.array(s['endpoints_m']);line=np.array([endpoints[0],endpoints[1],endpoints[1]]);row=dict(family=s['family'],temporal_face=s['temporal_face'],temporal_stations=s['temporal_stations'],other_stations=s['other_stations'])
  for kind,m in meshes.items():
   result=triangle_pairs(np.broadcast_to(line,(len(m.faces),3,3)).copy(),m.triangles);j=int(result['distance_m'].argmin());row[kind]=dict(actual_complete_segment_surface_distance_m=float(result['distance_m'][j]),point_on_contact_segment_m=result['first'][j].tolist(),point_on_parent_surface_m=result['second'][j].tolist(),parent_source_face=j)
  rows.append(row)
 blood=[]
 for piece in a['pieces']:
  s=a['segments'][piece['segment']]
  if s['family']=='temporal-centerline'and(piece['fluid_inside']or piece['classification']=='parent-fluid-boundary'):blood.append(s['temporal_stations'][0]+piece['high'])
 report=dict(schema='complete-new-contact-segment-parent-clearance/v1',no_fit_or_geometry_modified=True,segment_margins=rows,temporal_centerline_stable_fluid_exit_station_parameter=max(blood),temporal_centerline_stable_outer_exit_station_parameter=a['temporal_centerline_stable_outer_exit_station_parameter'],method='Actual segment/triangle distance against every complete native M3 classifier face; not endpoint or centroid sampling.',elapsed_seconds=time.monotonic()-start);atomic_json(P/'complete-segment-margin-proof.json',report);print(json.dumps(report))

if __name__=='__main__':main()
