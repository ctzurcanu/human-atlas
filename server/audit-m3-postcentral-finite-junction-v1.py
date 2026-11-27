"""Independent local wall/three-port patency checks; no installation."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,topology,winding_contains
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3'
P=B/'m3-postcentral-port-feasibility-v1';C=P/'finite-junction-educational-v1';OUT=C/'intrinsic-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);generation=json.loads((C/'generation-complete.json').read_text());z=np.load(P/'source-port-charts.npz');loaded={};audits={}
 # Recover exact source coordinates after the native millimetre round-trip;
 # only arithmetic roundoff is snapped, never anatomical proximity.
 original=[]
 for role in ['parent','branch']:
  q=np.load(C/(role+'-derived-paired-chart.npz'));original.extend([q['outer_positions'],q['inner_positions']])
 bank=CommonPositionBank(np.concatenate(original),tolerance=2e-12)
 for name,path,genus in [('wall',C/'wall-raw.npz',2),('outer',C/'outer-union.npz',0),('fluid',C/'fluid-union.npz',0)]:
  q=np.load(path);faces=np.asarray(q['faces']);mapped=np.array([bank.insert(p)for p in q['vertices']],int);mesh=Mesh(bank.positions,mapped[faces],FRAME);loaded[name]=(mesh,q['original_id'],q['source_face_id']);mesh.save(OUT/(name+'.npz'));audits[name]=audit(mesh,require_closed=True,expected_components=1,expected_genus=[genus]);atomic_json(OUT/(name+'-audit.json'),audits[name])
  assert audits[name]['passed'],(name,audits[name]['failures'])
 wall,origins,faceids=loaded['wall'];fluid=loaded['fluid'][0];outer=loaded['outer'][0]
 # Provenance separates classifier-only fluid disks and internal root fans
 # from the three permitted annular tissue end rims.
 forbidden=[];port_rows=[];ports=[]
 for cell in generation['classifiers']:
  own=origins==cell['original_id'];junc=np.flatnonzero(own&np.isin(faceids,cell['junction_faces']));forbidden.extend(junc.tolist())
  if cell['role'].endswith('fluid'):forbidden.extend(np.flatnonzero(own&np.isin(faceids,cell['disk_faces'])).tolist())
  else:
   for mouth in cell['mouths']:
    if mouth['kind']!='external-outer-classifier':continue
    ids=np.flatnonzero(own&np.isin(faceids,mouth['newFaceIds']));assert len(ids)>0;rim=Mesh(wall.vertices,wall.faces[ids],FRAME);t=topology(rim);loops=boundary_loops(rim.faces)
    assert len(t['components'])==1 and t['components'][0]['euler']==0 and len(loops)==2
    port_rows.append(dict(role=cell['role'],mouth=mouth['mouth'],faces=len(ids),topology=t,boundary_cycles=len(loops),wall_face_ids=ids.tolist()))
 assert len(port_rows)==3 and not forbidden
 # Source centerline polyline is used as a continuous open-path witness.
 # Test complete segments against all candidate triangles with an exact
 # binary-rational predicate; endpoint winding alone would miss a barrier.
 paths=[z['parent_vertices'][z['parent_station_rings'][33:42]].mean(1),z['branch_vertices'][z['branch_root_station_rings'][176:181]].mean(1)]
 centerproof={name:winding_contains(mesh,np.concatenate(paths),boundary_tolerance_m=1e-11)for name,mesh in [('fluid',fluid),('wall',wall)]}
 assert all(centerproof['fluid']['inside']) and not any(centerproof['fluid']['boundary_ambiguous']);assert not any(centerproof['wall']['inside']) and not any(centerproof['wall']['boundary_ambiguous'])
 exact_tests=0;hits=[];lo=wall.triangles.min(1);hi=wall.triangles.max(1)
 for path_index,path in enumerate(paths):
  for segment,(first,last)in enumerate(zip(path[:-1],path[1:])):
   candidates=np.flatnonzero(np.all(hi>=np.minimum(first,last),axis=1)&np.all(lo<=np.maximum(first,last),axis=1))
   for fid in candidates:
    exact_tests+=1
    if exact_segment_triangle_intersection(first,last,wall.triangles[fid])['intersects']:hits.append(dict(path=path_index,segment=segment,wall_face=int(fid)))
 assert not hits
 # The branch root center and parent station 36 share their actual source
 # point within storage noise; a tiny explicit connecting segment also must
 # stay within blood and avoid wall.
 first,last=paths[1][-1],paths[0][3];connector=winding_contains(fluid,np.array([first,last,(first+last)/2]),boundary_tolerance_m=1e-11);assert all(connector['inside']) and not any(connector['boundary_ambiguous'])
 for fid in np.flatnonzero(np.all(hi>=np.minimum(first,last),axis=1)&np.all(lo<=np.maximum(first,last),axis=1)):
  exact_tests+=1;assert not exact_segment_triangle_intersection(first,last,wall.triangles[fid])['intersects']
 # Literal source outer end controls survive the construction. Native normal
 # attributes of the unchanged continuation are preserved in separate source
 # archives; a joined full-surface normal continuity check remains pending.
 ports=[('parent',33,z['parent_station_rings'][33]),('parent',41,z['parent_station_rings'][41]),('branch',176,z['branch_root_station_rings'][176])]
 port_controls=[];keys={tuple(v)for v in wall.vertices[np.unique(wall.faces)]}
 for role,station,ring in ports:
  q=z[role+'_vertices'][ring];literal=all(tuple(p)in keys for p in q);assert literal;port_controls.append(dict(role=role,station=station,outer_controls=len(q),literal_outer_controls_survive=literal))
 sourceproof=[]
 for role,rows in zip(['parent','branch'],generation['source_rows']):
  mesh=Mesh.load(C/(role+'-source-untouched.npz'));assert mesh.geometry_sha==rows['geometry_sha256'] and mesh.attributes_sha==rows['attributes_sha256'];sourceproof.append(dict(role=role,literal_positions_indices_and_attributes_unchanged=True,geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha))
 proof=dict(schema='local-finite-M3-postcentral-intrinsic/v1',passed=True,installed=False,scope='Only the local three-port finite patch. M3 and postcentral continuations remain source outer sheets; the secondary postcentral crossing is unchanged and unresolved.',wall_faces=len(wall.faces),wall_volume_m3=audits['wall']['topology']['oriented_triangle_integral_m3'],wall_genus=2,fluid_connected_components=1,outer_connected_components=1,three_annular_ports=port_rows,port_outer_controls=port_controls,temporary_fluid_disks_in_wall=0,internal_root_classifier_faces_in_wall=0,continuous_path=dict(polyline_segments=12,exact_triangle_tests=exact_tests,barrier_hits=hits,parent_branch_connector_length_m=float(np.linalg.norm(last-first)),fluid_winding=centerproof['fluid'],wall_winding=centerproof['wall'],connector_fluid_winding=connector),source_preservation=sourceproof,maximum_source_roundoff_snap_m=bank.maximum_snap,new_intersection_point_snaps=bank.new_snaps,remaining=['Actual wall/fluid render and sections.','Full source 192+96-face replacement topology, literal indices, packed normal continuity and shared interface ownership.','Separate local-envelope cortex/secondary-course context checks.','Full vessel finite-wall and second postcentral junction construction.'],elapsed_seconds=time.monotonic()-start)
 atomic_json(OUT/'intrinsic-proof.json',proof)
 spec=dict(frame=FRAME,wall_seconds=40,resolution=900,meshes=[dict(name='Local M3-postcentral finite tissue wall — educational thickness',mesh=str(OUT/'wall.npz'),color=[.72,.15,.12,1.]),dict(name='Connected blood classifier — three open annular mouths',mesh=str(OUT/'fluid.npz'),color=[.2,.4,.85,.24])],views=[dict(name='front',direction=[0,0,1]),dict(name='back',direction=[0,0,-1]),dict(name='side',direction=[1,0,0])]);atomic_json(OUT/'render-spec.json',spec)
 print(json.dumps({k:proof[k]for k in ['passed','wall_faces','wall_volume_m3','wall_genus','fluid_connected_components','temporary_fluid_disks_in_wall','maximum_source_roundoff_snap_m','elapsed_seconds']}))

if __name__=='__main__':main()
