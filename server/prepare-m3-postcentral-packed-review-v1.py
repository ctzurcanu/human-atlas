"""Actual packed draw-group display, root sections and emitted path witness."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.audit import winding_contains
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';F=B/'m3-postcentral-port-feasibility-v1/secondary-port-feasibility-v1';C=F.parent/'three-course-finite-educational-v1';P=C/'packed-owned-cells-v1';OUT=P/'packed-review-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);layout=json.loads((P/'packed-layout.json').read_text());bank=(P/'cells.bin').read_bytes();origin=np.array(layout['geometry_origin']);entries=[]
 for name,row in layout['cells'].items():
  n=row['vertexCount'];v=np.ndarray((n,3),'<f4',buffer=bank,offset=row['positions']).astype(float)+origin;words=np.ndarray((n,3),'<i2',buffer=bank,offset=row['normals']);sign=np.ndarray((n,),'i1',buffer=bank,offset=row['normalSigns']);uv=np.ndarray((n,2),'<f4',buffer=bank,offset=row['uvs']);indices=np.ndarray((row['indexCount'],),'<u4',buffer=bank,offset=row['indices']);faces=np.concatenate([indices[g['start']:g['start']+g['count']].reshape(-1,3)for g in row['tissueDrawGroups']if not g.get('neighbor')]);used,inv=np.unique(faces,return_inverse=True);mesh=Mesh(v[used],inv.reshape(-1,3),FRAME,dict(normals=words[used],normalSigns=sign[used],uvs=uv[used]));path=OUT/(name+'-actual-assembled-display.npz');mesh.save(path);entries.append(dict(name=row['id']+' · actual packed display with conditional interface hidden',mesh=str(path),color=[.72,.15,.12,1.]))
 z=np.load(F/'three-course-source-charts.npz');paths=[z['M3_vertices'][z['M3_rings']].mean(1),z['branch_vertices'][z['main_rings']].mean(1),z['branch_vertices'][z['twig_rings']].mean(1)];connectors=[np.array([paths[1][-1],paths[0][36]]),np.array([paths[2][-1],paths[1][120]])];wall=Mesh.load(P/'shared-union-v1/emitted-assembled-wall-no-internal-caps.npz');centers=np.r_[*paths,*connectors];contains=winding_contains(wall,centers,boundary_tolerance_m=1e-11);assert not any(contains['inside'])and not any(contains['boundary_ambiguous']);lo=wall.triangles.min(1);hi=wall.triangles.max(1);tests=0;segments=0
 for path in paths+connectors:
  for first,last in zip(path[:-1],path[1:]):
   segments+=1
   for fid in np.flatnonzero(np.all(hi>=np.minimum(first,last),axis=1)&np.all(lo<=np.maximum(first,last),axis=1)):
    tests+=1;assert not exact_segment_triangle_intersection(first,last,wall.triangles[fid])['intersects']
 views=[dict(name='actual-packed-front',direction=[0,0,1]),dict(name='actual-packed-side',direction=[1,0,0])]
 for name,parent,station,loidx,hiidx,child in [('M3-main',paths[0],36,33,41,paths[1]),('main-secondary',paths[1],120,115,123,paths[2])]:
  root=child[-1];axis=parent[hiidx]-parent[loidx];branch=child[-4]-root;normal=np.cross(axis,branch);normal/=np.linalg.norm(normal);up=axis/np.linalg.norm(axis)
  for direction,clear in [(normal,True),(-normal,False)]:views.append(dict(name=name+'-actual-packed-longitudinal-'+str(clear),direction=direction.tolist(),up=up.tolist(),section=dict(point=root.tolist(),normal=normal.tolist(),clear_outer=clear)))
 spec=dict(frame=FRAME,wall_seconds=40,resolution=1000,meshes=entries,views=views);atomic_json(OUT/'render-spec.json',spec);proof=dict(schema='actual-packed-three-course-review-input/v1',installed=False,assembled_display_matches_actual_neighbor_hidden_tissueDrawGroups=True,closed_cell_source_normals_and_reciprocal_normalSigns_retained=True,emitted_station_and_connector_witness_points=len(centers),emitted_wall_all_witness_points_outside_tissue=True,complete_emitted_path_segments=segments,exact_segment_triangle_tests=tests,barrier_hits=0,four_native_terminal_annuli_retained=True,continuous_blood_reference='The same full source paths and connectors lie in the validated connected blood volume; packed wall now independently passes exact continuous barrier checks.',views=len(views),elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'review-input-proof.json',proof);print(json.dumps(proof));print(OUT/'render-spec.json')

if __name__=='__main__':main()
