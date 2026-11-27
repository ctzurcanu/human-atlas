"""Diagnostic root/parent-axis sections with actual source-derived shading."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.algorithms.finite_tissue_boolean import barycentric
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';C=P/'finite-junction-educational-v1';OUT=C/'longitudinal-section-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(P/'source-port-charts.npz');q=np.load(C/'wall-raw.npz');g=json.loads((C/'generation-complete.json').read_text());wall=Mesh.load(C/'intrinsic-v1/wall.npz');origins=q['original_id'];faceids=q['source_face_id'];rows={r['original_id']:r for r in g['classifiers']};sources={r:Mesh.load(C/(r+'-source-untouched.npz'))for r in ['parent','branch']};inputs={};input_normals={}
 for path in sorted(C.glob('*-make-*-input.npz')):
  role=path.name.split('-make-',1)[1].split('-input',1)[0];r=np.load(path);positions=r['positions_mm']/1000.+np.array(g['origin_m']);f=r['faces'];inputs[role]=(positions,f);normals=np.zeros_like(positions);count=next(x['source_faces']for x in rows.values()if x['role']==role);fs=f[:count];t=positions[fs];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
  for corner in range(3):np.add.at(normals,fs[:,corner],n)
  length=np.linalg.norm(normals,axis=1);normals=np.divide(normals,length[:,None],out=np.zeros_like(normals),where=length[:,None]>0);input_normals[role]=normals
 result=[];kinds={'native_outer_interpolation':0,'derived_inner_interpolation':0,'annular_port_or_numerical_extension_face_normal':0}
 for i,(oid,fid)in enumerate(zip(origins,faceids)):
  row=rows[int(oid)];role=row['role'];component=role.split('-')[0];points=wall.triangles[i]
  if int(fid)<row['source_faces']and role.endswith('outer'):
   sid=int(z[component+'_local_patch_face_ids'][int(fid)]);source=sources[component];f=source.faces[sid];weights=barycentric(points,source.vertices[f]);values=source.attributes['normals'][f].astype(float)
   if 'normalSigns'in source.attributes:values*=source.attributes['normalSigns'][f,None]
   n=weights@values;kinds['native_outer_interpolation']+=1
  elif int(fid)<row['source_faces']and role.endswith('fluid'):
   vp,fs=inputs[role];f=fs[int(fid)];weights=barycentric(points,vp[f]);n=-(weights@input_normals[role][f]);kinds['derived_inner_interpolation']+=1
  else:
   t=points;n=np.broadcast_to(np.cross(t[1]-t[0],t[2]-t[0]),(3,3)).copy();kinds['annular_port_or_numerical_extension_face_normal']+=1
  length=np.linalg.norm(n,axis=1);assert np.isfinite(n).all()and np.all(length>0);result.append(n/length[:,None])
 normals=np.array(result);atomic_npz(OUT/'wall-with-source-and-derived-corner-normals.npz',vertices=wall.vertices,faces=wall.faces,corner_normals=normals,frame=np.array(FRAME));root=z['branch_vertices'][z['branch_root_station_rings'][180]].mean(0);parent_axis=z['parent_vertices'][z['parent_station_rings'][41]].mean(0)-z['parent_vertices'][z['parent_station_rings'][33]].mean(0);branch_axis=z['branch_vertices'][z['branch_root_station_rings'][176]].mean(0)-root;normal=np.cross(parent_axis,branch_axis);normal/=np.linalg.norm(normal);up=parent_axis/np.linalg.norm(parent_axis)
 views=[dict(name='longitudinal-root-cut-near',direction=normal.tolist(),up=up.tolist(),section=dict(point=root.tolist(),normal=normal.tolist(),clear_outer=True)),dict(name='longitudinal-root-cut-far',direction=(-normal).tolist(),up=up.tolist(),section=dict(point=root.tolist(),normal=normal.tolist(),clear_outer=False))]
 spec=dict(frame=FRAME,wall_seconds=40,resolution=1000,meshes=[dict(name='Local finite tissue wall with source outer / derived inner corner normals',mesh=str(OUT/'wall-with-source-and-derived-corner-normals.npz'),color=[.72,.15,.12,1.])],views=views);atomic_json(OUT/'section-render-spec.json',spec);atomic_json(OUT/'section-proof.json',dict(scope='Actual longitudinal diagnostic cut through literal branch root and M3 local axis; display cut is deliberately uncapped.',source_outer_normals='Interpolated original native packed normal words/signs.',inner_normals='Area-weighted normals of the actually derived inner station mesh; not measured native inner normals.',rim_normals='Actual annular/classifier-extension output face normals.',normal_face_counts=kinds,point_m=root.tolist(),plane_normal=normal.tolist(),parent_axis=up.tolist(),elapsed_seconds=time.monotonic()-start));print(OUT/'section-render-spec.json')

if __name__=='__main__':main()
