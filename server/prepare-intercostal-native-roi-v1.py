"""Literal face-AABB review crops; no source cuts or camera-product changes."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/posterior-intercostal-native-v1/t10-left-spinal-contract-v1'

def crop_spec(spec,points,padding,outdir,*,label):
    started=time.monotonic();outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=False);points=np.asarray(points,float);padding=np.asarray(padding,float);low=points.min(0)-padding;high=points.max(0)+padding;rows=[];meshes=[];owned=set();combined=[]
    for i,item in enumerate(spec['meshes']):
        src=Mesh.load(item['mesh']);tri=src.triangles;keep=np.all(tri.max(1)>=low,axis=1)&np.all(tri.min(1)<=high,axis=1);face_ids=np.flatnonzero(keep)
        if not len(face_ids):continue
        faces=src.faces[face_ids];ids=np.unique(faces);mapping=np.full(len(src.vertices),-1,int);mapping[ids]=np.arange(len(ids));attrs={k:v[ids]for k,v in src.attributes.items()};part=src.source['source_id'];native_ids=np.asarray(src.source['literal_source_face_ids'])[face_ids]
        for fid in native_ids:
            key=(part,int(fid));assert key not in owned,'Duplicate native chart/context face';owned.add(key)
        crop=Mesh(src.vertices[ids],mapping[faces],src.frame,attrs,dict(source_id=part,input_mesh=str(Path(item['mesh']).resolve()),input_geometry_sha256=src.geometry_sha,input_attributes_sha256=src.attributes_sha,source_geometry_sha256=src.source['source_geometry_sha256'],source_attributes_sha256=src.source['source_attributes_sha256'],literal_source_face_ids=native_ids.tolist(),roi_low_m=low.tolist(),roi_high_m=high.tolist(),literal_intersecting_triangle_aabbs_only=True,no_vertex_motion_or_added_faces=True))
        assert np.array_equal(crop.triangles,src.triangles[face_ids]);assert all(crop.attributes[k].dtype==src.attributes[k].dtype and crop.attributes[k].tobytes()==src.attributes[k][ids].tobytes()for k in attrs)
        target=outdir/f'{i:02}-{Path(item["mesh"]).name}';crop.save(target);meshes.append({**item,'mesh':str(target),'name':item.get('name',part)+' — literal ROI source faces'});combined.append(crop.bounds)
        rows.append(dict(source_id=part,input_mesh=str(item['mesh']),input_geometry_sha256=src.geometry_sha,input_attributes_sha256=src.attributes_sha,input_face_count=len(src.faces),retained_input_face_ids=face_ids.tolist(),retained_native_source_face_ids=native_ids.tolist(),retained_source_vertex_ids=ids.tolist(),output_mesh=str(target),output_geometry_sha256=crop.geometry_sha,output_attributes_sha256=crop.attributes_sha,attribute_dtypes={k:str(v.dtype)for k,v in attrs.items()},normal_and_uv_words_bit_exact=True,no_vertex_motion=True,no_added_faces=True,actual_output_bounds_m=crop.bounds.tolist()))
    midpoint=points.mean(0);render=dict(frame=spec['frame'],resolution=900,wall_seconds=40,meshes=meshes,views=[dict(name=label+'-oblique',direction=[1,.15,-.35]),dict(name=label+'-posterior',direction=[0,0,-1]),dict(name=label+'-takeoff-axial',direction=[0,-1,0],section=dict(normal=[0,1,0],point=points[0].tolist(),clear_outer=True)),dict(name=label+'-recipient-axial',direction=[0,-1,0],section=dict(normal=[0,1,0],point=points[-1].tolist(),clear_outer=True))],roi_only=True,no_new_anatomical_geometry=True,whole_view_request_retained=True)
    atomic_json(outdir/'render-spec.json',render);atomic_json(outdir/'binding.json',dict(schema='literal-native-mesh-roi-review/v1',label=label,roi_landmarks_m=points.tolist(),roi_padding_m=padding.tolist(),requested_bounds_m=[low.tolist(),high.tolist()],actual_retained_face_bounds_m=[np.min(np.array(combined)[:,0],axis=0).tolist(),np.max(np.array(combined)[:,1],axis=0).tolist()],selection='Keep every original triangle whose AABB intersects the requested ROI; compact only unused controls. Large retained original triangles may extend beyond the ROI.',all_fields_normals_uv_words_preserved=True,no_duplicated_native_source_face=True,no_vertex_motion=True,no_mesh_capping=True,source_records=rows,elapsed_seconds=time.monotonic()-started))
    return dict(render_spec=str(outdir/'render-spec.json'),meshes=len(meshes),faces=len(owned),requested_bounds=[low.tolist(),high.tolist()],elapsed_seconds=time.monotonic()-started)

if __name__=='__main__':
    c=json.loads((B/'contract.json').read_bytes());s=json.loads((B/'render-spec.json').read_bytes());print(json.dumps(crop_spec(s,[c['takeoff_chart']['point_m'],c['recipient_chart']['point_m']],[.012,.014,.016],B/'native-roi-v1',label='T10-left-native-chart-ROI')))
