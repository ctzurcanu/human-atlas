"""Two unfinished source-bound Blender diagnostics; no acceptance or packing."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader


def main():
    began=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'unfinished-diagnostic-review-v23';assert not O.exists();O.mkdir()
    binding=json.loads((D/'source-binding-v19/evidence.json').read_bytes());source=ROOT/binding['snapshot'];expected=binding['snapshot_sha256'];assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();wall=Mesh.load(D/'finite-wall.npz');blood=Mesh.load(D/'closed-blood-classifier.npz');frame=wall.frame
    plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());unresolved=json.loads((D/'conservative-envelopes-v21/evidence.json').read_bytes())['unique_unresolved_physical_objects'];objects=['Mandible','Occipital bone','Sphenoid bone','Stylopharyngeus muscle.l','Ascending pharyngeal artery.l']
    pids=[p['id']for p in atlas['parts']if p['id'].startswith('ZA:')and p.get('sourceId')in objects and not p.get('suppressed')]
    pids+=['ZA:Submandibular gland.l','ZA:Parotid gland.l','ZA:Hypoglossal nerve (XII).l','ZA:Superior pharyngeal constrictor.l','TA98TONGUE:male:styloglossus:left','ZA:Stylohyoid muscle.l','ZA:Medial pterygoid muscle.l','ZA:Facial vein.l','ZA:Submandibular duct.l','ZA:Nasopharynx','ZA:Oropharynx','ZA:Laryngopharynx','ZA:Palatine tonsil.l','ZA:Palatopharyngeus muscle.l']
    missing=[pid for pid in pids if pid not in parts];assert not missing,missing
    meshes=[(wall,'HELD v17 actual finite wall · native parent + four branches',[.88,.19,.16,1.],None),(blood,'Diagnostic blood classifier · temporary end caps are NOT tissue',[.2,.35,.65,1.],None)];source_rows=[]
    for pid in dict.fromkeys(pids):
        part=parts[pid];mesh=reader.read_part(atlas,part,frame=frame,catalogue_path=source,catalogue_sha=expected);path=O/(pid.replace(':','_').replace('/','_')+'.npz');mesh.save(path)
        color=([.78,.72,.57,1.]if part.get('system')=='skeletal'else[.92,.7,.13,1.]if part.get('system')=='nervous'else[.2,.43,.82,1.]if part.get('system')=='venous'else[.78,.12,.11,1.]if part.get('system')=='arterial'else[.58,.31,.38,1.]if part.get('system')=='muscular'else[.42,.72,.57,1.])
        if pid in {'ZA:Nasopharynx','ZA:Oropharynx','ZA:Laryngopharynx'}:color=[.4,.72,.9,1.]
        meshes.append((mesh,pid+' · literal native diagnostic context',color,pid));source_rows.append(dict(part_id=pid,sourceId=part.get('sourceId'),source_record=part,geometry_sha256=mesh.geometry_sha,mesh=str(path.relative_to(ROOT)),mesh_sha256=digest_file(path),source_binding=mesh.source))
    whole=[];focused=[];crops=[];low=np.array([-.018,1.607,-.005]);high=np.array([.052,1.689,.065])
    for index,(mesh,name,color,pid)in enumerate(meshes):
        path=(D/('finite-wall.npz'if index==0 else'closed-blood-classifier.npz'))if index<2 else ROOT/next(r['mesh']for r in source_rows if r['part_id']==pid)
        whole.append(dict(mesh=str(path),name=name,color=color))
        tri=mesh.triangles;mask=np.all(tri.max(1)>=low,axis=1)&np.all(tri.min(1)<=high,axis=1)
        if not mask.any():continue
        crop=Mesh(mesh.vertices,mesh.faces[mask],mesh.frame,mesh.attributes,dict(mesh.source,diagnostic_crop_only=True,crop_is_uncapped=True));cp=O/('focused-'+str(index)+'.npz');crop.save(cp);focused.append(dict(mesh=str(cp),name=name+' · uncapped local crop',color=color));crops.append(dict(original_mesh_sha256=digest_file(path),original_geometry_sha256=mesh.geometry_sha,focused_mesh=str(cp.relative_to(ROOT)),focused_mesh_sha256=digest_file(cp),focused_geometry_sha256=crop.geometry_sha,original_triangle_indices=np.flatnonzero(mask).tolist(),no_triangle_or_position_reconstruction=True))
    views=[];sections=[]
    for name,fraction in [('ascending-palatine',.87),('tonsillar',.82)]:
        route=next(r for r in plan['routes']if r['name']==name);z=np.load(ROOT/route['guide']);points=z['points'];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))];k=int(np.argmin(abs(arc-arc[-1]*fraction)));tangent=z['analytic_tangents'][k];normal=np.cross(tangent,[0,1,0])
        if np.linalg.norm(normal)<1e-5:normal=np.cross(tangent,[1,0,0])
        normal/=np.linalg.norm(normal);view=dict(name='UNFINISHED '+name+' actual longitudinal local section · uncapped native mucosal patches',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=points[k].tolist(),clear_outer=True));views.append(view);sections.append(dict(branch=name,guide_sha256=digest_file(ROOT/route['guide']),guide_station=k,point_m=points[k].tolist(),plane_normal=normal.tolist(),tangent_m=tangent.tolist(),section_is_uncapped=True))
    atomic_json(O/'whole-spec.json',dict(frame=frame,resolution=1000,wall_seconds=60,meshes=whole,views=[dict(name='UNFINISHED v17 whole emitted wall/blood and fixed native context · medial oblique',direction=[-1,.08,.22]),dict(name='UNFINISHED v17 parent and four branches · anterior oblique',direction=[.35,.05,1])]))
    atomic_json(O/'focused-spec.json',dict(frame=frame,resolution=1100,wall_seconds=60,meshes=focused,views=[dict(name='UNFINISHED fixed cervical course · actual cropped wall and source context',direction=[-1,.12,.18])]+views))
    evidence=dict(schema='held-facial-diagnostic-review/v22',source_catalogue_sha256=expected,source_snapshot=str(source.relative_to(ROOT)),source_snapshot_sha256=digest_file(source),wall_geometry_sha256=wall.geometry_sha,wall_archive_sha256=digest_file(D/'finite-wall.npz'),blood_geometry_sha256=blood.geometry_sha,blood_archive_sha256=digest_file(D/'closed-blood-classifier.npz'),surface_ancestry_archive_sha256=digest_file(D/'surfaces-and-ancestry.npz'),native_context=source_rows,unresolved_physical_sources=unresolved,actual_native_pharyngeal_patches=['ZA:Nasopharynx','ZA:Oropharynx','ZA:Laryngopharynx'],diagnostic_crops=crops,focus_bounds_m=[low.tolist(),high.tolist()],sections=sections,diagnostic_only=True,unfinished=True,anatomical_acceptance=False,installable=False,owner_packing_performed=False,feasibility_waived=False,fabricated_airspace=False,native_geometry_moved=False,original_user_blender_scene_must_be_preserved=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'source-binding.json',evidence);print(json.dumps(dict(whole_spec=str(O/'whole-spec.json'),focused_spec=str(O/'focused-spec.json'),native_meshes=len(source_rows),focused_meshes=len(focused),wall_geometry_sha256=wall.geometry_sha,elapsed_seconds=evidence['elapsed_seconds'])))


if __name__=='__main__':main()
