"""Source replay, shared canonical coordinates, literal outer attributes and ports."""
from pathlib import Path
import json,sys,time
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank,barycentric
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';C=P/'three-course-finite-educational-v1';OUT=C/'assembly-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);g=json.loads((C/'generation-complete.json').read_text());z=np.load(F/'three-course-source-charts.npz');source=[Mesh.load(P/'finite-junction-educational-v1/parent-source-untouched.npz'),Mesh.load(P/'finite-junction-educational-v1/branch-source-untouched.npz')];prefix=np.r_[source[0].vertices,source[1].vertices];offset=len(source[0].vertices);rows={x['original_id']:x for x in g['classifiers']};charts={r:np.load(C/(r+'-derived-paired-chart.npz'))for r in ['M3','main','twig']};bank=CommonPositionBank(np.r_[prefix,*[q['inner_positions']for q in charts.values()]],tolerance=2e-12);meshes={};archives={};raw={}
 for role,path in [('wall','wall-raw.npz'),('outer','outer-union.npz'),('fluid','fluid-union.npz')]:
  q=np.load(C/path);mapped=np.array([bank.insert(v)for v in q['vertices']],int);raw[role]=(mapped[q['faces']],q['original_id'],q['source_face_id']);archives[role]=q
 positions=bank.positions
 # Preserve all old index corners outside explicit domains. A source face
 # must survive literally once, not just have similar area or a nearby fit.
 allowed={'M3':set(z['M3_root_patch_face_ids'].tolist()),'main':set(np.r_[z['main_root_to_M3_patch_face_ids'],z['main_secondary_patch_face_ids']].tolist()),'twig':set(z['twig_root_patch_face_ids'].tolist())};proof=[];fullouter=[];seen=defaultdict(list);localouter=[]
 for tri,oid,fid in zip(*raw['outer']):
  row=rows[int(oid)]
  if int(fid)>=row['source_faces']:continue
  course=row['course'];sid=row['source_face_ids'][int(fid)];index=0 if course=='M3'else 1;old=source[index].faces[sid]+(0 if index==0 else offset);points=positions[tri];sq=positions[old];matches=np.all(points[:,None,:]==sq[None,:,:],axis=2);literal=np.all(matches.sum(1)==1)and len(set(np.argmax(matches,axis=1)))==3
  if literal:
   tri=old[np.argmax(matches,axis=1)]
  seen[(course,sid)].append((tri,bool(literal)));fullouter.append(tri)
  if sid in allowed[course]:localouter.append(tri)
 for course in ['M3','main','twig']:
  ids=next(x['source_face_ids']for x in rows.values()if x['role']==course+'-outer');outside=[int(i)for i in ids if int(i)not in allowed[course]];assert all(len(seen[(course,i)])==1 and seen[(course,i)][0][1]for i in outside)
  proof.append(dict(course=course,source_faces=len(ids),explicit_modifiable_faces=len(allowed[course]),protected_outer_faces=len(outside),every_protected_face_present_once_with_literal_source_indices_and_positions=True))
 outer_side=Mesh(positions,np.array(fullouter),FRAME);ot=topology(outer_side);assert len(ot['components'])==1 and ot['open_edges']==48 and not ot['nonmanifold_edges']and not ot['inconsistent_winding_edges']and not ot['duplicate_faces']and not ot['zero_area_faces'];outer_side.save(OUT/'source-outer-assembly-repaired-both-roots.npz');Mesh(positions,np.array(localouter),FRAME).save(OUT/'both-roots-local-outer.npz')
 # Native outer source shading is transported by its literal source face.
 # Derived inner shading comes only from the actual derived inner chart.
 inputs={};smooth={}
 for item in g['native_inputs']:
  q=np.load(C/item['checkpoint_input']);v=q['positions_mm']/1000.+np.array(g['origin_m']);f=q['faces'];n=np.zeros_like(v);count=next(x['source_faces']for x in rows.values()if x['role']==item['role']);fs=f[:count];t=v[fs];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
  for k in range(3):np.add.at(n,fs[:,k],normal)
  length=np.linalg.norm(n,axis=1);n=np.divide(n,length[:,None],out=np.zeros_like(n),where=length[:,None]>0);inputs[item['role']]=(v,f);smooth[item['role']]=n
 wallfaces=[];corner=[];exact_corners=0;shadecounts=defaultdict(int)
 for tri,oid,fid in zip(*raw['wall']):
  row=rows[int(oid)];course=row['course'];role=row['role'];points=positions[tri]
  if int(fid)<row['source_faces']and role.endswith('outer'):
   sid=row['source_face_ids'][int(fid)];index=0 if course=='M3'else 1;src=source[index];fs=src.faces[sid];sq=src.vertices[fs];weights=barycentric(points,sq);norm=src.attributes['normals'][fs].astype(float)
   if 'normalSigns'in src.attributes:norm*=src.attributes['normalSigns'][fs,None]
   n=weights@norm;shadecounts['native_outer_interpolation']+=1
   for k,p in enumerate(points):
    match=np.flatnonzero(np.all(sq==p,axis=1))
    if len(match)==1:
     j=int(match[0]);tri[k]=int(fs[j])+(0 if index==0 else offset);n[k]=norm[j];exact_corners+=1
  elif int(fid)<row['source_faces']and role.endswith('fluid'):
   vp,fs=inputs[role];face=fs[int(fid)];weights=barycentric(points,vp[face]);n=-(weights@smooth[role][face]);shadecounts['derived_inner_interpolation']+=1
  else:
   n=np.broadcast_to(np.cross(points[1]-points[0],points[2]-points[0]),(3,3)).copy();shadecounts['actual_terminal_rim_or_numerical_extension_face_normal']+=1
  lengths=np.linalg.norm(n,axis=1);assert np.isfinite(n).all()and np.all(lengths>0);corner.append(n/lengths[:,None]);wallfaces.append(tri)
 positions=bank.positions;wallfaces=np.array(wallfaces);wall=Mesh(positions,wallfaces,FRAME);wt=topology(wall);assert wt['closed_oriented_manifold']and len(wt['components'])==1 and wt['components'][0]['closed_orientable_genus']==3
 atomic_npz(OUT/'wall-with-source-and-derived-corner-normals.npz',vertices=positions,faces=wallfaces,corner_normals=np.array(corner),original_id=raw['wall'][1],source_face_id=raw['wall'][2],frame=np.array(FRAME));wall.save(OUT/'wall.npz')
 for role in ['outer','fluid']:
  faces,oid,fid=raw[role];Mesh(positions,faces,FRAME).save(OUT/(role+'.npz'));t=topology(Mesh(positions,faces,FRAME));assert t['closed_oriented_manifold']and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==0
 # Exactly four terminal annuli; no internal root fan or blood closure disk.
 forbidden=[];ports=[]
 for row in rows.values():
  own=raw['wall'][1]==row['original_id'];fid=raw['wall'][2];forbidden.extend(np.flatnonzero(own&np.isin(fid,row['junction_faces'])).tolist())
  if row['role'].endswith('fluid'):forbidden.extend(np.flatnonzero(own&np.isin(fid,row['disk_faces'])).tolist())
  else:
   for mouth in row['mouths']:
    if mouth['kind']!='external-outer-classifier':continue
    ids=np.flatnonzero(own&np.isin(fid,mouth['newFaceIds']));r=Mesh(positions,wallfaces[ids],FRAME);t=topology(r);loops=boundary_loops(r.faces);assert len(t['components'])==1 and t['components'][0]['euler']==0 and len(loops)==2;ports.append(dict(course=row['course'],mouth=mouth['mouth'],tissue_faces=len(ids),boundary_cycles=2))
 assert len(ports)==4 and not forbidden
 report=dict(schema='full-native-three-course-finite-source-assembly/v1',installed=False,source_replay=proof,source_prefix_original_vertices_and_attributes_retained=True,maximum_source_roundoff_snap_m=bank.maximum_snap,new_intersection_point_snaps=bank.new_snaps,outer_side_topology=ot,wall_topology=wt,native_terminal_annuli=ports,root_fan_or_fluid_disk_faces_in_wall=0,shading_faces=dict(shadecounts),literal_native_outer_normal_corners_copied=exact_corners,scope='Complete native M3 and both postcentral source courses. Educational derived finite wall; not all MCA branches or a clinical thickness model.',remaining=['Independent full self-crossing and connected-fluid path checks.','Local fixed-context checks for both junction domains.','Actual full surface/longitudinal sections before installation.'],elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'assembly-proof.json',report)
 root=z['branch_vertices'][z['main_rings'][180]].mean(0);main_axis=z['M3_vertices'][z['M3_rings'][41]].mean(0)-z['M3_vertices'][z['M3_rings'][33]].mean(0);branch_axis=z['branch_vertices'][z['main_rings'][176]].mean(0)-root;plane=np.cross(main_axis,branch_axis);plane/=np.linalg.norm(plane);up=main_axis/np.linalg.norm(main_axis)
 spec=dict(frame=FRAME,wall_seconds=40,resolution=1000,meshes=[dict(name='Complete native three-course finite wall — derived educational thickness',mesh=str(OUT/'wall-with-source-and-derived-corner-normals.npz'),color=[.72,.15,.12,1.])],views=[dict(name='whole-front',direction=[0,0,1]),dict(name='whole-side',direction=[1,0,0]),dict(name='first-root-longitudinal-section',direction=plane.tolist(),up=up.tolist(),section=dict(point=root.tolist(),normal=plane.tolist(),clear_outer=True))]);atomic_json(OUT/'render-spec.json',spec);print(json.dumps(dict(wall_faces=len(wallfaces),protected_faces=sum(x['protected_outer_faces']for x in proof),terminal_annuli=len(ports),wall_genus=3,elapsed_seconds=report['elapsed_seconds'])))

if __name__=='__main__':main()
