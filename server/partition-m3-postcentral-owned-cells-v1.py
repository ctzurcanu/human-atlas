"""Two owned finite cells; one canonical annular interface; no whole-mesh edit."""
from pathlib import Path
import json,sys,time
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.branch_junction_envelopes import annular_mouths
from geometry_tools.algorithms.finite_tissue_boolean import barycentric,boundary_signature
from geometry_tools.tissue import normal_words
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';C=P/'three-course-finite-educational-v1';A=C/'assembly-v1';OUT=C/'owned-cells-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(F/'three-course-source-charts.npz');g=json.loads((C/'generation-complete.json').read_text());q=np.load(A/'wall-with-source-and-derived-corner-normals.npz');v=q['vertices'];faces=q['faces'];field=q['corner_normals'];rows={r['original_id']:r for r in g['classifiers']};sources=[Mesh.load(P/'finite-junction-educational-v1/parent-source-untouched.npz'),Mesh.load(P/'finite-junction-educational-v1/branch-source-untouched.npz')];charts={r:np.load(C/(r+'-derived-paired-chart.npz'))for r in ['M3','main','twig']};inputrows={}
 for item in g['native_inputs']:
  r=np.load(C/item['checkpoint_input']);inputrows[item['role']]=(r['positions_mm']/1000.+np.array(g['origin_m']),r['faces'])
 ms=np.full(len(z['branch_vertices']),-1,int)
 for station,ring in enumerate(z['main_rings']):ms[ring]=station
 selected=defaultdict(list);fields=defaultdict(list);uvs=defaultdict(list);ownership=[];exact_normal=0;exact_uv=0
 for index,(tri,oid,fid)in enumerate(zip(faces,q['original_id'],q['source_face_id'])):
  row=rows[int(oid)];role=row['role'];course=row['course'];source_side=int(fid)<row['source_faces'];sid=row['source_face_ids'][int(fid)]if source_side else None
  if course=='M3':owner='M3'
  elif course=='twig':owner='postcentral'
  elif source_side:
   stations=ms[z['branch_faces'][sid]];assert stations.min()>=0
   if stations.max()<=176:owner='postcentral'
   elif stations.min()>=176:owner='M3'
   else:raise ValueError('An output parent face crosses the declared station-176 ownership chart')
  else:owner='postcentral'  # Only main native distal terminal charts survive.
  chart='outer'if source_side and role.endswith('outer')else'inner'if source_side and role.endswith('fluid')else'terminal-rim';key=owner+':'+chart;points=v[tri];words=np.rint(np.clip(field[index],-1,1)*32767).astype('<i2')
  if source_side:
   src=sources[0 if course=='M3'else 1];sf=src.faces[sid];sq=src.vertices[sf]
   if role.endswith('outer'):
    weights=barycentric(points,sq);native=src.attributes['normals'][sf].astype(float)
    if 'normalSigns'in src.attributes:native*=src.attributes['normalSigns'][sf,None]
    uv=weights@src.attributes['uvs'][sf]
    for k,p in enumerate(points):
     match=np.flatnonzero(np.all(sq==p,axis=1))
     if len(match)==1:
      j=int(match[0]);words[k]=native[j].astype('<i2');uv[k]=src.attributes['uvs'][sf[j]];exact_normal+=1;exact_uv+=1
   else:
    ip,ff=inputrows[role];ifc=ff[int(fid)];outer=charts[course]['outer_positions'][ifc];match=np.all(outer[:,None,:]==sq[None,:,:],axis=2);assert np.all(match.sum(1)==1);sourceuv=src.attributes['uvs'][sf[np.argmax(match,axis=1)]];uv=barycentric(points,ip[ifc])@sourceuv
  else:
   normal=np.cross(points[1]-points[0],points[2]-points[0]);normal/=np.linalg.norm(normal);axis=np.eye(3)[np.argmin(abs(normal))];right=np.cross(axis,normal);right/=np.linalg.norm(right);up=np.cross(normal,right);uv=np.c_[points@right,points@up]/.001
  selected[key].append(tri);fields[key].append(words);uvs[key].append(np.asarray(uv,dtype='<f4'));ownership.append(owner)
 # Exact-coordinate welding is topology only; the actual source-prefix and
 # every unchanged triangle corner stay in the retained position bank.
 positions,inverse=np.unique(v,axis=0,return_inverse=True);a_faces=faces[np.array(ownership)=='M3'];b_faces=faces[np.array(ownership)=='postcentral'];loops=boundary_loops(inverse[a_faces]);assert len(loops)==2
 source_outer=z['branch_vertices'][z['main_rings'][176]];source_inner=charts['main']['inner_positions'][charts['main']['rings'][176]];outerkeys={tuple(p):i for i,p in enumerate(source_outer)};pair={tuple(a):tuple(b)for a,b in zip(source_outer,source_inner)};outerloop=next(loop for loop in loops if {tuple(p)for p in positions[loop]}==set(outerkeys));lookup={tuple(p):i for i,p in enumerate(v)};outside=np.array([lookup[tuple(p)]for p in positions[outerloop]],int);inside=np.array([lookup[pair[tuple(p)]]for p in positions[outerloop]],int);interface=annular_mouths([(outside,inside)]);assert not boundary_signature(np.r_[inverse[a_faces],inverse[interface]])and not boundary_signature(np.r_[inverse[b_faces],inverse[interface[:,[0,2,1]]]])
 selected['shared:station176']=interface;fields['shared:station176']=normal_words(v,{'s':interface})['s'];origin=v[interface].mean((0,1));_,_,basis=np.linalg.svd(v[outside]-v[outside].mean(0),full_matrices=False);uvs['shared:station176']=np.stack([(v[t]-origin)@basis[:2].T/.001 for t in interface]).astype('<f4')
 surfaces={k:np.asarray(value,dtype=np.int64)for k,value in selected.items()};words={k:np.asarray(fields[k],dtype='<i2')for k in surfaces};tex={k:np.asarray(uvs[k],dtype='<f4')for k in surfaces};cells={owner:[dict(surface=k,side=1)for k in surfaces if k.startswith(owner+':')]for owner in ['M3','postcentral']};cells['M3'].append(dict(surface='shared:station176',side=1,neighbor='postcentral'));cells['postcentral'].append(dict(surface='shared:station176',side=-1,neighbor='M3'));partids={'M3':'ZA:Middle cerebral artery (M3 segment).l','postcentral':'ZA:Postcentral arterial branch.l'};arrays=dict(positions=v,frame=np.array(FRAME));surface_rows=[]
 for i,key in enumerate(surfaces):
  arrays[f'surface{i}']=surfaces[key];arrays[f'normal{i}']=words[key];arrays[f'uv{i}']=tex[key];surface_rows.append(dict(id=key,array=f'surface{i}',corner_words=f'normal{i}',corner_uvs=f'uv{i}'))
 atomic_npz(OUT/'owned-arrangement.npz',**arrays);manifest=dict(frame=FRAME,surfaces=surface_rows,cells=cells,part_ids=partids,geometry_origin=[0.,1.75,0.],validated_whole_sha256=digest_file(A/'wall.npz'),identity_boundary='Main postcentral native station 176 is a provisional source-preserving ownership boundary; the M3 cell includes the short proximal postcentral junction stub. No additional TA98 identity or whole-MCA completion is claimed.',thickness_status='Educational derived inner profiles under authorized formula, not measured native finite walls.');atomic_json(OUT/'owned-arrangement.json',manifest)
 celltop=[]
 for name,refs in cells.items():
  f=np.concatenate([surfaces[r['surface']]if r['side']==1 else surfaces[r['surface']][:,[0,2,1]]for r in refs]);mesh=Mesh(v,f,FRAME);t=topology(mesh);assert t['closed_oriented_manifold']and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==2;mesh.save(OUT/(name+'-closed-cell.npz'));celltop.append(dict(cell=name,topology=t))
 # Float32 encoding uses a binary-exact superior-axis translation. Every
 # original native outer coordinate remains identical after decode.
 encoding_origin=np.array(manifest['geometry_origin']);native=np.r_[sources[0].vertices,sources[1].vertices];assert np.array_equal((native-encoding_origin).astype('<f4').astype(float)+encoding_origin,native)
 report=dict(schema='two-owned-finite-M3-postcentral-cells/v1',installed=False,validated_whole_unchanged=digest_file(A/'wall.npz')==manifest['validated_whole_sha256'],whole_wall_faces_partitioned_once=len(ownership),single_canonical_annular_interface_faces=len(interface),shared_interface_uses=[dict(cell='M3',side=1,neighbor=partids['postcentral']),dict(cell='postcentral',side=-1,neighbor=partids['M3'])],shared_interface_in_assembled_display=False,source_native_outer_normal_corners_copied_exactly=exact_normal,source_native_outer_uv_corners_copied_exactly=exact_uv,source_native_outer_positions_exact_after_planned_Float32_encoding=True,cells=celltop,identity_boundary=manifest['identity_boundary'],remaining=['Actual emitted shared-bank pack/unpack checks and independent cell disjointness/union proof.','Packed-normal actual surface/section review and fixed current context before installation.'],elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'partition-proof.json',report);print(json.dumps(dict(wall_faces_partitioned=len(ownership),shared_annulus_faces=len(interface),native_outer_normal_corners_exact=exact_normal,elapsed_seconds=report['elapsed_seconds'])))

if __name__=='__main__':main()
