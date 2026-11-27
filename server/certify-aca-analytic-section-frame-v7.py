"""Complete one analytic tangent section-frame repair; never fit or install."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
import open3d as o3d
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_sections import arclength,path_tangents,section_metrics,bend_diagnostics
W=ROOT/'.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1';O=W/'completed-feasible-course-v7';O.mkdir(exist_ok=False);started=time.monotonic();origin=np.array([0.,1.725,.028]);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();info=json.loads((W/'source-inputs.json').read_text());cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert sha(cat)==info['source_sha256']
def save(name,obj):(O/name).write_text(json.dumps(obj,indent=2,default=lambda x:x.tolist()if hasattr(x,'tolist')else str(x))+'\n')
def scene(m):
 sc=o3d.t.geometry.RaycastingScene(nthreads=1);sc.add_triangles(o3d.core.Tensor(((m.vertices-origin)*1000).astype(np.float32)),o3d.core.Tensor(m.faces.astype(np.uint32)));return sc
manifest=json.loads((W/'closed-sign-tissues-v4/closed-tissue-manifest.json').read_text());signed=[]
for row in manifest['rows']:
 if not row['closed']:continue
 m=Mesh.load(row['path']);assert m.geometry_sha==row['geometry_sha256'];signed.append((row['id'],scene(m),m.bounds))
context=json.loads((W/'active-fixed-context-v1/context-contract.json').read_text());template=np.load(W/'superior-course-candidate-v2-clamped/source-and-candidate-fields.npz');seed=np.load(W/'nonvascular-graph-seed-v5/nonvascular-lipschitz-polyline-seeds.npz');checks={};isometry={};fields_out={};meshes={};allpass=True
for side,suffix in [('left','l'),('right','r')]:
 proof=json.loads((W/f'graph-trust-v6/{side}/solve-proof.json').read_text());assert proof['completed']and proof['elapsed_seconds']<32;assert proof['minimum_sampled_clearance_mm']>=.2 and proof['maximum_middle_curvature_support']<=.85
 z=np.load(W/f'graph-trust-v6/{side}/diagnostic-checkpoint.npz');original=template[side+'_source_fields'];old=original[0].mean(1);rings=template[side+'_rings'];fixed=np.r_[np.arange(3),np.arange(23,37)];initial=seed[side+'_centres'];knots=np.array([2,5,8,10,13,16,20,23]);free=np.array([1,2,3,5,6]);p=initial[knots].copy();p[free]+=z['x'].reshape(-1,3)/1000;parameter=arclength(initial[2:24]);beforeT=path_tangents(old);curve=CubicSpline(parameter[knots-2],p,axis=0,bc_type=((1,beforeT[2]),(1,beforeT[23])));centres=old.copy();centres[2:24]=curve(parameter);centres[[2,13,23]]=np.array([old[2],[.0035 if suffix=='l'else-.0035,1.723,.0285],old[23]]);assert np.array_equal(centres,z['centres']),'Centreline must remain exactly unchanged'
 afterT=beforeT.copy();analytic=curve(parameter,1);analytic/=np.linalg.norm(analytic,axis=1)[:,None];afterT[2:24]=analytic;cross=np.cross(beforeT,afterT);dot=np.sum(beforeT*afterT,axis=1);rot=np.broadcast_to(np.eye(3),(37,3,3)).copy()
 for i,(v,c)in enumerate(zip(cross,dot)):
  assert c>-.999999;x,y,zz=v;skew=np.array([[0,-zz,y],[zz,0,-x],[-y,x,0]]);rot[i]+=skew+skew@skew/(1+c)
 fields=np.einsum('sij,ksvj->ksvi',rot,original-old[None,:,None])+centres[None,:,None];fields[:,fixed]=original[:,fixed];rot[fixed]=np.eye(3);assert np.array_equal(fields[:,fixed],original[:,fixed]);a,b=section_metrics(original),section_metrics(fields);area=float(abs(b['area']/a['area']-1).max());shape=float(abs(b['shape_distances']-a['shape_distances']).max());assert area<1e-11 and shape<1e-12
 row=next(r for r in info['sources']if r['part_id']==f'ZA:Anterior cerebral artery.{suffix}');source=Mesh.load(W/row['mesh_file']);v=source.vertices.copy();v[rings]=fields[0];attrs={k:a.copy()for k,a in source.attributes.items()};n=np.einsum('sij,svj->svi',rot,source.attributes['normals'][rings].astype(float));n/=np.linalg.norm(n,axis=2)[...,None];n=np.rint(n*32767).astype(np.int16);n[fixed]=source.attributes['normals'][rings[fixed]];attrs['normals'][rings]=n;m=Mesh(v,source.faces,source.frame,attrs,dict(source_id=row['part_id'],source_catalogue_sha256=sha(cat),finite_wall=False,analytic_section_orientation_repair=True));check=audit(m,require_closed=False,expected_components=1,max_pairs=250000,max_tests=1000000);checks[side]=check
 vv=[];ff=[];offset=0
 for entry in context['neighbors']:
  if any(word in entry['id']for word in ['Optic nerve','Optic chiasm','Olfactory nerve','Straight gyrus']):continue
  cm=Mesh.load(entry['path'])
  if entry['id']==f'ZA:Middle cerebral artery (M1-segment).{suffix}':
   table=np.load(W/f'native-carotid-terminal-cohort-v7b/mca-{suffix}-stations.npz');whole=Mesh.load(W/f'native-carotid-terminal-cohort-v7b/mca-{suffix}.npz');station=np.zeros(len(whole.vertices),int);station[table['rings']]=np.arange(len(table['rings']))[:,None];keep=station[whole.faces].min(1)[np.asarray(entry['source_face_ids'])]>=18;cm=Mesh(cm.vertices,cm.faces[keep],cm.frame)
   if not len(cm.faces):continue
  vv.append(cm.vertices);ff.append(cm.faces+offset);offset+=len(cm.vertices)
 other=scene(Mesh(np.concatenate(vv),np.concatenate(ff),source.frame));station=np.zeros(len(v),int);station[rings]=np.arange(37)[:,None];faces=source.faces[(station[source.faces].min(1)>=2)&(station[source.faces].max(1)<=23)];edges=np.unique(np.sort(np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]],axis=1),axis=0);used=np.unique(faces);samples=np.r_[v[used],v[edges].mean(1),v[faces].mean(1)];distance=[other.compute_distance(o3d.core.Tensor(((samples-origin)*1000).astype(np.float32)),nthreads=1).numpy().astype(float)]
 for pid,sc,box in signed:
  d=np.full(len(samples),1.);near=np.all(samples>=box[0]-.0003,axis=1)&np.all(samples<=box[1]+.0003,axis=1)
  if near.any():d[near]=sc.compute_signed_distance(o3d.core.Tensor(((samples[near]-origin)*1000).astype(np.float32)),nthreads=1,nsamples=3).numpy().astype(float)
  distance.append(d)
 distance=np.array(distance);support=bend_diagnostics(centres,fields,tangents=afterT)['curvature_support'];samplepass=bool(distance.min()>=.2);check['sampled_signed_clearance_mm']=float(distance.min());check['sampled_clearance_passed']=samplepass;check['analytic_middle_curvature_support_max']=float(np.asarray(support)[2:24].max());allpass=allpass and check['passed']and samplepass
 m.save(O/f'aca-{side}.npz');meshes[side]=m
 isometry[side]=dict(maximum_relative_section_area_error=area,maximum_pairwise_shape_error_m=shape,source_fixed_profiles_normals_UVs_literal=True,centreline_literal_v6=True,source_calibre_held=True,optimizer_energy_converged=False,geometry_frame_stage_complete=True,analytic_curve_endpoint_tangent_matches_native=True)
 for key,val in dict(rings=rings,source_fields=original,fields=fields,centres=centres,rotations=rot,source_normals=source.attributes['normals'],new_normals=attrs['normals'],analytic_tangents=afterT).items():fields_out[side+'_'+key]=val
 save('actual-geometry-gates.json',dict(complete=False,checks=checks,passed=False,installed=False));print(side,'self',check['self_intersections']['defect_pairs'],'clearance',float(distance.min()),'areaerror',area)
np.savez_compressed(O/'source-and-candidate-fields.npz',**fields_out);save('section-isometry-proof.json',isometry);save('actual-geometry-gates.json',dict(complete=True,checks=checks,passed=allpass,installed=False,anatomical_acceptance=False,source_catalogue_unchanged=sha(cat)==info['source_sha256'],elapsed_seconds=time.monotonic()-started,intended_native_terminal_root_union_required=True));save('candidate-manifest.json',dict(complete=True,source_catalogue_sha256=sha(cat),frame=meshes['left'].frame,passed_self_and_sampled_constraints=allpass,installed=False,geometry_edit='Source-to-analytic-C1-tangent proper section rotations only; literal centres and held section calibre; fixed0..2/23..36 source fields/normalwords/UVs unchanged.',full_context_union_gate_required=True,meshes=[dict(part_id=f'ZA:Anterior cerebral artery.{s}',candidate=f'aca-{side}.npz')for side,s in [('left','l'),('right','r')]]));print('passed',allpass,'elapsed',time.monotonic()-started)
