"""Explicit contained cardiac tissue fields; no normal-offset assumptions.

Coating constraints are independent of regional ownership, so chamber or
septal boundaries cannot change layer allocation. All calculation tetrahedra
remain inside the audited wall. Regional extents are provisional source-
guided constructions and require the actual rendered/section review.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='tissue-fields-v30');ap.add_argument('--calculation',default='wall-calculation-v26');ap.add_argument('--edge-sampled',action='store_true');args=ap.parse_args()
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';out=base/args.stage;out.mkdir(exist_ok=False);source=base/'whole-wall-v25/whole-wall.npz';m=Mesh.load(source);calc=base/args.calculation;assert json.loads((calc/'construction.json').read_text())['passed'];z=np.load(calc/'tetrahedra.npz');p=z['positions'];original=z['tetrahedra']
if args.edge_sampled:t=original
else:
 centres=p[original].mean(1);centre_ids=np.arange(len(original))+len(p);p=np.r_[p,centres];t=np.concatenate([np.c_[centre_ids,original[:,face]]for face in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]])
assert len(t)<1000000;o=np.array([0.,1.375,.03125]);charts=np.load(base/'wall-regions-v27/charts.npz');inner=charts['endocardial_faces'];owner=charts['face_owner'];q=m.triangles;normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];dist={};closest={};normals={}
for kind,condition in [('endo',inner),('epi',~inner)]:
 for chamber in range(4):
  mask=np.flatnonzero(condition&(owner==chamber));assert len(mask);scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces[mask].astype('uint32')));hit=scene.compute_closest_points(o3d.core.Tensor((p-o).astype('float32')),nthreads=1);points=hit['points'].numpy().astype(float)+o;distance=np.linalg.norm(p-points,axis=1);dist[kind,chamber]=np.rint(distance/.000002)*.000002;closest[kind,chamber]=points;normals[kind,chamber]=normal[mask[hit['primitive_ids'].numpy()]]
de=np.column_stack([dist['endo',i]for i in range(4)]);dp=np.column_stack([dist['epi',i]for i in range(4)]);DE=de.min(1);DP=dp.min(1);span=DE+DP;thickness=np.minimum(.00004,.15*span);fe=thickness-DE;fp=thickness-DP;assert np.all(fe+fp<=1e-16)
def septum_guidance(i,j,limit):
 a=closest['endo',i];b=closest['endo',j];direction=b-a;length=np.linalg.norm(direction,axis=1);unit=direction/np.maximum(length,1e-30)[:,None];projection=np.einsum('ij,ij->i',p-a,unit);opposition=np.einsum('ij,ij->i',normals['endo',i],normals['endo',j]);g=np.minimum.reduce([.5*(limit-de[:,i]-de[:,j]),projection,length-projection,.003*(-.3-opposition)]);g[length<.001]=np.minimum(g[length<.001],-.01);return g
baseline=-DE;g_ias=septum_guidance(0,2,.012);g_ivs=septum_guidance(1,3,.025);g_av=septum_guidance(0,3,.010);core=np.c_[-de,baseline+g_ias,baseline+g_ivs,baseline+g_av];core=np.rint(core/.000001)*.000001;source_count=int(z['original_positions'])if args.edge_sampled else len(z['positions']);candidate=np.flatnonzero((np.argmax(core,axis=1)==5)&(np.arange(len(p))>=source_count));assert len(candidate),'No source-guided interventricular septal tissue'
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());valve=next(x for x in a['parts']if x['id']=='ZA:Non-coronary leaflet');vm=AtlasReader().read_part(a,valve,frame=m.frame,catalogue_path=path);apex=m.vertices[np.argmin(m.vertices[:,1]-.3*m.vertices[:,0])];direction=apex-vm.vertices.mean(0);landmark=vm.vertices.mean(0)+direction/np.linalg.norm(direction)*.008;centre_index=int(candidate[np.argmin(np.linalg.norm(p[candidate]-landmark,axis=1))]);membrane_center=p[centre_index];membrane=.003-np.linalg.norm(p-membrane_center,axis=1)
fields=[fe,fp];field_names=['endocardium-internal-allocation','epicardium-external-allocation'];zero=[-1,-1]
def plane(name,values,priority=1):
 index=len(fields);fields.append(values);field_names.append(name);zero.append(priority);return index
core_planes={};endo_planes={};epi_planes={}
for i in range(7):
 for j in range(i+1,7):core_planes[i,j]=plane(f'myocardial-region:{i}:{j}',core[:,i]-core[:,j])
for kind,values,table in [('endo',de,endo_planes),('epi',dp,epi_planes)]:
 for i in range(4):
  for j in range(i+1,4):table[i,j]=plane(f'{kind}-chamber:{i}:{j}',values[:,j]-values[:,i])
membrane_plane=plane('membranous-versus-muscular-septum',membrane,priority=-1);names=[];cells=[]
def comparisons(index,count,table):return [(table[min(index,j),max(index,j)],1 if index<j else-1)for j in range(count)if j!=index]
for i,key in enumerate(['right-atrium','right-ventricle','left-atrium','left-ventricle','interatrial-septum','interventricular-muscular-septum','atrioventricular-septum']):
 refs=[(0,-1),(1,-1),*comparisons(i,7,core_planes)]
 if i==5:refs.append((membrane_plane,-1))
 names.append('myocardium:'+key);cells.append(refs)
names.append('fibrous:interventricular-membranous-septum');cells.append([(0,-1),(1,-1),*comparisons(5,7,core_planes),(membrane_plane,1)])
for kind,table,allocation in [('endocardium',endo_planes,0),('epicardium',epi_planes,1)]:
 for i,key in enumerate(['right-atrium','right-ventricle','left-atrium','left-ventricle']):names.append(kind+':'+key);cells.append([(allocation,1),*comparisons(i,4,table)])
atomic_npz(out/'fields.npz',positions=p,tetrahedra=t,fields=np.column_stack(fields),zero_side_priority=np.array(zero),origin=o);atomic_json(out/'cells.json',dict(names=names,constraints=cells,field_names=field_names));report=dict(prepared=True,source_geometry_sha256=m.geometry_sha,source_archive_sha256=digest_file(source),source_calculation_sha256=digest_file(calc/'tetrahedra.npz'),source_charts_sha256=digest_file(base/'wall-regions-v27/charts.npz'),positions=len(p),calculation_tetrahedra=len(t),all_internal_calculation_edges_sampled=args.edge_sampled,cells=names,planes=len(fields),coating_target_m=.00004,coating_maximum_local_span_fraction=.15,distance_field_quantum_m=.000002,regional_score_quantum_m=.000001,membranous_region_center_m=membrane_center.tolist(),membranous_guide_landmark_m=landmark.tolist(),membranous_region_radius_m=.003,membranous_guide_residual_m=float(np.linalg.norm(membrane_center-landmark)),native_atrium_ventricle_score_winners=np.bincount(np.argmax(core,axis=1),minlength=7).tolist(),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
