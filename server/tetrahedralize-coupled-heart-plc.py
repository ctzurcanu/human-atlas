"""One PLC keeps all three finite tissues conforming when edges refine."""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
import open3d as o3d
import meshpy.tet as tet
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.annular_pipe_envelopes import face_components,signed_volume
ap=argparse.ArgumentParser();ap.add_argument('--source',default='chart-coats-v75');ap.add_argument('--stage',default='tissue-domains-v81');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/args.source;out=b/args.stage;out.mkdir(exist_ok=False);z=np.load(source/'cells.npz');p=z['positions'];o=np.array([0.,1.375,.03125]);bank={};regions=[];source_volumes=[];names=['myocardium','endocardium','epicardium']
for domain,name in enumerate(names):
 f=z['cell_'+str(domain)];source_volumes.append(signed_volume(p,f));nc,labels=face_components(f);scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((p-o).astype('float32')),o3d.core.Tensor(f.astype('uint32')))
 for face in f:
  key=tuple(sorted(map(int,face)));bank.setdefault(key,[]).append((domain,face))
 for component in range(nc):
  ff=f[labels==component]
  if signed_volume(p,ff)<0:continue
  q=p[ff];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);area=np.linalg.norm(normal,axis=1);ids=np.argsort(area)[-32:];normal=normal[ids]/area[ids,None];centres=q[ids].mean(1);samples=np.concatenate([centres-normal*eps for eps in [.0000001,.0000003,.000001,.000003,.00001]]);d=scene.compute_signed_distance(o3d.core.Tensor((samples-o).astype('float32')),nthreads=1,nsamples=3).numpy();best=int(np.argmin(d));assert d[best]<-.00000004,(name,component,d.min());point=samples[best];regions.append([*((point-o)*1000),domain+1,-1.])
 assert any(r[3]==domain+1 for r in regions);print(name,'material seeds',sum(r[3]==domain+1 for r in regions),flush=True)
facets=[];markers=[]
for uses in bank.values():
 assert len(uses)<=2
 if len(uses)==2:
  f0,f1=uses[0][1],uses[1][1];assert any(np.array_equal(np.roll(f0,k),f1[[0,2,1]])for k in range(3))
 facets.append(uses[0][1]);owners=tuple(sorted(row[0]for row in uses));markers.append({(0,):1,(1,):2,(2,):3,(0,1):4,(0,2):5,(1,2):6}[owners])
facets=np.array(facets);used,inv=np.unique(facets,return_inverse=True);p=p[used];facets=inv.reshape(-1,3);info=tet.MeshInfo();info.set_points(((p-o)*1000).tolist());info.set_facets(facets.tolist(),markers=markers);info.regions.resize(len(regions))
for i,r in enumerate(regions):info.regions[i]=r
options=tet.Options('pA');options.nomergefacet=1;options.nomergevertex=1;options.steinerleft=100000;result=tet.build(info,options=options,attributes=True);v=np.array(result.points)/1000+o;t=np.array(result.elements,dtype='int64');attributes=np.asarray(result.element_attributes).reshape(len(t),-1);assert attributes.shape[1]==1;domains=attributes[:,0].astype('int64')-1;assert np.isin(domains,[0,1,2]).all();assert np.max(np.linalg.norm(v[:len(p)]-p,axis=1))<1e-12;v[:len(p)]=p;q=v[t];det=np.linalg.det(np.transpose(q[:,1:]-q[:,:1],(0,2,1)));assert np.all(det!=0);negative=det<0;t[negative]=t[negative][:,[0,2,1,3]];volumes=np.array([abs(det[domains==i]).sum()/6 for i in range(3)]);assert np.all(abs(volumes-source_volumes)<np.maximum(1e-15,abs(np.array(source_volumes))*1e-8)),(volumes,source_volumes);atomic_npz(out/'tetrahedra.npz',positions=v,tetrahedra=t,domains=domains,original_position_count=np.array(len(p)),origin=o);report=dict(passed=True,source_sha256=digest_file(source/'cells.npz'),source_points_literal=True,shared_interfaces_in_one_plc=True,source_volumes_mL=(np.array(source_volumes)*1e6).tolist(),calculation_volumes_mL=(volumes*1e6).tolist(),positions=len(v),added_steiner_points=len(v)-len(p),tetrahedra=len(t),material_seed_count=len(regions),source_facets=len(facets),elapsed_seconds=time.monotonic()-start,construction_only=True,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
