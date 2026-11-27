"""Conforming calculation domains from actual finite heart tissues."""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
import meshpy.tet as tet
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.annular_pipe_envelopes import face_components,signed_volume
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='tissue-domains-v76');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'chart-coats-v75';out=b/args.stage;out.mkdir(exist_ok=False);z=np.load(source/'cells.npz');p=z['positions'];o=np.array([0.,1.375,.03125]);global_points=list(p);all_tets=[];all_domains=[];rows=[]
for domain,name in enumerate(['myocardium','endocardium','epicardium']):
 f=z['cell_'+str(domain)];used,inv=np.unique(f,return_inverse=True);local_f=inv.reshape(-1,3);local_p=p[used];info=tet.MeshInfo();info.set_points(((local_p-o)*1000).tolist());info.set_facets(local_f.tolist(),markers=list(range(1,len(local_f)+1)));nc,label=face_components(local_f);holes=[]
 for component in range(nc):
  ff=local_f[label==component]
  if signed_volume(local_p,ff)>=0:continue
  q=local_p[ff];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);area=np.linalg.norm(normal,axis=1);i=int(np.argmax(area));normal=normal[i]/area[i];point=q[i].mean(0)+normal*.000001;holes.append(((point-o)*1000).tolist())
 if holes:info.set_holes(holes)
 options=tet.Options('pQYY');options.nomergefacet=1;options.nomergevertex=1;options.nobisect=2;result=tet.build(info,options=options);v=np.asarray(result.points)/1000+o;t=np.asarray(result.elements,dtype='int64');assert np.max(np.linalg.norm(v[:len(used)]-local_p,axis=1))<1e-12;v[:len(used)]=local_p;mapping=np.r_[used,np.arange(len(global_points),len(global_points)+len(v)-len(used))];global_points.extend(v[len(used):]);q=v[t];det=np.linalg.det(np.transpose(q[:,1:]-q[:,:1],(0,2,1)));assert np.all(det!=0);neg=det<0;t[neg]=t[neg][:,[0,2,1,3]];volume=float(abs(det).sum()/6);reference=signed_volume(local_p,local_f);assert abs(volume-reference)<max(1e-15,abs(reference)*1e-9),(name,volume,reference)
 faces=np.concatenate([t[:,ii]for ii in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]]);keys,counts=np.unique(np.sort(faces,axis=1),axis=0,return_counts=True);assert counts.max()==2;boundary=keys[counts==1];expected=np.sort(local_f,axis=1);expected=expected[np.lexsort(expected.T[::-1])];atomic_npz(out/('domain-'+str(domain)+'-calculation.npz'),positions=v,tetrahedra=t,source_faces=local_f,boundary_faces=boundary,global_source_ids=used);assert np.array_equal(boundary,expected),(name,'boundary facets changed',len(boundary),len(expected))
 all_tets.append(mapping[t]);all_domains.append(np.full(len(t),domain,dtype='int64'));row=dict(name=name,domain=domain,source_triangles=len(f),tetrahedra=len(t),holes=len(holes),source_faces_literal=True,source_points_literal=True,tissue_volume_mL=volume*1e6,reference_volume_mL=reference*1e6);rows.append(row);atomic_json(out/'progress.json',dict(domains=rows,elapsed_seconds=time.monotonic()-start),overwrite=True);print(row,flush=True)
atomic_npz(out/'tetrahedra.npz',positions=np.array(global_points),tetrahedra=np.concatenate(all_tets),domains=np.concatenate(all_domains),original_position_count=np.array(len(p)),origin=o);report=dict(passed=True,source_sha256=digest_file(source/'cells.npz'),source_geometry_sha256=json.loads((source/'construction.json').read_text())['source_geometry_sha256'],domains=rows,positions=len(global_points),tetrahedra=sum(len(t)for t in all_tets),all_source_boundaries_literal=True,elapsed_seconds=time.monotonic()-start,construction_only=True,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
