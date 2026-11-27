"""Source-bound finite coats with explicit taper at singular source vertices.

No nodal extrapolation into the myocardium. Singular carrier stars remain
uncoated and explicit; the literal outer anatomy is never moved.
"""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
import open3d as o3d
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
from geometry_tools.algorithms.canonical_wall_rims import canonical_wall_rims
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='chart-coats-v58');ap.add_argument('--source',default='whole-wall-v25');ap.add_argument('--charts',default='wall-regions-v27');ap.add_argument('--finite-chart-ends',action='store_true');ap.add_argument('--exclude-diagnosed-core',default=None);args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);m=Mesh.load(b/args.source/'whole-wall.npz');o=np.array([0.,1.375,.03125]);p=m.vertices;f=m.faces;q=m.triangles;cross=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);area=np.linalg.norm(cross,axis=1);normal=cross/area[:,None];sums=np.zeros_like(p);inc=[[]for _ in p]
for k in range(3):
 np.add.at(sums,f[:,k],cross)
 for face,vertex in enumerate(f[:,k]):inc[vertex].append(face)
directions=-sums/np.linalg.norm(sums,axis=1)[:,None];bad=[];corrected=[];carrier_height=area/np.linalg.norm(q-np.roll(q,1,axis=1),axis=2).max(1);minimum_height=np.full(len(p),np.inf)
if args.exclude_diagnosed_core:
 from scipy.spatial import cKDTree
 old=b/args.exclude_diagnosed_core;metadata=json.loads((old/'construction.json').read_text());assert metadata['source_geometry_sha256']==m.geometry_sha;cm=Mesh.load(old/'cell-0.npz');cp=json.loads((old/'core-audit.json').read_text());pairs=cp['self_intersections']['examples'];assert len(pairs)==cp['self_intersections']['defect_pairs'];vertices=np.unique(cm.faces[np.array(pairs).ravel()]);closest=cKDTree(p).query(cm.vertices[vertices])[1];support=np.any(np.isin(f,closest),axis=1);bad.extend(np.unique(f[support]).tolist())
for k in range(3):np.minimum.at(minimum_height,f[:,k],carrier_height)
for i,faces in enumerate(inc):
 n=normal[faces]
 if i in bad:directions[i]=0;continue
 if np.all(n@directions[i]<-.02):continue
 lp=linprog([0.,0.,0.,-1.],A_ub=np.c_[n,np.ones(len(n))],b_ub=np.zeros(len(n)),bounds=[(-1,1)]*3+[(0,1)],method='highs')
 if not lp.success or np.any(n@lp.x[:3]>=-1e-10):bad.append(i);directions[i]=0
 else:
  unit=lp.x[:3]/np.linalg.norm(lp.x[:3])
  if np.min(-n@unit)<.03:bad.append(i);directions[i]=0
  else:directions[i]=unit;corrected.append(i)
scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((p-o).astype('float32')),o3d.core.Tensor(f.astype('uint32')));epsilon=.000001;initial=np.flatnonzero(np.linalg.norm(directions,axis=1)>0);rays=np.c_[p[initial]-o+directions[initial]*epsilon,directions[initial]].astype('float32');initial_hit=scene.cast_rays(o3d.core.Tensor(rays),nthreads=1)['t_hit'].numpy();unsupported=initial[(~np.isfinite(initial_hit))|(initial_hit<.00001)|(minimum_height[initial]<.00001)];bad=sorted(set(bad)|set(map(int,unsupported)));uncoated_faces=np.any(np.isin(f,bad),axis=1);fixed=np.array(bad,dtype=int)if args.finite_chart_ends else np.unique(f[uncoated_faces]);directions[fixed]=0;live=np.linalg.norm(directions,axis=1)>0;hit=np.zeros(len(p));hit[initial]=initial_hit;hit=hit[live];assert np.isfinite(hit).all()and np.all(hit>0);height=carrier_height;local=np.full(len(p),np.inf)
for k in range(3):np.minimum.at(local,f[:,k],height)
distance=np.zeros(len(p));distance[live]=np.minimum(.00004,np.minimum(.05*(hit+epsilon),.05*local[live]));inner=((p-o+directions*distance[:,None]).astype('float32').astype(float)+o);charts=np.load(b/args.charts/'charts.npz');endo=charts['endocardial_faces'];owner=charts['face_owner'];v=np.r_[p,inner];cells={};names=['myocardium']
def rim_for(selected):
 edges=np.concatenate([selected[:,[0,1]],selected[:,[1,2]],selected[:,[2,0]]]);canonical=np.sort(edges,axis=1);unique,first,counts=np.unique(canonical,axis=0,return_index=True,return_counts=True);assert counts.max()<=2;boundary=edges[first[counts==1]];return canonical_wall_rims(boundary,len(p))
if args.finite_chart_ends:
 labels=np.where(uncoated_faces,2,np.where(endo,0,1));edge_faces={}
 for i,face in enumerate(f):
  for k in range(3):edge_faces.setdefault(tuple(sorted((int(face[k]),int(face[(k+1)%3])))),[]).append(i)
 cycles=[]
 for vertex,faces in enumerate(inc):
  first=faces[0];order=[first];previous=-1;current=first
  while True:
   face=f[current];other=[int(k)for k in face if k!=vertex];neighbors=[next(k for k in edge_faces[tuple(sorted((vertex,j)))]if k!=current)for j in other];nxt=next(k for k in neighbors if k!=previous)
   if nxt==first:break
   assert nxt not in order;order.append(nxt);previous,current=current,nxt
  cycles.append(np.array(order))
 for iteration in range(32):
  changes=[]
  for vertex,cycle in enumerate(cycles):
   values=labels[cycle]
   if any(np.sum((values==kind)&(np.roll(values,1)!=kind))>1 for kind in range(3)):changes.extend(cycle.tolist())
  if not changes:break
  labels[np.unique(changes)]=2
 else:raise ValueError('Finite chart boundary conditioning did not converge')
 uncoated_faces=labels==2;assert area[uncoated_faces].sum()<.05*area.sum(),'More than 5% uncoated source support'
 selected=f[~uncoated_faces];rims=rim_for(selected);cells['myocardium']=np.r_[f[uncoated_faces],selected+len(p),rims[:,[0,2,1]]]
 for k,kind in enumerate(['endocardium','epicardium']):
  selected=f[labels==k];rims=rim_for(selected);cells[kind]=np.r_[selected,selected[:,[0,2,1]]+len(p),rims];names.append(kind)
else:
 cells['myocardium']=f+len(p)
 for kind,mask in [('endocardium',endo),('epicardium',~endo)]:
  for j,key in enumerate(['right-atrium','right-ventricle','left-atrium','left-ventricle']):
   selected=f[mask&(owner==j)];rims=rim_for(selected);name=kind+':'+key;cells[name]=np.r_[selected,selected[:,[0,2,1]]+len(p),rims];names.append(name)
print('Source-bound coats',len(bad),'unsupported vertices;',int(uncoated_faces.sum()),'uncoated source faces',flush=True)
r=normalize_shared_boundary(v,cells,o);atomic_json(out/'encoding.json',r['evidence']);whole=normalize_shared_boundary(r['positions'],{'whole':np.concatenate(list(r['cells'].values()))},o);exterior=Mesh(whole['positions'],whole['cells']['whole'],m.frame);assert exterior.surface_sha==m.surface_sha,'Finite tissue exterior differs from the original whole wall';arrays={'positions':r['positions'],'outer_source_vertices':np.array(len(p))};rows=[]
for i,name in enumerate(names):
 faces=r['cells'][name];arrays['cell_'+str(i)]=faces;ids,inv=np.unique(faces,return_inverse=True);mesh=Mesh(r['positions'][ids],inv.reshape(-1,3),m.frame);mesh.save(out/('cell-'+str(i)+'.npz'))
 if i==0:
  proof=audit(mesh,require_closed=True,max_pairs=2000000,max_tests=16000000);atomic_json(out/'core-audit.json',proof);print('Core',proof['passed'],proof['failures'],'crossings',proof['self_intersections']['defect_pairs'],flush=True)
 rows.append(dict(name=name,triangles=len(faces),geometry_sha256=mesh.geometry_sha))
atomic_npz(out/'cells.npz',**arrays);atomic_npz(out/'fields.npz',outer=p,inner=inner,faces=f,uncoated_faces=uncoated_faces,fixed_source_controls=fixed);report=dict(constructed=True,source_geometry_sha256=m.geometry_sha,cells=rows,names=names,nominal_coat_thickness_m=.00004,minimum_live_thickness_m=float(distance[live].min()),singular_source_vertices=bad,taper_source_vertices=fixed.tolist(),uncoated_source_faces=int(uncoated_faces.sum()),uncoated_surface_area_mm2=float(area[uncoated_faces].sum()*5e5),source_surface_area_mm2=float(area.sum()*5e5),positions_encoded_once=True,corrected_direction_count=len(corrected),core_audit_passed=proof['passed'],elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
