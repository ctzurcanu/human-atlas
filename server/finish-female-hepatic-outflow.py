"""Attach the three owned hepatic vein trees to the native female caval ports."""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.annular_bridge import zipper,arc,attach_bridge
from geometry_tools.algorithms.finite_vessel_wall import finite_wall
from geometry_tools.algorithms.finite_normal_field import inward_field
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/'female-outflow-v262';out.mkdir(exist_ok=False);start=time.monotonic();p=b/'final-package-v246/female.json';a=json.loads(p.read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame='model:local-female-ta98:assembled-metres'
def compact(m):
 used=np.unique(m.faces);v,inv=np.unique(m.vertices[used],axis=0,return_inverse=True);ix=np.full(len(m.vertices),-1,int);ix[used]=inv;return Mesh(v,ix[m.faces],frame)
def get(pid):return compact(r.read_part(a,parts[pid],frame=frame,catalogue_path=p))
branches=[];caps=[]
whole=get('TA98HEPATIC:female:hepatic_vein')
cava_source=get('LOCAL:female:vena_cava:1797');native_ports=[loop for loop in boundary_loops(cava_source.faces)if 1.23<cava_source.vertices[loop,1].mean()<1.26];native_centers=np.array([cava_source.vertices[loop].mean(0)for loop in native_ports]);root_centers=[]
for term in ['A12.3.09.006','A12.3.09.007','A12.3.09.008']:
 m=get('TA98HEPATIC:female:region:'+term);root_centers.append(np.concatenate([m.vertices[x]for x in boundary_loops(m.faces)]).mean(0))
_,port_assignment=linear_sum_assignment(np.linalg.norm(np.asarray(root_centers)[:,None]-native_centers[None],axis=2))
for term in ['A12.3.09.006','A12.3.09.007','A12.3.09.008']:
 m=get('TA98HEPATIC:female:region:'+term);ls=boundary_loops(m.faces);assert len(ls)==2;root=np.concatenate([m.vertices[x]for x in ls]).mean(0)
 j=port_assignment[len(branches)];xyz=cava_source.vertices[native_ports[j]];_,_,axes=np.linalg.svd(xyz-native_centers[j]);normal=axes[-1]
 if np.dot(normal,m.vertices[np.unique(m.faces)].mean(0)-native_centers[j])<0:normal=-normal
 offset=float(max(root@normal+.003,native_centers[j]@normal+.003))
 from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank,split_closed_cell
 from scipy.sparse import coo_matrix
 from scipy.sparse.csgraph import connected_components
 bank=SharedCutBank(whole.vertices);cut=split_closed_cell(bank,whole.faces,normal,offset,'hepatic-root:'+term,allow_holes=True);vv=np.asarray(bank.positions);ff=cut.positive_faces;ancestry=cut.positive_ancestry;e=np.r_[ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]];graph=coo_matrix((np.ones(len(e)*2),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(vv),len(vv))).tocsr();_,labels=connected_components(graph,directed=False);regions={tuple(sorted(map(tuple,t)))for t in m.triangles};named=np.array([tuple(sorted(map(tuple,t)))in regions for t in whole.triangles]);scores={int(k):int(np.count_nonzero((labels[ff[:,0]]==k)&(ancestry>=0)&named[np.maximum(ancestry,0)]))for k in np.unique(labels[ff[:,0]])};chosen=max(scores,key=scores.get);take=labels[ff[:,0]]==chosen;cap=np.flatnonzero(ancestry[take]<0);closed=Mesh(vv,ff[take],frame);assert len(boundary_loops(closed.faces[cap]))==2;closed.save(out/(term+'-tree.npz'));check=audit(closed,require_closed=True,max_pairs=500000,max_tests=8000000);atomic_json(out/(term+'-tree-check.json'),check);print(term,check['passed'],check['failures'],flush=True)
 if not check['passed']:raise ValueError('Owned vein root needs repair; no installation')
 branches.append(closed);caps.append(cap)
cava=get('LOCAL:female:vena_cava:1797');loops=boundary_loops(cava.faces);inner,field=inward_field(cava.vertices,cava.faces,.0002);outer=cava.vertices.copy();portloops=[loop for loop in loops if 1.23<cava.vertices[loop,1].mean()<1.26];assert len(portloops)==3
for loop in portloops:
 xyz=outer[loop];c=xyz.mean(0);_,_,axes=np.linalg.svd(xyz-c);n=axes[-1];rad=xyz-c;rad-=np.outer(rad@n,n);outer[loop]=c+rad*1.04;inner[loop]=c+rad*.94
wall,ls,cs=finite_wall(outer,inner,cava.faces);cv=wall['PairedPositions'].astype(float);cf=wall['FiniteFaces'].astype(int);cm=Mesh(cv,cf,frame);cm.save(out/'cava.npz');check=audit(cm,require_closed=True,max_pairs=500000,max_tests=8000000);atomic_json(out/'cava-check.json',check);print('cava',check['passed'],check['failures'],flush=True)
if not check['passed']:raise ValueError('Native caval wall failed')
ids=[]
for loop in portloops:
 j=next(j for j,ring in enumerate(ls)if set(ring)==set(loop));keys={tuple(sorted(map(int,t)))for t in cs[j]};ids.append(np.array([i for i,t in enumerate(cf)if tuple(sorted(map(int,t)))in keys]))
centers=np.array([m.vertices[np.unique(m.faces[cap])].mean(0)for m,cap in zip(branches,caps)]);targets=np.array([cv[np.unique(cf[i])].mean(0)for i in ids]);assignment=port_assignment;rows=[]
for k,(m,cap,j)in enumerate(zip(branches,caps,assignment)):
 q=attach_bridge(m.vertices,m.faces,cap,cv,cf,ids[j],stations=20,controls=32);new=Mesh(q['vertices'],q['faces'],frame);origin=np.array([0,1.25,0]);new=Mesh((new.vertices-origin).astype('<f4').astype(float)+origin,new.faces,frame);new.save(out/('branch-'+str(k)+'.npz'));check=audit(new,require_closed=True,max_pairs=500000,max_tests=8000000);atomic_json(out/('branch-'+str(k)+'-check.json'),check);print('joined',k,check['passed'],check['failures'],flush=True)
 if not check['passed']:raise ValueError('Caval collar failed')
 rows.append(dict(term=['A12.3.09.006','A12.3.09.007','A12.3.09.008'][k],branch=str(out/('branch-'+str(k)+'.npz')),branch_interface_faces=q['interface_faces'].tolist(),cava_interface_faces=ids[j].tolist(),native_port_center_m=targets[j].tolist()))
atomic_json(out/'outflow.json',dict(passed=True,branches=rows,cava=str(out/'cava.npz'),native_cava_preserved_outside_hepatic_ports=True,field=field,lumen_disks_added=0,elapsed_seconds=time.monotonic()-start))
