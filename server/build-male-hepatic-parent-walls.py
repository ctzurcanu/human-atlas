"""Complete native arterial/caval parent walls at literal hepatic mouths."""
from pathlib import Path
from itertools import combinations
import argparse,importlib.util,json,sys,time
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim
from geometry_tools.tissue import pack,normal_words
ap=argparse.ArgumentParser();ap.add_argument('--fit-root-midline',action='store_true');ap.add_argument('--reference-mouth',action='store_true');ap.add_argument('--replay-course',action='store_true');ap.add_argument('--annular-discovery',action='store_true');ap.add_argument('--role',choices=['artery','portal','hepatic_vein']);ap.add_argument('--stage');ap.add_argument('--bounded-field',action='store_true');args=ap.parse_args()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/(args.stage or ('male-parent-walls-v230'if args.annular_discovery else'male-parent-walls-v229'if args.replay_course else'male-parent-walls-v228'if args.reference_mouth else'male-parent-walls-v227'if args.fit_root_midline else'male-parent-walls-v225'));out.mkdir(exist_ok=False);start=time.monotonic()
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cat.read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame='model:local-ta98:assembled-metres';results={}
def annular_patch(network,center,normal):
 from collections import defaultdict
 from scipy.sparse import coo_matrix
 from scipy.sparse.csgraph import connected_components
 points,inv=np.unique(network.vertices,axis=0,return_inverse=True);f=inv[network.faces];tri=points[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);fn/=np.linalg.norm(fn,axis=1)[:,None]
 selected=np.flatnonzero((abs(fn@normal)>.75)&(np.linalg.norm(tri.mean(1)-center,axis=1)<.012));sf=f[selected];edges=defaultdict(list)
 for i,t in enumerate(sf):
  for x,y in zip(t,np.roll(t,-1)):edges[tuple(sorted((int(x),int(y))))].append(i)
 pairs=[p for p in edges.values()if len(p)==2];graph=coo_matrix((np.ones(len(pairs)*2),(np.r_[[p[0]for p in pairs],[p[1]for p in pairs]],np.r_[[p[1]for p in pairs],[p[0]for p in pairs]])),shape=(len(sf),len(sf))).tocsr();_,labels=connected_components(graph,directed=False);choices=[]
 for label in set(labels):
  take=np.flatnonzero(labels==label)
  if len(take)<12:continue
  try:rings=boundary_loops(sf[take])
  except ValueError:continue
  if len(rings)!=2:continue
  cc=points[np.unique(sf[take])].mean(0);order=np.argsort([np.linalg.norm(points[loop]-cc,axis=1).mean()for loop in rings]);outer,inner=rings[order[-1]],rings[order[0]]
  if len(outer)!=len(inner):continue
  _,ii=linear_sum_assignment(np.linalg.norm(points[outer,None]-points[inner][None],axis=2));choices.append((np.linalg.norm(cc-center),selected[take],points[outer],points[inner][ii]))
 if not choices:raise ValueError('No actual annular mouth found near the declared native interface')
 _,face_ids,outer,inner=min(choices,key=lambda x:x[0]);return face_ids,outer,inner
for role,pid in [('artery','ZA:Coeliac trunk'),('hepatic_vein','ZA:Inferior vena cava (abdominal part)'),('portal','TA98HEPATIC:male:portal-distal')]:
 if args.role and role!=args.role:continue
 parent=r.read_part(a,parts[pid],frame=frame,catalogue_path=cat);network=r.read_part(a,parts['TA98HEPATIC:male:'+role],frame=frame,catalogue_path=cat)
 used=np.unique(parent.faces);v,inv=np.unique(parent.vertices[used],axis=0,return_inverse=True);mapping=np.full(len(parent.vertices),-1,int);mapping[used]=inv;f=mapping[parent.faces];loops=boundary_loops(f)
 nv=np.unique(network.vertices[np.unique(network.faces)],axis=0);tree=cKDTree(nv);dist=[tree.query(v[loop])[0].max()for loop in loops];j=int(np.argmin(dist));loop=loops[j];assert dist[j]<.002
 # Recover the two actual wall faces on either side of the source mid-surface.
 root_midpoints=v[loop]
 if args.reference_mouth:
  reference_id='ZA:Common hepatic artery'if role=='artery'else'ZA:Inferior vena cava (thoracic part)'
  reference=r.read_part(a,parts[reference_id],frame=frame,catalogue_path=cat);used2=np.unique(reference.faces);rv,ri=np.unique(reference.vertices[used2],axis=0,return_inverse=True);rm=np.full(len(reference.vertices),-1,int);rm[used2]=ri;rf=rm[reference.faces];rl=boundary_loops(rf)
  ring=min(rl,key=lambda ring:np.linalg.norm(rv[ring].mean(0)-v[loop].mean(0)));assert len(ring)==len(loop)
  _,order=linear_sum_assignment(np.linalg.norm(v[loop,None]-rv[ring][None],axis=2));root_midpoints=rv[ring][order]
  if args.replay_course:
   from geometry_tools.hepatic_course_field import HepaticCourseField
   spec=importlib.util.spec_from_file_location('hepatic_builder',ROOT/'server/build-male-hepatic-connections.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
   if role=='artery':
    field=HepaticCourseField(a,b/'hepatic-contained-course-v163/field.npz',b/'separated-courses-v105/evidence.json');root_midpoints=field(module.warp(root_midpoints,role,b))
 discovered=None
 if args.annular_discovery:
  _,_,axis=np.linalg.svd(v[loop]-v[loop].mean(0));face_ids,op,ip=annular_patch(network,v[loop].mean(0),axis[-1]);assert len(op)==len(loop)
  _,order=linear_sum_assignment(np.linalg.norm(v[loop,None]-op[None],axis=2));discovered=(face_ids,op[order],ip[order]);pairs=np.stack([op[order],ip[order]],axis=1);residual=np.linalg.norm(op[order]-v[loop],axis=1).tolist()
 else:pairs=[];residual=[]
 for point in ([]if discovered else root_midpoints):
  dd,ii=tree.query(point,k=12);choices=[]
  for x,y in combinations(range(12),2):
   if max(dd[x],dd[y])>.003 or min(dd[x],dd[y])<1e-8:continue
   err=np.linalg.norm((nv[ii[x]]+nv[ii[y]])*.5-point);choices.append((err,x,y))
  err,x,y=min(choices);pairs.append([nv[ii[x]],nv[ii[y]]]);residual.append(err)
 # A native source mismatch needs an explicit bounded root fit rather than
 # being declared a literal source correspondence. Only the end controls move.
 limit=.012 if discovered else .0005 if args.fit_root_midline else 2e-7
 assert max(residual)<limit,('Unmatched native centreline mouth',role,max(residual))
 sums=np.zeros_like(v)
 for i,oldid in enumerate(used):sums[inv[i]]+=parent.attributes['normals'][oldid]
 normal=sums/np.linalg.norm(sums,axis=1)[:,None];pairs=np.asarray(pairs);thickness=np.median(np.linalg.norm(pairs[:,0]-pairs[:,1],axis=1))*.5
 fields=np.array([v.copy(),v-2*thickness*normal])if discovered else np.array([v+thickness*normal,v-thickness*normal])
 if args.bounded_field:
  from geometry_tools.algorithms.finite_normal_field import inward_field
  fields[1],field_evidence=inward_field(v,f,2*thickness);atomic_json(out/(role+'-field.json'),field_evidence)
 if args.replay_course:
  delta=root_midpoints.mean(0)-v[loop].mean(0);center=v[loop].mean(0);length=min(np.linalg.norm(v[ring].mean(0)-center)for k,ring in enumerate(loops)if k!=j);t=np.clip(np.linalg.norm(v-center,axis=1)/length,0,1);weight=1-t*t*(3-2*t);fields+=weight[None,:,None]*delta
 for i,vid in enumerate(loop):
  order=np.argsort((pairs[i]-v[vid])@normal[vid]);fields[0,vid]=pairs[i,order[-1]];fields[1,vid]=pairs[i,order[0]]
 n=len(v);vertices=fields.reshape(-1,3);root_lookup={tuple(point):i for i,point in enumerate(vertices)};root_set=set(map(tuple,pairs.reshape(-1,3)));root_faces=[];root_ids=[]
 effective=network.attributes['normals'].astype(np.int64)*network.attributes['normalSigns'][:,None];root_words=[]
 for fi,tri in enumerate(network.faces):
  xyz=network.vertices[tri]
  if all(tuple(p)in root_set for p in xyz):
   root_faces.append([root_lookup[tuple(p)]for p in xyz[[0,2,1]]]);root_words.append(-effective[tri[[0,2,1]]]);root_ids.append(fi)
 assert len(root_faces)==2*len(loop),('Expected an annular tissue interface',role,len(root_faces),len(loop))
 surfaces={'outer':f,'inner':f[:,[0,2,1]]+n,'hepatic-mouth':np.asarray(root_faces,int)};cells={'parent':[dict(surface=sid,side=1,neighbor=None)for sid in surfaces]}
 for k,ring in enumerate(loops):
  if k==j:continue
  sid='native-mouth:'+str(k);surfaces[sid]=annular_rim(ring,ring+n);cells['parent'].append(dict(surface=sid,side=1,neighbor=None))
 words=normal_words((vertices-[0,1.25,0]).astype('<f4').astype(float),surfaces);words['hepatic-mouth']=np.asarray(root_words,'<i2')
 result=pack(vertices,surfaces,cells,{'parent':pid},frame,out/role,origin=[0,1.25,0],corner_words=words,max_pairs=500000,max_tests=8000000,checkpoint_directory=out/(role+'-checks'))
 result['parts']['parent']['tissueDrawGroups'][2]['neighbor']='TA98HEPATIC:male:'+role
 results[role]=dict(passed=result['passed'],parent=pid,bank=result['bank'],bank_sha256=result['bank_sha256'],part=result['parts']['parent'],network_interface_faces=root_ids,parent_interface_faces=list(range(len(f)*2,len(f)*2+len(root_faces))),maximum_mid_surface_residual_m=max(residual),half_wall_thickness_m=thickness)
 atomic_json(out/(role+'-result.json'),results[role]);print(role,result['passed'],result['checks']['failures'],flush=True)
 if not result['passed']:raise ValueError('Native parent wall failed; not installed')
atomic_json(out/'parents.json',dict(passed=True,interfaces=results,elapsed_seconds=time.monotonic()-start))
