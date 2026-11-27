"""Adapt the female's native recipients to the fitted liver's actual mouths."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.annular_bridge import attach_bridge,find_annular_patch,trim_named_inlet
from geometry_tools.algorithms.finite_vessel_wall import finite_wall
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank,split_closed_cell
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
from geometry_tools.algorithms.finite_normal_field import inward_field
ap=argparse.ArgumentParser();ap.add_argument('--role',choices=['artery','portal','biliary'],required=True);ap.add_argument('--native-course',action='store_true');ap.add_argument('--local-wall-storage',action='store_true');ap.add_argument('--bounded-field',action='store_true');args=ap.parse_args();role=args.role
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/(('female-recipient-v237-'if args.bounded_field else'female-recipient-v234-'if args.local_wall_storage else'female-recipient-v232-'if args.native_course else'female-recipient-v231-')+role);out.mkdir(exist_ok=False);start=time.monotonic()
male_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';female_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(female_path.read_text());male=json.loads(male_path.read_text());parts={p['id']:p for p in a['parts']};mp={p['id']:p for p in male['parts']};r=AtlasReader();frame='construction:female-liver:unscaled-donor-metres'
fit=json.loads((b/'final-context-v217/fit.json').read_text())['chosen'];old=a['liverTransferIntegration'];scale=old['scale']*fit['additional_scale'];translation=np.array(old['translation_m'])*fit['additional_scale']+(1-fit['additional_scale'])*np.array(fit['center_m'])+fit['translation_m']
def compact(m):
 used=np.unique(m.faces);v,inv=np.unique(m.vertices[used],axis=0,return_inverse=True);mapping=np.full(len(m.vertices),-1,int);mapping[used]=inv;return Mesh(v,mapping[m.faces],frame)
def native(pid):
 m=r.read_part(a,parts[pid],frame='model:local-female-ta98:assembled-metres',catalogue_path=female_path);return compact(Mesh((m.vertices-translation)/scale,m.faces,frame))
def child_mesh():
 if role=='biliary':return compact(Mesh.load(b/('common-hepatic-mouth-v235/child.npz'if args.bounded_field else'right-caudate-duct-v226/emitted-wall.npz')))
 return compact(r.read_part(male,mp['TA98HEPATIC:male:'+role],frame=frame,catalogue_path=male_path))
def clip_open(m,limit):
 bank=SharedCutBank(m.vertices);result=[];source=[]
 for fi,t in enumerate(m.faces):
  poly=[];d=m.vertices[t,1]-limit
  for k in range(3):
   x,y=int(t[k]),int(t[(k+1)%3]);dx,dy=d[k],d[(k+1)%3]
   if dx<=0:poly.append(x)
   if dx*dy<0:poly.append(bank.intersection(x,y,np.array([0.,1.,0.]),limit,'native-portal-upper'))
  if poly:
   j=poly.index(min(poly));poly=poly[j:]+poly[:j]
  for k in range(1,len(poly)-1):result.append([poly[0],poly[k],poly[k+1]]);source.append(fi)
 clipped=compact(Mesh(np.array(bank.positions),np.asarray(result,int),frame))
 if args.native_course:
  from scipy.sparse import coo_matrix
  from scipy.sparse.csgraph import connected_components
  loops=boundary_loops(clipped.faces);target=(np.array([-.010,1.185,.030])-translation)/scale
  choices=[loop for loop in loops if abs(clipped.vertices[loop,1].mean()-limit)<1e-8];ring=min(choices,key=lambda loop:np.linalg.norm(clipped.vertices[loop].mean(0)-target))
  f=clipped.faces;edges=np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]];graph=coo_matrix((np.ones(len(edges)*2),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(clipped.vertices),len(clipped.vertices))).tocsr();_,labels=connected_components(graph,directed=False);take=labels[f[:,0]]==labels[ring[0]];clipped=compact(Mesh(clipped.vertices,f[take],frame))
 return clipped
def finite_parent(m,thickness,choose):
 tri=m.triangles;fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal=np.zeros_like(m.vertices)
 for j in range(3):np.add.at(normal,m.faces[:,j],fn)
 length=np.linalg.norm(normal,axis=1);assert np.all(length>0);normal/=length[:,None]
 origin=np.array([0,1.25,0])if args.local_wall_storage else np.zeros(3)
 if args.bounded_field:
  inner,field=inward_field(m.vertices,m.faces,thickness);atomic_json(out/'wall-field.json',field)
 else:inner=m.vertices-thickness*normal
 wall,loops,caps=finite_wall(m.vertices-origin,inner-origin,m.faces);pv=wall['PairedPositions'].astype(float)+origin;pf=wall['FiniteFaces'].astype(int);j=choose(loops,m.vertices)
 capkeys={tuple(sorted(map(int,t)))for t in caps[j]};ids=np.array([i for i,t in enumerate(pf)if tuple(sorted(map(int,t)))in capkeys]);assert len(ids)
 return Mesh(pv,pf,frame),ids
child=child_mesh()
if role=='artery':
 pid='LOCAL:female:proper_hepatic_artery:1718'if args.native_course else'LOCAL:female:common_hepatic_artery:1722';source=native(pid)
 if args.bounded_field:source=clip_open(source,(1.1875-translation[1])/scale)
 parent,pc=finite_parent(source,.00015/scale,lambda loops,v:int(np.argmax([v[loop,1].mean()for loop in loops]))if args.native_course else int(np.argmax([v[loop,2].mean()for loop in loops])))
 cc=json.loads((b/'male-parent-walls-v230/artery-result.json').read_text())['network_interface_faces']
 if args.native_course:
  from scipy.sparse import coo_matrix
  from scipy.sparse.csgraph import connected_components
  bank=SharedCutBank(child.vertices);cut=split_closed_cell(bank,child.faces,[0,1,0],1.25390625,'female-native-proper',allow_holes=True);v=np.asarray(bank.positions);nf=cut.negative_faces
  edges=np.r_[nf[:,[0,1]],nf[:,[1,2]],nf[:,[2,0]]];graph=coo_matrix((np.ones(len(edges)*2),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),len(v))).tocsr();_,labels=connected_components(graph,directed=False)
  old_root=set(cc);where=np.flatnonzero(np.isin(cut.negative_ancestry,list(old_root)));assert len(where);discard=labels[nf[where[0],0]];drop=labels[nf[:,0]]==discard
  removed_keys={tuple(sorted(map(int,t)))for t in nf[drop]};retained=np.r_[cut.positive_faces,nf[~drop]];owners={}
  for i,t in enumerate(retained):owners.setdefault(tuple(sorted(map(int,t))),[]).append(i)
  keep=np.array([len(owners[tuple(sorted(map(int,t)))])==1 for t in retained]);final=retained[keep];capkeys={tuple(sorted(map(int,t)))for t in cut.interface_faces if tuple(sorted(map(int,t)))in removed_keys};cc=np.array([i for i,t in enumerate(final)if tuple(sorted(map(int,t)))in capkeys]);assert len(cc)
  child=compact(Mesh(v,final,frame));child.save(out/'intrahepatic-artery.npz');atomic_json(out/'arterial-native-course.json',dict(native_common_artery_preserved=True,native_proper_artery_extended=True,donor_extrahepatic_stem_removed=True,other_caudad_branches_preserved=True,cut_y_donor_m=1.25390625))
elif role=='portal':
 pid='LOCAL:female:portal_vein:1795';limit=(1.185-translation[1])/scale;source=clip_open(native(pid),limit)
 parent,pc=finite_parent(source,.00025/scale,lambda loops,v:int(np.argmin([abs(v[loop,1].mean()-limit)for loop in loops])))
 cc=find_annular_patch(child.vertices,child.faces,[-.0058649,1.242,.032422],[0,1,0])
else:
 pid='LOCAL:female:bile_duct:1031';source=native(pid)
 # Cut above the native cystic confluence. The preserved wall below this
 # annular cut carries the existing pancreatic/duodenal connections.
 cut_y=1.2075;limit=(cut_y-translation[1])/scale;bank=SharedCutBank(source.vertices);split=split_closed_cell(bank,source.faces,[0,1,0],limit,'female-common-hepatic',allow_holes=True)
 parent=Mesh(np.asarray(bank.positions),split.negative_faces,frame);pc=np.flatnonzero(split.negative_ancestry<0)
 if args.native_course:
  lower=r.read_part(male,mp['TA98PART:FJ3079:A05.8.01.061'],frame=frame,catalogue_path=male_path);row=mp['TA98PART:FJ3079:A05.8.01.061'];indices=np.r_[tuple(np.arange(g['start']//3,(g['start']+g['count'])//3)for g in row['tissueDrawGroups']if g.get('neighbor')=='TA98HEPATIC:male:biliary')];keys={tuple(sorted(map(tuple,tri)))for tri in lower.triangles[indices]};cc=np.array([i for i,tri in enumerate(child.triangles)if tuple(sorted(map(tuple,tri)))in keys]);assert len(cc)==len(indices),('Source-bound biliary mouth changed',len(cc),len(indices))
 else:cc=find_annular_patch(child.vertices,child.faces,[-.024,1.2518,.038],[0,1,0])
 # Carry all shared native cystic/duct controls through the same positive
 # monotone Y map. Actual gallbladder attachment coordinates remain fixed.
 cystic=r.read_part(a,parts['LOCAL:female:cystic_duct:authored'],frame='model:local-female-ta98:assembled-metres',catalogue_path=female_path);source_row=parts['LOCAL:female:cystic_duct:authored'];gall_faces=np.r_[tuple(np.arange(g['start']//3,(g['start']+g['count'])//3)for g in source_row['tissueDrawGroups']if g.get('neighbor','').startswith('LOCAL:female:GALLWALL:'))];freeze=float(cystic.triangles[gall_faces,:,1].max());target=float(child.vertices[np.unique(child.faces[np.asarray(cc)]),1].mean()*scale+translation[1]-.002);shift=cut_y-target
 assert 0<shift<(cut_y-freeze)/1.5,'Monotone duct map would fold'
 def warp(m):
  w=m.vertices*scale+translation;t=np.clip((w[:,1]-freeze)/(cut_y-freeze),0,1);w[:,1]-=shift*t*t*(3-2*t);return Mesh((w-translation)/scale,m.faces,frame)
 parent=warp(parent)
 for name in ['LOCAL:female:cystic_duct:authored','CYSTICFOLD:female:A05.8.02.012']:
  warped=warp(native(name));warped.save(out/(name.replace(':','_')+'.npz'))
 atomic_json(out/'native-duct-map.json',dict(freeze_y_m=freeze,cut_y_m=cut_y,shift_y_m=shift,minimum_jacobian=1-1.5*shift/(cut_y-freeze),gallbladder_attachments_literal=True))
source.save(out/'native-source.npz');parent.save(out/'native-parent-before-bridge.npz');child.save(out/'child.npz')
q=attach_bridge(parent.vertices,parent.faces,pc,child.vertices,child.faces,cc)
origin=np.array([0,1.25,0]);f,_,conform=conform_quantized_surfaces(q['vertices'],{'parent':q['faces']},origin);emitted=(q['vertices']-origin).astype('<f4').astype(float)+origin;new=Mesh(emitted,f['parent'],frame);new.save(out/'parent.npz')
check=audit(new,require_closed=True,max_pairs=1000000,max_tests=12000000);atomic_json(out/'parent-check.json',check)
atomic_json(out/'bridge.json',dict(passed=check['passed'],parent_id=pid,role=role,parent=str(out/'parent.npz'),child=str(out/'child.npz'),parent_interface_faces=q['interface_faces'].tolist(),child_interface_faces=list(map(int,cc)),source_controls_literal=q['parent_source_controls_literal']and q['child_source_controls_literal'],curve=q['path'],cycles=q['source_cycle_controls'],conform=conform,scale=scale,translation_m=translation.tolist(),lumen_disks_added=0,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(role=role,passed=check['passed'],failures=check['failures'],elapsed_seconds=time.monotonic()-start)),flush=True)
if not check['passed']:raise ValueError('Female recipient failed; no installation')
