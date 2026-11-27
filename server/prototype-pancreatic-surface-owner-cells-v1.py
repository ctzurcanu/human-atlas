#!/usr/bin/env python3
"""One source-fixed five-owner annular seam prototype; no atlas installation."""
import json,time,sys
from pathlib import Path
from itertools import combinations
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_rim_bridge import bridge_boundaries
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.tissue import validate_cells

def known_plane_annulus(v, a, b, normal):
 import manifold3d
 from collections import Counter
 from geometry_tools.algorithms.tube_rim_bridge import _edge_chain
 n=np.asarray(normal,float);n/=np.linalg.norm(n);origin=v[a].mean(0)
 loops=[np.asarray(a)[::-1],np.asarray(b)[::-1]]
 def project(n):
  u=np.cross(np.eye(3)[np.argmin(abs(n))],n);u/=np.linalg.norm(u);w=np.cross(n,u)
  return [(v[q]-origin)@np.array([u,w]).T for q in loops]
 area=lambda q:float((q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1)).sum()/2)
 xy=project(n)
 if area(xy[0])<0:n=-n;xy=project(n)
 if not area(xy[0])>0 or not area(xy[1])<0:raise ValueError('Known-plane annulus needs an opposed inner hole')
 residual=float(abs((v[np.r_[a,b]]-origin)@n).max())
 if residual>2e-12:raise ValueError('Annular controls are outside the declared common cut plane')
 index=np.r_[loops[0],loops[1]];tri=np.asarray(manifold3d.triangulate(xy,epsilon=1e-12,allow_convex=False),int);faces=index[tri]
 expected=Counter((int(y),int(x))for q in [a,b]for x,y in zip(q,np.roll(q,-1)))
 if _edge_chain(faces)!=expected:raise ValueError('Known-plane triangulation changed a literal boundary chain')
 points=v[faces];actual=float(np.linalg.norm(np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0]),axis=1).sum()/2);target=area(xy[0])+area(xy[1])
 if abs(actual-target)>max(1e-16,target*1e-9):raise ValueError('Known-plane annular area mismatch')
 return faces,dict(method='declared-plane polygon with one inner hole',plane_residual_m=residual,area_m2=actual,source_controls_literal=True,lumen_disk_created=False)

def main():
 start=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion/pancreatic-public-cohort-v1/finite-network-construction-v1';plan=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else B/'surface-owner-seam-plan-v1.json';p=json.loads(plan.read_text());source=B/p['source'];z=np.load(source);v=z['vertices'];f=z['faces'];owners=z['owners'];diagnostic=Path(p.get('diagnostic',str(B/'surface-owner-diagnostic-v1/proof.json')));d=json.loads(diagnostic.read_text());assert d['passed'];O=Path(p.get('output',str(B/'surface-owner-cells-prototype-v1')));O.mkdir(exist_ok=False);status={'schema':'pancreatic-five-owner-surface-cells/v1','plan_sha256':digest_file(plan),'source_sha256':digest_file(source),'diagnostic_sha256':digest_file(diagnostic),'passed':False,'charts':{},'anatomical_acceptance':False,'installed':False};checkpoint=O/'checkpoint.json'
 def save():
  status['elapsed_seconds']=time.monotonic()-start;atomic_json(checkpoint,status,overwrite=checkpoint.exists())
 pairs=[tuple(sorted(x))for x in d['expected_owner_pairs']];loops=d['directed_seam_loops'];charts={};cells={name:[]for name in p['expected_genus']};surfaces={};uses={name:[]for name in cells}
 for name in cells:
  sid='free:'+name;surfaces[sid]=f[owners==name];cells[name].append(surfaces[sid]);uses[name].append({'surface':sid,'side':1});assert len(surfaces[sid])
 # Existing DP enumerates the second ring's phases. Ring order does not change
 # the required opposed boundary chain; put the smaller ring second.
 for pair in sorted(pairs,key=lambda x:np.prod([len(loops['|'.join(x)+'|'+s]['cycles'][0])for s in ['outer','inner']])*min(len(loops['|'.join(x)+'|outer']['cycles'][0]),len(loops['|'.join(x)+'|inner']['cycles'][0]))):
  key='|'.join(pair);a=np.array(loops[key+'|outer']['cycles'][0]);b=np.array(loops[key+'|inner']['cycles'][0]);first,second=(a,b)if len(b)<=len(a)else(b,a);status['current_stage']='annulus:'+key;save()
  reuse=Path(p['reused_annuli_directory'])/(key.replace('|','-')+'-annulus.npz')if p.get('reused_annuli_directory')else None
  if reuse is not None and reuse.exists():
   old=np.load(reuse);seam=old['faces'];assert np.array_equal(old['vertices'],v)and np.array_equal(old['outer_boundary'],a)and np.array_equal(old['inner_boundary'],b)
   e=dict(method='literal completed source annulus reuse',source=str(reuse),sha256=digest_file(reuse),source_controls_literal=True,lumen_disk_created=False)
  elif pair==('dorsal','inferior')and p.get('known_plane_normal')is not None:seam,e=known_plane_annulus(v,a,b,p['known_plane_normal'])
  else:seam,e=bridge_boundaries(v,first,second)
  sid='seam:'+key;charts[sid]=seam;surfaces[sid]=seam;cells[pair[0]].append(seam);cells[pair[1]].append(seam[:,::-1]);uses[pair[0]].append({'surface':sid,'side':1,'neighbor':pair[1]});uses[pair[1]].append({'surface':sid,'side':-1,'neighbor':pair[0]});status['charts'][sid]={'pair':pair,'proof':e,'source_outer_controls':a,'source_inner_controls':b};atomic_npz(O/(key.replace('|','-')+'-annulus.npz'),vertices=v,faces=seam,outer_boundary=a,inner_boundary=b);save()
 meshes={name:Mesh(v,np.concatenate(ff),str(z['frame']))for name,ff in cells.items()};arrays={'vertices':v,'original_external_faces':f,'original_external_source_face_ancestry':z['whole_source_face_ancestry'],'original_external_source_original_ids':z['source_original_ids'],'original_external_source_physical_face_ids':z['source_physical_face_ids']}
 for name,m in meshes.items():arrays['cell_'+name+'_faces']=m.faces;m.save(O/(name+'-cell.npz'))
 for i,(sid,ff)in enumerate(surfaces.items()):arrays['surface_'+str(i)]=ff
 atomic_npz(O/'shared-arrangement-before-gates.npz',**arrays);status['cell_uses']=uses;status['surface_arrays']={sid:'surface_'+str(i)for i,sid in enumerate(surfaces)};status['all_source_external_faces_literal_and_once']=bool(sum(len(ff)for sid,ff in surfaces.items()if sid.startswith('free:'))==len(f));status['canonical_reciprocal_chart_identity']=True;status['source_positions_prefix_literal']=True;status['current_stage']='all5cells saved; numerical gates';save()
 audits={}
 for name,m in meshes.items():
  audits[name]=audit(m,require_closed=True,expected_components=1,expected_genus=p['expected_genus'][name]);status['audits']=audits;status['current_stage']='audited:'+name;save()
  if not audits[name]['passed']:status['reason']='cell-audit:'+name;atomic_json(O/'proof.json',status);return
 quality=validate_cells(meshes);status['finite_cell_quality']=quality;save();origin=v.mean(0);solids={name:make_solid(m.vertices-origin,m.faces)[0]for name,m in meshes.items()};joined=None
 for solid in solids.values():joined=solid if joined is None else joined+solid
 original=make_solid(v-origin,f)[0];extra=joined-original;missing=original-joined;status['union_replay']={'extra_empty':bool(extra.is_empty()),'missing_empty':bool(missing.is_empty()),'extra_volume_m3':float(extra.volume()),'missing_volume_m3':float(missing.volume()),'original_volume_m3':float(original.volume()),'union_volume_m3':float(joined.volume())};status['passed']=bool(quality['passed']and extra.is_empty()and missing.is_empty()and status['all_source_external_faces_literal_and_once']);status['current_stage']='complete';save();atomic_json(O/'proof.json',status);print('PASS',status['passed'],'audits',{n:q['failures']for n,q in audits.items()},'quality',quality['failures'],'union',status['union_replay'],flush=True)
if __name__=='__main__':main()
