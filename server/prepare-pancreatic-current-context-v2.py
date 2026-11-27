#!/usr/bin/env python3
"""Freeze a new source-bound context scope only after network geometry passes."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.finite_tissue_boolean import make_solid

def main():
 p=argparse.ArgumentParser();p.add_argument('--network-dir',required=True);p.add_argument('--output',required=True);args=p.parse_args();N=Path(args.network_dir).resolve();D=Path(args.output).resolve();proofpath=N/'proof.json';proof=json.loads(proofpath.read_text());assert proof['passed']is True,'Network proof has not passed';wallpath=N/'wall.npz';z=np.load(wallpath);v=z['vertices'];f=z['faces'];tags=z['source_original_ids'];faces=z['source_face_ids'];names=json.loads(str(z['source_id_map_json']));inputs=proof['inputs'];native_names=['nativeouter','nativefluid'];native_tags={name:int(next(k for k,val in names.items()if val==name))for name in native_names}
 for name,tag in native_tags.items():assert inputs[name]['manifold_original_id']==tag
 D.mkdir(exist_ok=False);mask=~np.isin(tags,list(native_tags.values()));footprints={}
 for name,tag in native_tags.items():
  row=inputs[name];path=Path(row['path']);assert digest_file(path)==row['sha256'];q=np.load(path);ov=q[row['vertices_field']];of=q['faces'];tri=ov[of];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2;take=np.flatnonzero(tags==tag);emitted=v[f[take]];ar=np.linalg.norm(np.cross(emitted[:,1]-emitted[:,0],emitted[:,2]-emitted[:,0]),axis=1)/2;ret=np.bincount(faces[take],weights=ar,minlength=len(of));partial=np.flatnonzero((ret>0)&(abs(ret-area)>np.maximum(1e-16,area*1e-8)));mask[take[np.isin(faces[take],partial)]]=True;footprints[name]={'source_input_face_ids':partial,'emitted_fragments_selected':int(np.sum(np.isin(faces[take],partial))),'source_sha256':row['sha256'],'qualifier':'Partial emitted source facets only; complete source-face losses need independent root-footprint replay.'}
 subset=Mesh(v,f[mask],'model:local-ta98:assembled-metres');low,high=subset.bounds;A=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';atlas=json.loads(A.read_text());V=ROOT/'.local-models/geometry-workbench/circulation-census-v1/aortic-blood-installed-v1/resolved-male.json';view=json.loads(V.read_text());assert digest_file(A)==view['sourceSha256'];assert view['view']['mode']=='ordinary-local-atlas';visibility={q['id']:q for q in view['effectivePartVisibility']};parts={q['id']:q for q in atlas['parts']};owners=[]
 force=['ZA:Abdominal aorta','ZA:Coeliac trunk','ZA:Superior mesenteric artery','ZA:Pancreatic duct','ZA:Accessory pancreatic duct','TA98PART:FJ3079:A05.8.02.013','TA98HEPATIC:male:portal','TA98HEPATIC:male:biliary','ZA:Splenic vein','TA98PART:ZA:Splenic artery:A12.2.12.051'];replaced='TA98PART:ZA:Splenic artery:trunk'
 for key,part in parts.items():
  if not visibility[key]['available']or key==replaced:continue
  bounds=np.array(part.get('bounds',[[np.inf]*3,[-np.inf]*3]));overlap=not(np.any(bounds[0]>high)or np.any(bounds[1]<low))
  if not overlap and key not in force:continue
  name=(key+' '+part.get('name','')).lower();cat='other-active'
  if key.startswith('TA98PART:ZA:Pancreas:'):cat='allowed-pancreatic-parenchyma'
  elif key=='TA98PART:ZA:Splenic artery:A12.2.12.051':cat='declared-hilar-interface-peer'
  elif 'duct'in name or key=='TA98HEPATIC:male:biliary':cat='duct-critical'
  elif part.get('system')=='venous':cat='vein-critical'
  elif part.get('system')=='skeletal'or'bone'in part.get('material','').lower():cat='skeleton-critical'
  elif any(w in name for w in ['stomach','duodenum','jejun','ileum','colon','rectum','oesoph','esoph']):cat='gastrointestinal-critical'
  elif any(w in name for w in ['liver','spleen','diaphragm']):cat='organ-critical'
  elif part.get('system')=='arterial':cat='artery-critical'
  stage='critical-duct-vein'if cat in ['duct-critical','vein-critical','declared-hilar-interface-peer']else('skeleton-organs'if cat in ['skeleton-critical','organ-critical']else('artery-gi'if cat in ['artery-critical','gastrointestinal-critical']else('allowed-pancreas'if cat=='allowed-pancreatic-parenchyma'else'other-fixed')))
  owners.append({'id':key,'name':part.get('name'),'category':cat,'stage':stage,'defaultHidden':visibility[key]['defaultHidden'],'record_bounds_metres':bounds,'forced_critical_even_when_AABB_disjoint':key in force and not overlap})
 for key in force:assert any(q['id']==key for q in owners),'Mandatory target unavailable: '+key
 atomic_npz(D/'new-wall-face-scope.npz',vertices=v,faces=f[mask],whole_wall_face_ids=np.flatnonzero(mask),source_original_ids=tags[mask],source_face_ids=faces[mask],frame=np.array(subset.frame));scope={'schema':'pancreatic-new-wall-current-context-scope/v2','catalogue_sha256':digest_file(A),'visibility_snapshot_sha256':digest_file(V),'wall_npz_sha256':digest_file(wallpath),'wall_proof_sha256':digest_file(proofpath),'network_directory':str(N),'whole_wall_faces':len(f),'selected_new_and_root_footprint_faces':int(mask.sum()),'excluded_unchanged_native_parent_faces':int((~mask).sum()),'source_id_names':names,'root_footprint_ancestry':footprints,'scope_bounds_metres':[low,high],'owners':owners,'replaced_native_parent_excluded':replaced,'foreign_artery_exemptions':[],'allowed_parenchyma_rule':'Only actual pancreas finite bulk cells; no pancreatic-node/duct exception.','anatomical_acceptance':False};atomic_json(D/'scope-v2-before-tests.json',scope)
 whole=Mesh.load(N/'outer.npz',frame=subset.frame);native_row=inputs['nativeouter'];q=np.load(native_row['path']);native=Mesh(q[native_row['vertices_field']],q['faces'],subset.frame);origin=np.r_[whole.vertices,native.vertices].mean(0);ws,_=make_solid(whole.vertices-origin,whole.faces);ns,_=make_solid(native.vertices-origin,native.faces);new=ws-ns;out=new.to_mesh64();vv=np.asarray(out.vert_properties)[:,:3]+origin;ff=np.asarray(out.tri_verts);atomic_npz(D/'new-occupied-envelope.npz',vertices=vv,faces=ff,frame=np.array(subset.frame));atomic_json(D/'new-occupied-envelope-proof.json',{'schema':'pancreatic-new-occupied-comparison-envelope/v2','whole_outer_sha256':digest_file(N/'outer.npz'),'native_outer_sha256':native_row['sha256'],'operation':'passed candidate outer minus unchanged native outer; comparison volume only','new_volume_m3':float(new.volume()),'status':str(new.status()),'temporary_comparison_only':True,'source_or_atlas_modified':False});print('PREPARED',len(owners),Counter(q['stage']for q in owners),'newfaces',int(mask.sum()),flush=True)
if __name__=='__main__':main()
