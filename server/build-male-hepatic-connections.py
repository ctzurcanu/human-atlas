"""Prepare the local male hepatic connections; installation remains a reviewed job.

Finite communicating network walls and owned, optional TA98 surface regions
are kept distinct. Existing liver tissue is deliberately outside this revision.
"""
from pathlib import Path
import argparse,copy,json,sys,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_bytes,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
from geometry_tools.spatial import BVH
from geometry_tools.distance import barycentric
from geometry_tools.tissue import normal_words,decode_packed
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells

ROLES=['artery','portal','hepatic_vein','biliary']
COMMON='TA98PART:FJ3079:A05.8.01.061'
PORTAL_DISTAL='TA98HEPATIC:male:portal-distal'
NAMES={'artery':'Hepatic arterial tree','portal':'Hepatic portal vein · liver branches and junction','hepatic_vein':'Hepatic veins and caval junction','biliary':'Hepatic ducts and common hepatic junction'}
SYSTEMS={'artery':'arterial','portal':'venous','hepatic_vein':'venous','biliary':'digestive'}
PROTOTYPES={'artery':'ZA:Common hepatic artery','portal':'ZA:Hepatic portal vein','hepatic_vein':'ZA:Hepatic veins','biliary':COMMON}
LIMITS='Local macroscopic communicating network. Named branch surface regions are selectable; regional boundaries remain provisional. Liver parenchymal section channels and microscopic anatomy are not certified by this revision.'

def flat_fallback(v,faces,words,origin):
 q=(v-origin).astype('<f4').astype(float)[faces];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);bad=np.einsum('tki,ti->tk',words.astype(float),n)<=0
 if np.any(bad):
  flat=np.rint(n/np.linalg.norm(n,axis=1)[:,None]*32767).astype('<i2');words[bad]=np.broadcast_to(flat[:,None],words.shape)[bad]
 return words

def extras(b,out,atlas,cat):
 reader=AtlasReader();byid={p['id']:p for p in atlas['parts']};names=['common_lower','portal_lower'];paths=[b/'closed-biliary-fork-v87/common-hepatic-fixed-lower.npz',b/'explicit-portal-connection-v86/portal-fixed-distal.npz'];oldids=[COMMON,'ZA:Hepatic portal vein'];partids={'common_lower':COMMON,'portal_lower':PORTAL_DISTAL};vlist=[];surfaces={};cells={};words={};offset=0;group_neighbors={};origin=np.array([0,1.25,0]);source_refs={}
 for name,path,oldid in zip(names,paths,oldids):
  m=Mesh.load(path);oldpart=byid[oldid];old=reader.read_part(atlas,oldpart,frame=m.frame,catalogue_path=cat);source_refs[name]=oldid;local=normal_words((m.vertices-origin).astype('<f4').astype(float),{name:m.faces})[name];labels=np.zeros(len(m.faces),int);neighbors={};oldlabels=np.zeros(len(old.faces),int)
  for i,g in enumerate(oldpart.get('tissueDrawGroups',[])):
   oldlabels[g['start']//3:(g['start']+g['count'])//3]=i;neighbors[i]=g.get('neighbor')
  lookup={tuple(sorted(map(tuple,old.vertices[f]))):i for i,f in enumerate(old.faces)};bv=BVH(old);effective=old.attributes['normals'].astype(np.int64)
  if 'normalSigns'in old.attributes:effective=effective*old.attributes['normalSigns'][:,None]
  n=len(m.vertices)//2
  for i,f in enumerate(m.faces):
   p=m.vertices[f];key=tuple(sorted(map(tuple,p)));sourceface=lookup.get(key);root=name=='common_lower' and np.any(f<n) and np.any(f>=n) and np.all(abs(m.vertices[f[f<n],1]-1.2518)<1e-10)
   if root:labels[i]=999;neighbors[999]='TA98HEPATIC:male:biliary';continue
   if sourceface is None:sourceface=bv.closest_point(p.mean(0))['face_id']
   labels[i]=int(oldlabels[sourceface]);indices=old.faces[sourceface];t=old.vertices[indices];weights=barycentric(p,np.repeat(t[None],3,axis=0))[0];normal=weights@effective[indices].astype(float);normal/=np.linalg.norm(normal,axis=1)[:,None];local[i]=np.rint(normal*32767).astype('<i2')
   for corner,point in enumerate(p):
    dist=np.linalg.norm(t-point,axis=1);j=int(dist.argmin())
    if dist[j]<1e-12:local[i,corner]=effective[indices[j]].astype('<i2')
  if name=='common_lower':
   native_labels={label for label,peer in neighbors.items()if peer in ['FJ3080','TA98PART:FJ3079:A05.8.02.013']};keep=~np.isin(labels,list(native_labels));point_ids={tuple(point):i for i,point in enumerate(m.vertices)};replacement_faces=[];replacement_words=[];replacement_labels=[]
   for label in sorted(native_labels):
    original_faces=np.flatnonzero(oldlabels==label)
    for fi in original_faces:
     f=old.faces[fi];replacement_faces.append([point_ids[tuple(point)]for point in old.vertices[f]]);replacement_words.append(effective[f].astype('<i2'));replacement_labels.append(label)
   m=Mesh(m.vertices,np.r_[m.faces[keep],np.asarray(replacement_faces,int)],m.frame);local=np.r_[local[keep],np.asarray(replacement_words,dtype='<i2')];labels=np.r_[labels[keep],np.asarray(replacement_labels,int)]
  local=flat_fallback(m.vertices,m.faces,local,origin);cells[name]=[]
  for label in sorted(set(labels)):
   sid=name+':'+str(label);take=labels==label;surfaces[sid]=m.faces[take]+offset;words[sid]=local[take];cells[name].append(dict(surface=sid,side=1,neighbor=None));group_neighbors[sid]=neighbors.get(label)
  vlist.append(m.vertices);offset+=len(m.vertices)
 v=np.concatenate(vlist);packed=pack_finite_cells(v,surfaces,words,cells,partids,origin);target=out/'native-parents.bin';atomic_bytes(target,packed.data);meshes=decode_packed(packed,'model:local-ta98:assembled-metres',origin);checks={}
 for name,m in meshes.items():
  m.save(out/(name+'.npz'));checks[name]=audit(m,require_closed=name=='common_lower',max_pairs=1000000,max_tests=8000000)
  for g,use in zip(packed.cells[name]['tissueDrawGroups'],cells[name]):
   if group_neighbors[use['surface']]:g['neighbor']=group_neighbors[use['surface']]
 assert all(q['passed']for q in checks.values()),{n:q['failures']for n,q in checks.items()}
 meta=dict(bank=str(target),bank_sha256=digest_file(target),parts=packed.cells,checks=checks,passed=True,source_refs=source_refs)
 atomic_json(out/'native-parents.json',meta);return meta

def warp(v,role,b):
 p=np.array(v,float,copy=True)
 if role=='hepatic_vein':return p
 glob=json.loads((b/'separated-courses-v105/evidence.json').read_text());a=glob['anchors'][role];center=np.array(a['center_m']);t=np.clip((np.linalg.norm(p-center,axis=1)-a['radius_m'])/.03,0,1);p+=(t*t*(3-2*t))[:,None]*np.array(glob['offsets_m'][role]);local=json.loads((b/'local-course-clearance-v106/evidence.json').read_text())
 for f in local['fields']:
  if f['role']!=role:continue
  t=np.clip((np.linalg.norm(p-np.array(f['fixed_anchor_center_m']),axis=1)-f['fixed_anchor_radius_m'])/f['blend_m'],0,1);w=t*t*(3-2*t)*np.exp(-.5*np.sum((p-f['center_m'])**2,axis=1)/f['sigma_m']**2);p+=w[:,None]*np.array(f['vector_m'])
 return p

def prepare(b,out,atlas,cat):
 from geometry_tools.open3d_backend import library
 import importlib.util
 helper=importlib.util.spec_from_file_location('hepatic_region_labels',ROOT/'server/hepatic-region-labels.py');labels_module=importlib.util.module_from_spec(helper);helper.loader.exec_module(labels_module)
 o3d=library();a=copy.deepcopy(atlas);byid={p['id']:p for p in a['parts']};reader=AtlasReader();records=[];region_records=[];regions={};banks={};termrows={};source_to_region={};rootconcepts={};origin=np.array([0,1.25,0]);allmeshes={}
 def chunk(path):
  path=Path(path);url='/local-models/'+str(path.relative_to(ROOT/'.local-models'));i=len(a['chunks']);a['chunks'].append(dict(url=url,bytes=path.stat().st_size));return i,url
 for role in ROLES:
  layout=json.loads((b/('packed-hepatic-'+role+'-v128/tissue.json')).read_text());assert layout['passed'];assert digest_file(layout['bank'])==layout['bank_sha256'];row=layout['parts'][role];ci,url=chunk(layout['bank']);prototype=byid[PROTOTYPES[role]];uid='TA98HEPATIC:male:'+role;root='HEPATIC:male:network:'+role;rootconcepts[role]=root;data=Path(layout['bank']).read_bytes();signature=hashlib.sha256(b''.join(data[row[k]:row[k]+row['vertexCount']*size]for k,size in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps(row['geometryOrigin']).encode()).hexdigest();p={**copy.deepcopy(row),'id':uid,'chunk':ci,'conceptId':root,'name':NAMES[role],'sourceId':'Hepatic connection:male:'+role,'system':SYSTEMS[role],'material':prototype['material'],'groups':['Liver','Hepatic connections'],'pairedTissue':True,'sectionSolid':True,'tissueWallAssembly':'hepatic-connections:male:'+role,'sectionAssembly':'hepatic-connections:male:'+role,'hepaticConnectionRole':role,'hepaticSourceModel':'local-ta98','hepaticAttributeSignature':signature,'license':'Local review only; Z-Anatomy and BodyParts3D source credits retained','provenance':dict(label='Local hepatic connection assembly',url='https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html',detail='Native named branch geometry and preserved native parent interfaces. Local macroscopic repair; original source variants remain hidden.'),'coverageLimitation':LIMITS,'description':LIMITS};p['tissueAppearance']=copy.deepcopy(prototype.get('tissueAppearance',dict(color=[x/255 for x in a['materials'][prototype['material']]['color'][:3]],pattern='solid',scale=.0005,relief=0)));records.append(p);banks[role]=dict(url=url,sha256=layout['bank_sha256']);full=reader.read_part({**a,'parts':[p]},p,frame='model:local-ta98:assembled-metres',catalogue_path=cat);allmeshes[role]=full
  source=b/('term-cells-'+role+'-v90');terms=json.loads((source/'evidence.json').read_text())['parts'];scene=o3d.t.geometry.RaycastingScene(nthreads=1);gid_to_term={}
  for term in terms:
   roi=Mesh.load(term['path']);vertices=warp(roi.vertices,role,b);gid=scene.add_triangles(o3d.core.Tensor((vertices-origin).astype('<f4')),o3d.core.Tensor(roi.faces.astype('<u4')));gid_to_term[gid]=term
  centers=full.vertices[full.faces].mean(1);hits=scene.compute_closest_points(o3d.core.Tensor((centers-origin).astype('<f4')),nthreads=1);gids,selection_check=labels_module.regularize_labels(full.vertices,full.faces,hits['geometry_ids'].numpy());atomic_json(out/(role+'-selection-check.json'),selection_check);regiondata=bytearray(data);group_runs=[];regions[role]=[]
  for gid,term in gid_to_term.items():
   selected=np.flatnonzero(gids==gid);assert len(selected),('No actual triangles for TA98 region',term['term']);rid='TA98HEPATIC:male:region:'+term['term'];regions[role].append(rid)
   while len(regiondata)%4:regiondata.append(0)
   offset=len(regiondata);indices=np.ndarray((row['indexCount']//3,3),'<u4',buffer=data,offset=row['indices'])[selected];regiondata.extend(indices.astype('<u4').tobytes());used=full.vertices[full.faces[selected].ravel()];rp={**copy.deepcopy(p),'id':rid,'conceptId':rid,'name':term['name'],'ta98Term':term['term'],'indices':offset,'indexCount':int(indices.size),'triangles':len(indices),'bounds':[used.min(0).tolist(),used.max(0).tolist()],'defaultHidden':True,'pairedTissue':False,'sectionSolid':False,'hepaticRegionOf':uid,'coverageLimitation':'Owned surface region of the communicating finite network. Regional boundary provisional; not a separate microscopic wall layer.','tissueDrawGroups':[dict(start=0,count=int(indices.size))]};region_records.append(rp);termrows[term['term']]=dict(name=term['name'],id=rid,role=role)
   for oldid in term['source_ids']:
    if not oldid.startswith('@'):source_to_region[oldid]=rid
  rp_path=out/(role+'-regions.bin');atomic_bytes(rp_path,bytes(regiondata));rci,rurl=chunk(rp_path)
  for rp in region_records:
   if rp['hepaticConnectionRole']==role:rp['chunk']=rci
  labels=[next(rid for rid in regions[role]if rid.endswith(gid_to_term[int(gid)]['term']))for gid in gids];begin=0
  for i in range(1,len(labels)+1):
   if i==len(labels)or labels[i]!=labels[begin]:group_runs.append(dict(start=begin*3,count=(i-begin)*3,replacement=labels[begin]));begin=i
  p['tissueDrawGroups']=group_runs
 # Retain the actual downstream native surfaces. Only the upper common duct
 # and the blind upper portal continuation are replaced by the new junctions.
 extra=json.loads((out/'native-parents.json').read_text());assert extra['passed'];eci,eurl=chunk(extra['bank']);lower={**copy.deepcopy(byid[COMMON]),**copy.deepcopy(extra['parts']['common_lower']),'chunk':eci,'name':'Common hepatic duct · retained lower part','description':'Preserved downstream native wall, joined to the hepatic duct network.','coverageLimitation':'Native downstream cystic/bile interfaces retained; the upper junction is represented by the communicating hepatic duct assembly.','hepaticRetainedParent':True,'hepaticSourceModel':'local-ta98','hepaticNativeInterfaces':['FJ3080','TA98PART:FJ3079:A05.8.02.013']};records.append(lower);portal_lower={**copy.deepcopy(byid['ZA:Hepatic portal vein']),**copy.deepcopy(extra['parts']['portal_lower']),'id':PORTAL_DISTAL,'chunk':eci,'conceptId':'HEPATIC:male:portal-distal','name':'Hepatic portal vein · retained lower trunk','description':'Preserved native outer reference below the repaired hepatic junction.','coverageLimitation':'Native downstream outer reference retained; no new wall-layer claim for the foregut confluence.','sectionSolid':False,'hepaticRetainedParent':True,'hepaticSourceModel':'local-ta98'};records.append(portal_lower)
 for record,key in [(lower,'common_lower'),(portal_lower,'portal_lower')]:
  for attribute in ['positions','normals','normalSigns','uvs','colors','indices']:
   if attribute not in extra['parts'][key]:record.pop(attribute,None)
 # Source geometries are retained literally. They are not additional active
 # organs. The optional region views index the new wall itself.
 source_to_region.update({'FJ2816':'TA98PART:ZA:Liver:caudate','FJ3094':'TA98HEPATIC:male:region:A05.8.01.069',
  'FJ2416':'TA98HEPATIC:male:region:A12.3.09.006',
  'FJ1913':'TA98HEPATIC:male:region:A12.3.12.002','FJ2405':'TA98HEPATIC:male:region:A12.3.12.002',
  'FJ1893':'TA98HEPATIC:male:region:A12.3.12.005','FJ2404':'TA98HEPATIC:male:region:A12.3.12.005'})
 retire=set(source_to_region)|{'ZA:Common hepatic artery','ZA:Proper hepatic artery','ZA:Hepatic portal vein','ZA:Hepatic veins','ZA:Inferior vena cava (thoracic part)'}
 for pid in retire:
  if pid in byid:byid[pid].update(suppressed=True,defaultHidden=True,reviewStatus='Retained source; replaced by communicating male hepatic assembly')
 a['parts']=[p for p in a['parts']if p['id']!=COMMON]+records+region_records
 a['defaultHidden']=sorted(set(a.get('defaultHidden',[]))|retire|{p['id']for p in region_records})
 a['concepts']=[c for c in a['concepts']if c['id']not in {p['conceptId']for p in region_records}|set(rootconcepts.values())]
 for role in ROLES:a['concepts'].append(dict(id=rootconcepts[role],name=NAMES[role],elements=['TA98HEPATIC:male:'+role],partOf='ZA:Liver',hepaticComponents=regions[role]))
 a['concepts'].append(dict(id='HEPATIC:male:portal-distal',name=portal_lower['name'],elements=[PORTAL_DISTAL]))
 parents=json.loads((ROOT/'app/data/ta98-parents.json').read_text())['parents'];groups=[]
 for term,r in termrows.items():
  ids=[]
  for other,s in termrows.items():
   t=other;seen=set()
   while t and t not in seen:
    if t==term:ids.append(s['id']);break
    seen.add(t);t=parents.get(t)
  if term=='A12.2.12.015':ids=['TA98HEPATIC:male:artery']
  elif term=='A12.3.12.001':ids=['TA98HEPATIC:male:portal',PORTAL_DISTAL]
  elif term=='A05.8.01.061':ids=['TA98HEPATIC:male:biliary',COMMON]
  elif term=='A12.3.09.001':ids=['ZA:Inferior vena cava (abdominal part)',r['id']]
  cid='HEPATIC:male:term:'+term;groups.append(dict(term=term,conceptId=cid,parts=ids));a['concepts'].append(dict(id=cid,name=r['name'],elements=ids,ta98Term=term,ta98Kind='exact',coverageLimitation=LIMITS))
  a['concepts'].append(dict(id=r['id'],name=r['name'],elements=[r['id']],ta98Term=term,ta98Kind='exact',partOf=rootconcepts[r['role']],coverageLimitation='Named surface region; provisional regional boundary.'))
 alias={'ZA:Common hepatic artery':['TA98HEPATIC:male:artery'],'ZA:Proper hepatic artery':[rid for rid in regions['artery']if not rid.endswith('A12.2.12.015')],'ZA:Hepatic portal vein':['TA98HEPATIC:male:portal',PORTAL_DISTAL],'ZA:Hepatic veins':['TA98HEPATIC:male:hepatic_vein'],'ZA:Inferior vena cava (thoracic part)':['TA98HEPATIC:male:region:A12.3.09.001'],'TA98BILIARY:A05.8.01.061':['TA98HEPATIC:male:biliary',COMMON]}
 for c in a['concepts']:
  if c['id']in alias:c['elements']=alias[c['id']]
  else:c['elements']=list(dict.fromkeys(source_to_region.get(pid,pid)for pid in c.get('elements',[])))
 a['hepaticConnectionIntegration']=dict(schema='male-hepatic-connections/v1',sourceCatalogueSha256=digest_file(cat),wholePartIds=['TA98HEPATIC:male:'+r for r in ROLES],regionPartIds=[p['id']for p in region_records],termGroups=groups,banks=banks,regionAllocation='Disjoint index subsets of actual communicating wall; optional, hidden initially',remainingLimits=LIMITS,anatomicallyAccepted=False)
 a['version']+=' · male hepatic connections';a['triangles']=sum(p['indexCount']//3 for p in a['parts']if not p.get('suppressed'));target=out/'atlas.json';atomic_json(target,a);before={p['id']:p for p in atlas['parts']};changed=[p['id']for p in a['parts']if before.get(p['id'])!=p];atomic_json(out/'candidate.json',dict(proposed_catalogue=str(target),expected_catalogue_sha256=digest_file(cat),changed_parts=changed,whole_parts=4,ta98_surface_regions=len(region_records),source_liver_tissue_unchanged=all(before[p['id']]==p for p in a['parts']if p['id'].startswith('TA98PART:ZA:Liver:')),external_deployment=False));print(json.dumps(dict(parts=len(a['parts']),changed=len(changed),regions=len(region_records),output=str(target))))

def main():
 p=argparse.ArgumentParser();p.add_argument('--workbench',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--stage',choices=['extras','candidate'],required=True);args=p.parse_args();b=args.workbench.resolve();out=args.output.resolve();out.mkdir(exist_ok=True);cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(cat)=='a24acee189983bc9b33899ebf0c9797db79683778a5574ebe5d3c461a830fa5c';a=json.loads(cat.read_text())
 if args.stage=='extras':assert not(out/'native-parents.json').exists();extras(b,out,a,cat)
 else:assert not(out/'atlas.json').exists();prepare(b,out,a,cat)
if __name__=='__main__':main()
