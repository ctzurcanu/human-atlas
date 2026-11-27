"""Pack passed hepatic networks with exclusive term views and reciprocal mouths."""
from pathlib import Path
import argparse,copy,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();sex=args.sex;start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/args.stage;out.mkdir(exist_ok=False);a=json.loads((b/'final-package-v246'/(sex+'.json')).read_text());parts={p['id']:p for p in a['parts']};cs={c['id']:c for c in a['concepts']};r=AtlasReader();frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres';S=a['liverTransferIntegration']['scale']if sex=='female'else 1.;T=np.array(a['liverTransferIntegration']['translation_m'])if sex=='female'else np.zeros(3);O=np.array([0,1.25,0.]);world_origin=O*S+T;interfaces=[];changed=[]
for p in a['parts']:
 if p.get('hepaticConnectionRole')or p.get('hepaticRetainedParent'):p['hepaticSourceModel']='local-female-ta98'if sex=='female'else'local-ta98'
def usebank(bank):
 url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'));i=len(a['chunks']);a['chunks'].append(dict(url=url,bytes=bank.stat().st_size));return i,url
def encode(pid,mesh,neighbors=None,labels=None):
 """Smooth native surfaces; exact shared charts use canonical signed words."""
 row=parts[pid];data=bytearray();v=mesh.vertices;f=mesh.faces.astype(int);surfaces,_,proof=conform_quantized_surfaces(v,{'mesh':f},O);local=(v-O).astype('<f4');v=local.astype(float)+O;stored=v[f];valid=np.linalg.norm(np.cross(stored[:,1]-stored[:,0],stored[:,2]-stored[:,0]),axis=1)>0;old_to_new=np.cumsum(valid)-1
 if labels is not None:labels=np.asarray(labels,object)[valid]
 if neighbors is not None:neighbors={peer:old_to_new[np.asarray(ids,int)[valid[np.asarray(ids,int)]]].tolist()for peer,ids in neighbors.items()}
 f=surfaces['mesh'];assert len(f)==int(valid.sum());tri=v[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);vn=np.zeros_like(v)
 for j in range(3):np.add.at(vn,f[:,j],fn)
 vn/=np.maximum(np.linalg.norm(vn,axis=1)[:,None],1e-30);words=np.rint(vn*32767).astype('<i2');sign=np.ones(len(v),dtype='i1');newlocal=local.tolist();newwords=words.tolist();newsign=sign.tolist();ff=f.copy();neighbors=neighbors or{}
 canonical=set(i for ids in neighbors.values()for i in ids)
 for fi,t in enumerate(f):
  if fi in canonical:
   order=np.lexsort((tri[fi,:,2],tri[fi,:,1],tri[fi,:,0]));p=tri[fi,order];n=np.cross(p[1]-p[0],p[2]-p[0]);n/=np.linalg.norm(n);w=np.rint(n*32767).astype('<i2');sg=1 if np.dot(n,fn[fi])>0 else-1
   for j,vid in enumerate(t):ff[fi,j]=len(newlocal);newlocal.append(local[vid].tolist());newwords.append(w.tolist());newsign.append(sg)
  else:
   bad=np.flatnonzero((words[t].astype(float)@fn[fi])<=0)
   for j in bad:
    ff[fi,j]=len(newlocal);newlocal.append(local[t[j]].tolist());newwords.append(np.rint(fn[fi]/np.linalg.norm(fn[fi])*32767).astype('<i2').tolist());newsign.append(1)
 def add(arr):
  data.extend(b'\0'*((-len(data))%4));off=len(data);data.extend(arr.tobytes());return off
 layout=dict(positions=add(np.array(newlocal,'<f4')),normals=add(np.array(newwords,'<i2')),normalSigns=add(np.array(newsign,'i1')),indices=add(ff.astype('<u4')),vertexCount=len(newlocal),indexCount=ff.size,triangles=len(ff),geometryOrigin=world_origin.tolist(),geometryScale=S)
 groups=[];owner=np.full(len(f),'',object)
 for peer,ids in neighbors.items():owner[np.asarray(ids,int)]=peer
 replacement=np.asarray(labels,object)if labels is not None else np.full(len(f),'',object)
 for i,(peer,rep)in enumerate(zip(owner,replacement)):
  if groups and groups[-1].get('neighbor','')==peer and groups[-1].get('replacement','')==rep:groups[-1]['count']+=3
  else:
   g=dict(start=i*3,count=3)
   if peer:g['neighbor']=peer
   if rep:g['replacement']=rep
   groups.append(g)
 bank=out/(pid.replace(':','_').replace('/','_')+'.bin');bank.write_bytes(data);ci,url=usebank(bank);row.update(layout,chunk=ci,tissueDrawGroups=groups,pairedTissue=True,sectionSolid=True,bounds=(np.array([v[np.unique(f)].min(0),v[np.unique(f)].max(0)])*S+T).tolist());row.pop('sourceOffset',None);row.pop('uvs',None);row.pop('colors',None);row.pop('gastricLumen',None);row.pop('tissueChart',None);row['reviewStatus']='local-patent-hepatic-interfaces';row['description']='Finite hepatic conduit; native distal circulation beyond the declared interface is outside this local repair.';changed.append(pid)
 if labels is not None:
  signature=hashlib.sha256(b''.join(data[layout[k]:layout[k]+len(newlocal)*n]for k,n in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps([world_origin.tolist(),S]).encode()).hexdigest();row['hepaticAttributeSignature']=signature
  for rid in sorted(set(replacement)-{''}):
   target=parts[rid];indices=ff[replacement==rid];off=add(indices.astype('<u4'));target.update({k:layout[k]for k in ['positions','normals','normalSigns','vertexCount','geometryOrigin','geometryScale']},indices=off,indexCount=indices.size,triangles=len(indices),chunk=ci,hepaticAttributeSignature=signature,tissueDrawGroups=[dict(start=0,count=indices.size)],sectionSolid=False,pairedTissue=True);target['bounds']=row['bounds'];target.pop('uvs',None);target.pop('sourceOffset',None);changed.append(rid)
  bank.write_bytes(data);a['chunks'][ci]['bytes']=len(data);a['hepaticConnectionIntegration']['banks'][row['hepaticConnectionRole']]=dict(url=url,sha256=digest_file(bank))
 return row
def read_current(pid):return r.read_part(a,parts[pid],frame=frame)
def donor(m):return Mesh((m.vertices-T)/S,m.faces,'construction:hepatic:donor-metres')if sex=='female'else m
def oldlabels(pid,new):
 old=donor(read_current(pid));row=parts[pid];labels=np.full(len(old.faces),'',object)
 for g in row['tissueDrawGroups']:
  if g.get('replacement'):labels[g['start']//3:(g['start']+g['count'])//3]=g['replacement']
 import open3d as o3d
 scene=o3d.t.geometry.RaycastingScene();scene.add_triangles(o3d.core.Tensor((old.vertices-O).astype('float32')),o3d.core.Tensor(old.faces.astype('uint32')));q=scene.compute_closest_points(o3d.core.Tensor((new.triangles.mean(1)-O).astype('float32')));result=labels[q['primitive_ids'].numpy()];assert np.all(result!='');return result
# Add the sixth term as a real branch view of the communicating biliary wall.
rid='TA98HEPATIC:'+sex+':region:A05.8.01.068';template=copy.deepcopy(parts['TA98HEPATIC:'+sex+':region:A05.8.01.069']);template.update(id=rid,conceptId='HEPATIC:'+sex+':term:A05.8.01.068',name='Right duct of caudate lobe',ta98Term='A05.8.01.068');a['parts'].append(template);parts[rid]=template;a['concepts'].append(dict(id=template['conceptId'],name=template['name'],elements=[rid],ta98Term='A05.8.01.068',ta98Kind='exact',partOf='HEPATIC:'+sex+':term:A05.8.01.064'));a['hepaticConnectionIntegration']['regionPartIds'].append(rid);a['hepaticConnectionIntegration']['termGroups'].append(dict(term='A05.8.01.068',conceptId=template['conceptId'],parts=[rid]))
role_neighbors={};network={};parents=[]
if sex=='male':
 for role in ['artery','hepatic_vein']:
  q=json.loads((b/'male-parent-walls-v230'/(role+'-result.json')).read_text());assert q['passed'];network[role]=read_current('TA98HEPATIC:male:'+role);role_neighbors[role]={q['parent']:q['network_interface_faces']};parent=Mesh.load(b/'male-parent-walls-v230'/role/(q['parent'].replace(':','_')+'.npz'));parents.append((q['parent'],parent,{'TA98HEPATIC:male:'+role:q['parent_interface_faces']}));interfaces.append(dict(first=q['parent'],second='TA98HEPATIC:male:'+role,first_faces=q['parent_interface_faces'],second_faces=q['network_interface_faces']))
 q=json.loads((b/'common-hepatic-mouth-v243/junction.json').read_text());assert q['passed'];network['biliary']=Mesh.load(b/'common-hepatic-mouth-v243/child.npz');pid='TA98PART:FJ3079:A05.8.01.061';role_neighbors['biliary']={pid:q['child_interface_faces']};parents.append((pid,Mesh.load(b/'common-hepatic-mouth-v243/parent.npz'),{'TA98HEPATIC:male:biliary':q['parent_interface_faces']}));interfaces.append(dict(first=pid,second='TA98HEPATIC:male:biliary',first_faces=q['parent_interface_faces'],second_faces=q['child_interface_faces']))
else:
 for role,stage in [('artery','female-inlets-v242-artery'),('portal','female-inlets-v242-portal'),('biliary','female-inlets-v249-biliary')]:
  q=json.loads((b/stage/'junction.json').read_text());assert q['passed'];pid=q['parent_id'];network[role]=Mesh.load(b/stage/'child.npz');role_neighbors[role]={pid:q['child_interface_faces']};parents.append((pid,Mesh.load(b/stage/'parent.npz'),{'TA98HEPATIC:female:'+role:q['parent_interface_faces']}));interfaces.append(dict(first=pid,second='TA98HEPATIC:female:'+role,first_faces=q['parent_interface_faces'],second_faces=q['child_interface_faces']))
 # Carry the one 18-micron native wall-control correction through the exact
 # shared neighbouring chart, rather than breaking an existing dependency.
 mapping=json.loads((b/'female-inlets-v249-biliary/shared-control-map.json').read_text());src=np.array(mapping['source']);dst=np.array(mapping['destination']);row=parts['LOCAL:female:bile_duct:1031'];neighbors={g['neighbor']for g in row['tissueDrawGroups']if g.get('neighbor')}
 for pid in sorted(neighbors):
  m=donor(read_current(pid));v=m.vertices.copy();keys={(tuple(((x-O).astype('<f4').astype(float)+O))):y for x,y in zip(src,dst)};moves=0
  for i,x in enumerate(v):
   key=tuple((x-O).astype('<f4').astype(float)+O)
   if key in keys:v[i]=keys[key];moves+=1
  if moves:parents.append((pid,Mesh(v,m.faces,m.frame),{}))
 for pid in ['TA98HEPATIC:female:portal-distal','TA98HEPATIC:female:common-hepatic-parent']:
  parts[pid].update(suppressed=True,defaultHidden=True,replacementParts=['LOCAL:female:portal_vein:1795']if pid.endswith('portal-distal')else['LOCAL:female:bile_duct:1031']);changed.append(pid)
  for marker in ['hepaticConnectionIntegration']:
   for g in a[marker]['termGroups']:g['parts']=[parts[pid]['replacementParts'][0]if x==pid else x for x in g['parts']]
   for c in a['concepts']:
    c['elements']=[parts[pid]['replacementParts'][0]if x==pid else x for x in c.get('elements',[])]
# Explicit duct branch ownership on the single shared wall bank.
for role,new in network.items():
 pid='TA98HEPATIC:'+sex+':'+role;labels=oldlabels(pid,new)
 if role=='biliary':
  path=np.array(json.loads((b/'right-caudate-duct-v226/construction.json').read_text())['waypoints']);points=new.triangles.mean(1);distance=np.full(len(points),np.inf);fraction=np.zeros(len(points))
  for i,(x,y)in enumerate(zip(path[:-1],path[1:])):
   t=np.clip((points-x)@(y-x)/np.dot(y-x,y-x),0,1);d=np.linalg.norm(points-(x+t[:,None]*(y-x)),axis=1);take=d<distance;distance[take]=d[take];fraction[take]=(i+t[take])/(len(path)-1)
  labels[(distance<.00085)&(fraction>.025)]=rid
  from importlib.util import spec_from_file_location,module_from_spec
  sp=spec_from_file_location('hepatic_regions',ROOT/'server/hepatic-region-labels.py');mod=module_from_spec(sp);sp.loader.exec_module(mod);values=sorted(set(labels));regular,proof=mod.regularize_labels(new.vertices,new.faces,np.array([values.index(x)for x in labels]));labels=np.asarray(values,object)[regular];atomic_json(out/'branch-labels.json',proof)
 encode(pid,new,role_neighbors[role],labels)
 for c in a['concepts']:
  if c.get('elements')==[pid]and'hepaticComponents'in c:c['hepaticComponents']=sorted(set(labels))
for pid,m,ng in parents:encode(pid,m,ng)
# Source curves/regions refer to the same bank prefix after adding the branch
# envelope to both the gross and functional section classifier banks.
envelope=Mesh.load(b/'right-caudate-duct-v226/outer.npz');positions=(envelope.vertices-O).astype('<f4');encoded=Mesh(positions.astype(float)+O,envelope.faces,'model:local-ta98:assembled-metres')
def update_section(marker,path,sourcehash=None):
 raw=bytearray(r.asset(marker['bankUrl']).read_bytes());desc=next(v for v in marker['exclusions']if v['id']=='biliary')
 for key,arr in [('positions',positions),('indices',envelope.faces.astype('<u4'))]:raw.extend(b'\0'*((-len(raw))%4));desc[key]=len(raw);raw.extend(arr.tobytes())
 desc.update(vertexCount=len(positions),indexCount=envelope.faces.size,geometrySha256=encoded.geometry_sha,origin=world_origin.tolist(),geometryScale=S);path.write_bytes(raw);ci,url=usebank(path);marker.update(bankUrl=url,bankSha256=digest_file(path))
 if sourcehash:marker['sourceBankSha256']=sourcehash
 return ci,url,digest_file(path)
ci,url,sha=update_section(a['liverSectionIntegration'],out/'liver-sections.bin')
for pid in a['liverIntegration']['partIds']:parts[pid]['chunk']=ci
a['liverIntegration'].update(bankUrl=url,bankSha256=sha);a['liverFeatureIntegration']['bankSha256']=sha;a['liverCompletionIntegration']['hostBankSha256']=sha
for g in a['liverFeatureIntegration']['termGroups']:g['surface']['sourceBankSha256']=sha
ci,_,_=update_section(a['liverCompletionIntegration']['functionalSections'],out/'functional-sections.bin',sha)
for pid in a['liverCompletionIntegration']['functionalPartIds']:parts[pid]['chunk']=ci
represented={g['term']for key in ['liverIntegration','liverFeatureIntegration','liverCompletionIntegration','hepaticConnectionIntegration']for g in a[key]['termGroups']if g['term'].startswith('A05.8.01.')};assert len(represented)==69
a['liverFinalIntegration']=dict(schema='liver-completion-interfaces/v1',representedTerms=sorted(represented),interfaces=interfaces,sourceCatalogueSha256=digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json'),fixed=['skeleton','diaphragm','stomach'],femaleFinalFit=str(b/'final-context-v217/fit.json')if sex=='female'else None,newBranch='A05.8.01.068',microscopicIllustrative=True,anatomicalAcceptance=False)
atomic_json(out/'atlas.json',a);atomic_json(out/'packaging.json',dict(passed=True,sex=sex,represented_liver_terms=len(represented),changed=changed,interfaces=interfaces,elapsed_seconds=time.monotonic()-start));print(sex,69,'terms packaged',len(changed),'updated meshes',flush=True)
