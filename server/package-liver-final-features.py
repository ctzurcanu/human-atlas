"""Package audited complementary liver cells and rebind actual surface witnesses."""
from pathlib import Path
import argparse,copy,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/args.stage;out.mkdir(exist_ok=False);reader=AtlasReader();mp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';fp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(mp.read_text());female=json.loads(fp.read_text());old=copy.deepcopy(a);byid={p['id']:p for p in a['parts']}
proc=json.loads((b/'final-processes-v223/processes.json').read_text());assert proc['passed'];packed=json.loads((b/'final-processes-v223/packed/tissue.json').read_text());assert packed['passed'];origin=np.array([0,1.25,0.]);oldbank=reader.asset(a['liverIntegration']['bankUrl']);data=bytearray(oldbank.read_bytes());offset=len(data);data.extend(Path(proc['bank']).read_bytes());records=copy.deepcopy(packed['parts'])
def append(arr,bank=data):
 bank.extend(b'\0'*((-len(bank))%4));off=len(bank);bank.extend(arr.tobytes());return off
names={'serosa':('Tunica serosa','A05.8.01.052'),'subserosa':('Tela subserosa','A05.8.01.053'),'fibrous-appendix':('Fibrous appendix of liver','A05.8.01.028'),'papillary-process':('Papillary process of caudate lobe','A05.8.01.031'),'caudate-process':('Caudate process','A05.8.01.032')}
colors={'serosa':[.63,.31,.29],'subserosa':[.68,.43,.36],'fibrous-appendix':[.67,.57,.45]}
ids=[]
for key,rec in records.items():
 pid=rec['id'];ids.append(pid)
 if pid not in byid:
  row=copy.deepcopy(byid['TA98PART:ZA:Liver:capsule']);row.update(id=pid,conceptId=pid,name=names[key][0],sourceId='Liver:'+key,ta98Term=names[key][1],tissueWallLayer=key,material='LIVER:male:'+key);a['parts'].append(row);byid[pid]=row
 row=byid[pid];row.update(rec)
 for k in ['positions','normals','indices','normalSigns']:row[k]+=offset
 row.update(defaultHidden=False,sectionSolid=key not in ['right','quadrate'],pairedTissue=True,coverageLimitation='Local finite tissue and gross regions. Regional process boundaries and histological thickness allocation are provisional; microscopic unit is illustrative.',reviewStatus='local-finite-liver-v2')
 row['description']=row['coverageLimitation'];row.pop('suppressed',None);row.pop('sourceOffset',None)
 if key in colors:row['tissueAppearance']={'color':colors[key],'pattern':'fibrous','relief':.000002,'scale':.0004};a['materials'][row['material']]={'color':[round(c*255)for c in colors[key]],'roughness':.4 if key=='serosa'else.7,'metalness':0}
host=byid['TA98PART:ZA:Liver:capsule'];host['liverReferenceTriangleCount']=proc['referenceSurfaceTriangleCount'];ancestry=np.array(proc['host_source_faces'],int)
# Source computational classifiers for the two unchanged hilar-open lobes are
# literal. The subdivided closed lobes use their actual new finite cells.
section=a['liverSectionIntegration'];kept={v['id']:v for v in section['volumes']};volumes=[]
for key in ['right','left','quadrate','caudate','fibrous-appendix','papillary-process','caudate-process']:
 if key in ['right','quadrate']:volumes.append(kept[key]);continue
 p=byid['TA98PART:ZA:Liver:'+key];m=Mesh.load(b/'final-processes-v223/packed'/(key+'.npz'))
 positions=(m.vertices-origin).astype('<f4');em=Mesh(positions.astype(float)+origin,m.faces,'model:local-ta98:assembled-metres')
 volumes.append(dict(id=key,partId=p['id'],positions=append(positions),indices=append(m.faces.astype('<u4')),vertexCount=len(positions),indexCount=m.faces.size,origin=origin.tolist(),geometrySha256=em.geometry_sha))
section['volumes']=volumes
bank=out/'liver.bin';bank.write_bytes(data);url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'));sha=digest_file(bank);chunk=len(a['chunks']);a['chunks'].append(dict(url=url,bytes=len(data)))
for pid in ids:byid[pid]['chunk']=chunk
li=a['liverIntegration'];li.update(schema='male-liver-lobes-capsule/v2',bankUrl=url,bankSha256=sha,partIds=ids,remainingLimits='69/69 liver terms represented: 63 organ-scale regions, tissues, ducts or bound references; six microscopic features in an illustrative unit. Process/functional boundaries and coat thicknesses remain provisional.')
groups=[]
termnames={'right':'A05.8.01.026','left':'A05.8.01.027','quadrate':'A05.8.01.029','caudate':'A05.8.01.030','capsule':'A05.8.01.054',**{k:v[1]for k,v in names.items()}}
groups.append(dict(term='A05.8.01.001',conceptId='ZA:Liver',parts=ids))
for key,term in termnames.items():
 pid='TA98PART:ZA:Liver:'+key;children=[pid]+(['TA98PART:ZA:Liver:fibrous-appendix']if key=='left'else['TA98PART:ZA:Liver:papillary-process','TA98PART:ZA:Liver:caudate-process']if key=='caudate'else[]);groups.append(dict(term=term,conceptId=pid,parts=children))
li['termGroups']=groups;cs={c['id']:c for c in a['concepts']}
for g in groups:
 c=cs.get(g['conceptId'])
 if c is None:c=dict(id=g['conceptId'],name=byid[g['conceptId']]['name'],partOf='ZA:Liver');a['concepts'].append(c);cs[c['id']]=c
 c.update(elements=g['parts'],ta98Term=g['term'],ta98Kind='exact')
for c in a['concepts']:
 if c['id']=='ZA:collection:Liver':c['elements']=ids
section.update(bankUrl=url,bankSha256=sha)
features=a['liverFeatureIntegration'];features['bankSha256']=sha;completion=a['liverCompletionIntegration'];completion['hostBankSha256']=sha
v=np.ndarray((host['vertexCount'],3),'<f4',buffer=data,offset=host['positions']).astype(float)+origin;f=np.ndarray((host['indexCount']//3,3),'<u4',buffer=data,offset=host['indices']);ref=f[:len(ancestry)]
import open3d as o3d
scene=o3d.t.geometry.RaycastingScene();scene.add_triangles(o3d.core.Tensor((v-origin).astype('float32')),o3d.core.Tensor(ref.astype('uint32')))
def witness(w):
 point=np.array(w['point']);q=scene.compute_closest_points(o3d.core.Tensor((point-origin)[None,:].astype('float32')));fi=int(q['primitive_ids'].numpy()[0]);uv=q['primitive_uvs'].numpy()[0].astype(float);weights=np.r_[1-uv.sum(),uv];weights=np.clip(weights,0,1);weights/=weights.sum();w.update(face=fi,barycentric=weights.tolist(),point=(weights@v[f[fi]]).tolist())
def rebind(x):
 if isinstance(x,dict):
  if x.get('partId')==host['id']and'face'in x and'point'in x:witness(x)
  for value in x.values():rebind(value)
 elif isinstance(x,list):
  for value in x:rebind(value)
for g in features['termGroups']:
 g['surface']['faces']=np.flatnonzero(np.isin(ancestry,g['surface']['faces'])).tolist();g['surface']['sourceBankSha256']=sha
rebind(features);rebind(completion)
# Functional buffers retain the same coordinates, with the classifier bank
# literally copied into their existing bank so both section paths agree.
fs=completion['functionalSections'];functional=bytearray(reader.asset(fs['bankUrl']).read_bytes());oldraw=oldbank.read_bytes()
for dest,source in zip(fs['exclusions'],section['exclusions']):
 for k,mult in [('positions',12),('indices',4)]:
  n=source['vertexCount']if k=='positions'else source['indexCount'];functional.extend(b'\0'*((-len(functional))%4));dest[k]=len(functional);functional.extend(data[source[k]:source[k]+n*mult])
 dest.update({k:source[k]for k in ['vertexCount','indexCount','origin','geometrySha256']})
funpath=out/'functional.bin';funpath.write_bytes(functional);funurl='/local-models/'+str(funpath.relative_to(ROOT/'.local-models'));fs.update(bankUrl=funurl,bankSha256=digest_file(funpath),sourceBankSha256=sha);ci=len(a['chunks']);a['chunks'].append(dict(url=funurl,bytes=len(functional)))
for pid in completion['functionalPartIds']:byid[pid]['chunk']=ci
atomic_json(out/'male.json',a)
# Transfer the completed donor assembly by one uniform scale and translation.
fit=json.loads((b/'final-context-v217/fit.json').read_text())['chosen'];oldtransfer=female['liverTransferIntegration'];S=oldtransfer['scale']*fit['additional_scale'];T=np.array(oldtransfer['translation_m'])*fit['additional_scale']+(1-fit['additional_scale'])*np.array(fit['center_m'])+fit['translation_m'];fb={p['id']:p for p in female['parts']};fc={c['id']:c for c in female['concepts']};chunks={c['url']:i for i,c in enumerate(female['chunks'])}
def remap(value):
 if isinstance(value,str):return value.replace(':male:',':female:').replace('male-liver','female-liver').replace('male-hepatic','female-hepatic').replace('ZA:Liver','LOCAL:female:liver')if not value.startswith('TA98PART:ZA:Liver:')else value
 if isinstance(value,list):return[remap(x)for x in value]
 if isinstance(value,dict):return{k:remap(x)for k,x in value.items()}
 return value
transferids=set(oldtransfer['partIds'])|set(ids)
for source in a['parts']:
 pid=remap(source['id'])
 if pid not in transferids:continue
 p=remap(copy.deepcopy(source));p['liverTransferSourcePartId']=source['id'];p['geometryScale']=S;p['geometryOrigin']=(np.array(source.get('geometryOrigin',[0,0,0]))*S+T).tolist();p['bounds']=(np.array(source['bounds'])*S+T).tolist();ch=a['chunks'][source['chunk']]
 if ch['url']not in chunks:chunks[ch['url']]=len(female['chunks']);female['chunks'].append(copy.deepcopy(ch))
 p['chunk']=chunks[ch['url']]
 if source.get('sourceOffset'):p['sourceOffset']=(np.array(source['sourceOffset'])*S).tolist()
 if pid in fb:fb[pid].clear();fb[pid].update(p)
 else:female['parts'].append(p);fb[pid]=p
 if p['material']in a['materials']:female['materials'][p['material']]=copy.deepcopy(a['materials'][p['material']])
for key in ['liverIntegration','liverFeatureIntegration','liverCompletionIntegration','liverSectionIntegration']:
 female[key]=remap(copy.deepcopy(a[key]))
def transform(x):
 if isinstance(x,dict):
  if'point'in x and'barycentric'in x:x['point']=(np.array(x['point'])*S+T).tolist()
  if'positions'in x and'indices'in x and'origin'in x:
   x['origin']=(np.array(x['origin'])*S+T).tolist();x['geometryScale']=S
  for val in x.values():transform(val)
 elif isinstance(x,list):
  for val in x:transform(val)
for key in ['liverFeatureIntegration','liverCompletionIntegration','liverSectionIntegration']:transform(female[key])
for g in female['liverIntegration']['termGroups']:
 c=fc.get(g['conceptId'])
 if c is None:c=dict(id=g['conceptId'],name=fb[g['conceptId']]['name'],partOf='LOCAL:female:liver');female['concepts'].append(c);fc[c['id']]=c
 c.update(elements=g['parts'],ta98Term=g['term'],ta98Kind='exact')
oldtransfer.update(scale=S,geometryScale=S,translation_m=T.tolist(),partIds=sorted(transferids),finalContextFit=str(b/'final-context-v217/fit.json'))
for p in female['parts']:
 if p.get('hepaticConnectionRole'):
  raw=reader.asset(female['chunks'][p['chunk']]['url']).read_bytes();signature=hashlib.sha256(b''.join(raw[p[k]:p[k]+p['vertexCount']*n]for k,n in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps([p.get('geometryOrigin'),p.get('geometryScale',1)]).encode()).hexdigest();p['hepaticAttributeSignature']=signature
# All witnesses are replayed from the emitted transformed arrays, eliminating
# floating point differences caused by changing operation order.
fh=fb['TA98PART:ZA:Liver:capsule'];fv=(v-origin)*S+np.array(fh['geometryOrigin'])
def replay(x):
 if isinstance(x,dict):
  if x.get('partId')==fh['id']and'face'in x and'barycentric'in x:x['point']=(np.array(x['barycentric'])@fv[f[x['face']]]).tolist()
  for value in x.values():replay(value)
 elif isinstance(x,list):
  for value in x:replay(value)
replay(female['liverFeatureIntegration']);replay(female['liverCompletionIntegration'])
atomic_json(out/'female.json',female);atomic_json(out/'features.json',dict(passed=True,source_male_sha256=digest_file(mp),source_female_sha256=digest_file(fp),new_tissue_parts=ids,source_geometry_bank_sha256=proc['bank_sha256'],reference_facets=len(ancestry),female_scale=S,female_translation_m=T.tolist(),elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(male=str(out/'male.json'),female=str(out/'female.json'),features_added=5,female_scale=S,seconds=time.monotonic()-start)))
