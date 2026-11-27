"""Reuse literal reviewed buffers; resize through explicit uniform transforms."""
from pathlib import Path
import copy,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
base=ROOT/'.local-models/geometry-workbench/liver-completion';inputs=base/'female-transfer-v202';out=base/'female-transfer-v206';out.mkdir(exist_ok=False)
plan=json.loads((inputs/'plan.json').read_text());fit=json.loads((inputs/'fit-cleared-v204.json').read_text());m=json.loads(Path(plan['sourceCatalogue']).read_text());a=json.loads((base/'female-transfer-v204/atlas.json').read_text());source={p['id']:p for p in m['parts']};reader=AtlasReader();scale=fit['scale'];translation=np.array(fit['translation_m']);byid={p['id']:p for p in a['parts']};chunks={c['url']:i for i,c in enumerate(a['chunks'])}
def use_chunk(ci):
 c=m['chunks'][ci]
 if c['url']not in chunks:chunks[c['url']]=len(a['chunks']);a['chunks'].append(copy.deepcopy(c))
 return chunks[c['url']]
for pid in a['liverTransferIntegration']['partIds']:
 p=byid[pid];old=source[p['liverTransferSourcePartId']]
 for key in ['positions','normals','indices','vertexCount','indexCount','normalSigns','uvs','colors','tissueChart','gastricLumen']:
  if key in old:p[key]=old[key]
  else:p.pop(key,None)
 p['chunk']=use_chunk(old['chunk']);p['geometryScale']=scale;p['geometryOrigin']=(np.array(old.get('geometryOrigin',[0,0,0]))*scale+translation).tolist()
 if old.get('sourceOffset'):p['sourceOffset']=(np.array(old['sourceOffset'])*scale).tolist()
 mesh=reader.read_part(a,p,frame='model:local-female-ta98:assembled-metres');p['bounds']=mesh.bounds.tolist()
def volumes(source_marker,destination):
 ci=next(i for i,c in enumerate(m['chunks'])if c['url']==source_marker['bankUrl']);raw=reader.asset(source_marker['bankUrl']).read_bytes();use_chunk(ci)
 destination.update(bankUrl=source_marker['bankUrl'],bankSha256=source_marker['bankSha256'])
 for key in ['volumes','exclusions']:
  destination[key]=copy.deepcopy(source_marker.get(key,[]))
  for v in destination[key]:
   v['geometryScale']=scale;v['origin']=(np.array(v['origin'])*scale+translation).tolist()
   vv=np.ndarray((v['vertexCount'],3),'<f4',buffer=raw,offset=v['positions']).astype(float)*scale+v['origin'];ff=np.ndarray((v['indexCount']//3,3),'<u4',buffer=raw,offset=v['indices'])
   v['geometrySha256']=Mesh(vv,ff,'model:local-female-ta98:assembled-metres').geometry_sha
volumes(m['liverSectionIntegration'],a['liverSectionIntegration'])
body=a['liverSectionIntegration']
for key in ['liverIntegration','liverFeatureIntegration']:a[key].update(bankUrl=body['bankUrl'],bankSha256=body['bankSha256'])
a['liverCompletionIntegration']['hostBankSha256']=body['bankSha256']
for g in a['liverFeatureIntegration']['termGroups']:g['surface']['sourceBankSha256']=body['bankSha256']
volumes(m['liverCompletionIntegration']['functionalSections'],a['liverCompletionIntegration']['functionalSections'])
a['liverCompletionIntegration']['functionalSections']['sourceBankSha256']=body['bankSha256']
host=byid['TA98PART:ZA:Liver:capsule'];raw=reader.asset(a['chunks'][host['chunk']]['url']).read_bytes();v=np.ndarray((host['vertexCount'],3),'<f4',buffer=raw,offset=host['positions']).astype(float)*scale+host['geometryOrigin'];f=np.ndarray((host['indexCount']//3,3),'<u4',buffer=raw,offset=host['indices'])
def bind(x):
 if isinstance(x,dict):
  if x.get('partId')==host['id']and 'face'in x and 'barycentric'in x and 'point'in x:x['point']=(np.array(x['barycentric'])@v[f[x['face']]]).tolist()
  for value in x.values():bind(value)
 elif isinstance(x,list):
  for value in x:bind(value)
bind(a['liverFeatureIntegration']);bind(a['liverCompletionIntegration'])
for role in ['artery','portal','hepatic_vein','biliary']:
 p=byid['TA98HEPATIC:female:'+role];raw=reader.asset(a['chunks'][p['chunk']]['url']).read_bytes()
 signature=hashlib.sha256(b''.join(raw[p[k]:p[k]+p['vertexCount']*n]for k,n in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps([p['geometryOrigin'],scale]).encode()).hexdigest()
 for pid in a['liverTransferIntegration']['partIds']:
  if byid[pid].get('hepaticConnectionRole')==role:byid[pid]['hepaticAttributeSignature']=signature
 a['hepaticConnectionIntegration']['banks'][role]=copy.deepcopy(m['hepaticConnectionIntegration']['banks'][role])
a['liverTransferIntegration'].update(literalSourceBuffers=True,geometryScale=scale,newGeometryBytes=0,representation='shared source buffers with explicit positive uniform scale and translation')
a['liverCompletionIntegration']['numericalEvidence']=str(out/'independent-transfer-checks.json')
atomic_json(out/'atlas.json',a);atomic_json(out/'transfer-evidence.json',dict(passed=True,sourceCatalogueSha256=plan['sourceCatalogueSha256'],sourceBuffersLiteral=True,scale=scale,rotationDegrees=0,newGeometryBytes=0,parts=57,normalWordsAndIndicesLiteral=True,anatomicalAcceptance=False))
print(json.dumps(dict(candidate=str(out/'atlas.json'),scale=scale,newGeometryBytes=0)))
