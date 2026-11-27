"""Transfer one complete reviewed liver assembly with shared packed arrays."""
from pathlib import Path
import argparse,copy,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_bytes,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh

ap=argparse.ArgumentParser();ap.add_argument('--version',default='v202');ap.add_argument('--fit');args=ap.parse_args()
inputs=ROOT/'.local-models/geometry-workbench/liver-completion/female-transfer-v202'
base=inputs.parent/('female-transfer-'+args.version)
if base!=inputs:base.mkdir(exist_ok=False)
plan=json.loads((inputs/'plan.json').read_text());fit=json.loads(Path(args.fit or inputs/'fit.json').read_text())
mp=Path(plan['sourceCatalogue']);fp=Path(plan['targetCatalogue'])
assert digest_file(mp)==plan['sourceCatalogueSha256']and digest_file(fp)==plan['targetCatalogueSha256']
m=json.loads(mp.read_text());a=json.loads(fp.read_text());reader=AtlasReader()
scale=fit['scale'];translation=np.array(fit['translation_m']);origin=scale*np.array([0.,1.25,0.])+translation
source={p['id']:p for p in m['parts']};ids=plan['movingPartIds'];selected=set(ids)
idmap={'ZA:Liver':'LOCAL:female:liver','TA98PART:FJ3079:A05.8.01.061':'TA98HEPATIC:female:common-hepatic-parent','ZA:Inferior vena cava (abdominal part)':'LOCAL:female:vena_cava:1797'}
def rename(s):
    if s in idmap:return idmap[s]
    if s=='local-ta98':return 'local-female-ta98'
    for x,y in [('TA98HEPATIC:male:','TA98HEPATIC:female:'),('TA98LIVER:male:','TA98LIVER:female:'),('HEPATIC:male:','HEPATIC:female:'),('hepatic-connections:male:','hepatic-connections:female:'),('liver:male:','liver:female:'),('male-liver-','female-liver-'),('male-hepatic-','female-hepatic-')]:s=s.replace(x,y)
    return s
def remap(x):
    if isinstance(x,str):return rename(x)
    if isinstance(x,list):return [remap(v)for v in x]
    if isinstance(x,dict):return {k:remap(v)for k,v in x.items()}
    return x
transferred=[remap(copy.deepcopy(source[pid]))for pid in ids];byid={p['id']:p for p in transferred}
class Bank:
    def __init__(self,ci):self.ci=ci;self.raw=reader.asset(m['chunks'][ci]['url']).read_bytes();self.data=bytearray();self.cache={}
    def put(self,offset,length,dtype,width=1,oldorigin=None):
        key=(offset,length,dtype,width,tuple(oldorigin)if oldorigin is not None else None)
        if key in self.cache:return self.cache[key]
        values=np.ndarray((length,width)if width>1 else(length,),dtype,buffer=self.raw,offset=offset)
        if oldorigin is not None:values=((values.astype(float)+oldorigin)*scale+translation-origin).astype('<f4')
        self.data.extend(b'\0'*((-len(self.data))%4));at=len(self.data);self.data.extend(values.tobytes());self.cache[key]=at;return at
banks={}
def bank(ci):
    if ci not in banks:banks[ci]=Bank(ci)
    return banks[ci]
for oldpid,p in zip(ids,transferred):
    old=source[oldpid];b=bank(old['chunk']);n=old['vertexCount'];oldorigin=np.array(old.get('geometryOrigin',[0.,0.,0.]))+np.array(old.get('sourceOffset',[0.,0.,0.]))
    p['positions']=b.put(old['positions'],n,'<f4',3,oldorigin);p['indices']=b.put(old['indices'],old['indexCount'],'<u4')
    for key,dtype,width in [('normals','<i2',3),('normalSigns','i1',1),('uvs','<f4',2),('colors','u1',3)]:
        if key in old:p[key]=b.put(old[key],n,dtype,width)
    p['geometryOrigin']=origin.tolist();p.pop('sourceOffset',None);p['bounds']=(np.array(old['bounds'])*scale+translation).tolist()
    p['liverTransferSourcePartId']=oldpid;p['liverTransferSourceCatalogueSha256']=plan['sourceCatalogueSha256'];p['liverTargetModel']='local-female-ta98'
    if 'hepaticSourceModel'in p:p['hepaticSourceModel']='local-female-ta98'

markers=['liverIntegration','liverFeatureIntegration','liverCompletionIntegration','liverSectionIntegration','hepaticConnectionIntegration','hepaticContainmentIntegration']
for key in markers:a[key]=remap(copy.deepcopy(m[key]))
def transform_volumes(oldmarker,newmarker):
    ci=next(i for i,c in enumerate(m['chunks'])if c['url']==oldmarker['bankUrl']);b=bank(ci)
    for key in ['volumes','exclusions']:
        for old,v in zip(oldmarker.get(key,[]),newmarker.get(key,[])):
            v['positions']=b.put(old['positions'],old['vertexCount'],'<f4',3,np.array(old['origin']))
            v['indices']=b.put(old['indices'],old['indexCount'],'<u4');v['origin']=origin.tolist()
            vv=np.ndarray((v['vertexCount'],3),'<f4',buffer=b.data,offset=v['positions']).astype(float)+origin
            ff=np.ndarray((v['indexCount']//3,3),'<u4',buffer=b.data,offset=v['indices']).copy()
            v['geometrySha256']=Mesh(vv,ff,'model:local-female-ta98:assembled-metres').geometry_sha
    return ci
body_ci=transform_volumes(m['liverSectionIntegration'],a['liverSectionIntegration'])
functional_ci=transform_volumes(m['liverCompletionIntegration']['functionalSections'],a['liverCompletionIntegration']['functionalSections'])
packed={}
for ci,b in banks.items():
    path=base/f'bank-{ci}.bin';atomic_bytes(path,bytes(b.data));url='/local-models/'+str(path.relative_to(ROOT/'.local-models'));newci=len(a['chunks']);a['chunks'].append(dict(url=url,bytes=len(b.data)));packed[ci]=dict(chunk=newci,url=url,sha256=digest_file(path),bytes=len(b.data))
for oldpid,p in zip(ids,transferred):p['chunk']=packed[source[oldpid]['chunk']]['chunk']
body=packed[body_ci];func=packed[functional_ci]
for key in ['liverIntegration','liverSectionIntegration','liverFeatureIntegration']:a[key].update(bankUrl=body['url'],bankSha256=body['sha256'])
a['liverCompletionIntegration']['hostBankSha256']=body['sha256']
a['liverCompletionIntegration']['functionalSections'].update(bankUrl=func['url'],bankSha256=func['sha256'],sourceBankSha256=body['sha256'])
for g in a['liverFeatureIntegration']['termGroups']:g['surface']['sourceBankSha256']=body['sha256']
host=byid['TA98PART:ZA:Liver:capsule'];hbank=banks[source[host['id']]['chunk']]
hv=np.ndarray((host['vertexCount'],3),'<f4',buffer=hbank.data,offset=host['positions']).astype(float)+origin
hf=np.ndarray((host['indexCount']//3,3),'<u4',buffer=hbank.data,offset=host['indices'])
def bind(x):
    if isinstance(x,dict):
        if x.get('partId')==host['id']and 'face'in x and 'barycentric'in x and 'point'in x:x['point']=(np.array(x['barycentric'])@hv[hf[x['face']]]).tolist()
        for v in x.values():bind(v)
    elif isinstance(x,list):
        for v in x:bind(v)
bind(a['liverFeatureIntegration']);bind(a['liverCompletionIntegration'])
for role in ['artery','portal','hepatic_vein','biliary']:
    p=byid['TA98HEPATIC:female:'+role];old=source['TA98HEPATIC:male:'+role];b=banks[old['chunk']]
    signature=hashlib.sha256(b''.join(b.data[p[k]:p[k]+p['vertexCount']*n]for k,n in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps(origin.tolist()).encode()).hexdigest()
    for part in transferred:
        if part.get('hepaticConnectionRole')==role:part['hepaticAttributeSignature']=signature
    a['hepaticConnectionIntegration']['banks'][role]=dict(url=packed[old['chunk']]['url'],sha256=packed[old['chunk']]['sha256'])

oldconcepts={c['id']:c for c in a['concepts']};retired=[]
for p in a['parts']:
    term=oldconcepts.get(p['conceptId'],{}).get('ta98Term','')
    oldbranch=p['id']in ['LOCAL:female:l_hepatic_artery:1717','LOCAL:female:r_hepatic_artery:1719','LOCAL:female:intermediate_hepatic_artery:1716','LOCAL:female:accessory_left_hepatic_artery:1715','LOCAL:female:hepatic_veins:1794']
    if (term.startswith('A05.8.01.')and term not in ['A05.8.01.061'])or p['id']=='HRA:VH_F_common_hepatic_duct'or oldbranch:
        p['suppressed']=True;p['defaultHidden']=True;retired.append(p['id'])
a['liverIntegration']['sourceParts']=retired
needed={source[pid]['conceptId']for pid in ids}
for key in ['liverIntegration','liverFeatureIntegration','liverCompletionIntegration','hepaticConnectionIntegration']:needed.update(g['conceptId']for g in m[key]['termGroups'])
concepts=[remap(copy.deepcopy(c))for c in m['concepts']if c['id']in needed or c.get('elements')and set(c['elements'])<=selected]
available=set(byid)|{p['id']for p in a['parts']}
for c in concepts:
    c['elements']=[pid for pid in c.get('elements',[])if pid in available]
    if c['id']=='LOCAL:female:liver':c['name']='Liver'
newcs={c['id']for c in concepts};a['concepts']=[c for c in a['concepts']if c['id']not in newcs]+concepts
a['parts']+=transferred;a['defaultHidden']=list(dict.fromkeys(a.get('defaultHidden',[])+retired+[p['id']for p in transferred if p.get('defaultHidden')]))
a['liverTransferIntegration']=dict(schema='male-to-female-liver/v1',sourceCatalogueSha256=plan['sourceCatalogueSha256'],targetPreviousCatalogueSha256=plan['targetCatalogueSha256'],scale=scale,translation_m=translation.tolist(),rotationDegrees=0,partIds=[p['id']for p in transferred],retiredPartIds=retired,fixedSkeleton=True,allPartsUseOneTransform=True,targetDistalJunctions='Native distal female vascular/biliary continuity is a separate check; proximity does not establish a joined wall.',localOnly=True)
a['liverCompletionIntegration']['numericalEvidence']=str(base/'transfer-evidence.json')
atomic_json(base/'atlas.json',a)
atomic_json(base/'transfer-evidence.json',dict(sourceCatalogueSha256=plan['sourceCatalogueSha256'],targetPreviousCatalogueSha256=plan['targetCatalogueSha256'],scale=scale,orientationUnchanged=True,partCount=len(transferred),indicesAndNormalWordsLiteral=True,packedSharedArrays=True,packedBytes=sum(x['bytes']for x in packed.values()),retiredParts=retired,allUnrelatedRecordsLiteral=True,anatomicalAcceptance=False))
print(json.dumps(dict(candidate=str(base/'atlas.json'),scale=scale,parts=len(transferred),retired=len(retired),banks=len(packed),bytes=sum(x['bytes']for x in packed.values()))))
