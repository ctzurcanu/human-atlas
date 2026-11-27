"""Transfer the coordinated valve apparatus and native root references together."""
from pathlib import Path
import argparse,copy,json,re,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
ap=argparse.ArgumentParser();ap.add_argument('--wall',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);mp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';fp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';male=json.loads(mp.read_text());a=json.loads(fp.read_text());source={p['id']:p for p in male['parts']};old={p['id']:p for p in a['parts']};mc={c['id']:c for c in male['concepts']};byc={c['id']:c for c in a['concepts']};marker=male['finiteCardiacIntegration'];wall=json.loads((b/args.wall/'emission.json').read_text());origin=np.asarray(wall['origin']);fit=json.loads((b/'female-fit-v68/fit.json').read_text());r=np.asarray(fit['rotation']);s=fit['scale'];translation=np.asarray(fit['translation_m']);reader=AtlasReader();mapping={pid:pid.replace('TA98HEART:male:','TA98HEART:female:')for pid in marker['wallPartIds']};auxids=marker['valvePartIds']+marker['papillaryPartIds']+[pid for pid in marker['rootPartIds']if pid!='ZA:Ascending aorta'];mapping.update({pid:'LOCAL:female:heart-template:'+re.sub('[^a-z0-9]+','-',pid.lower()).strip('-')for pid in auxids});mapping.update({'ZA:Right atrium':'LOCAL:female:r_atrium','ZA:Left atrium':'LOCAL:female:l_atrium','ZA:Right ventricle':'TA98HEART:female:right-ventricle','ZA:Left ventricle':'TA98HEART:female:left-ventricle','ZA:Ascending aorta':'LOCAL:female:aortic_arch:1513'});buffer=bytearray();aux=[];proofs=[]
def put(value):
 buffer.extend(b'\0'*((-len(buffer))%4));at=len(buffer);buffer.extend(value.tobytes());return at
for pid in auxids:
 p=copy.deepcopy(source[pid]);m=reader.read_part(male,p,frame='model:local-ta98:assembled-metres',catalogue_path=mp);v=((s*(m.vertices@r.T)+translation-origin).astype('<f4').astype(float)+origin);f=m.faces.copy();effective=m.attributes['normals'].astype(float)*m.attributes.get('normalSigns',np.ones(len(m.vertices)))[:,None];n=effective@r.T;n/=np.linalg.norm(n,axis=1)[:,None];words=np.rint(n*32767).astype('<i2');q=v[f];fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);bad=np.argwhere(np.einsum('fki,fi->fk',words[f].astype(float),fn)<=0)
 if len(bad):
  oldids=f[bad[:,0],bad[:,1]];extra=fn[bad[:,0]];extra/=np.linalg.norm(extra,axis=1)[:,None];v=np.r_[v,v[oldids]];words=np.r_[words,np.rint(extra*32767).astype('<i2')];f[bad[:,0],bad[:,1]]=np.arange(len(v)-len(oldids),len(v))
 mesh=Mesh(v,f,'model:local-female-ta98:assembled-metres',dict(normals=words,normalSigns=np.ones(len(v),'i1')));proof=audit(mesh,require_closed=False,max_pairs=1300000,max_tests=14000000);atomic_json(out/(str(len(aux))+'-audit.json'),proof);assert proof['passed'],(pid,proof['failures']);proofs.append(dict(sourcePartId=pid,partId=mapping[pid],geometrySha256=mesh.geometry_sha,passed=True,closed=proof['topology']['closed_oriented_manifold']));print('Female apparatus',len(aux)+1,len(auxids),pid,'passed',flush=True)
 for key in ['positions','normals','indices','normalSigns','uvs','colors','sourceOffset','geometryOrigin','geometryScale','tissueDrawGroups','gastricLumen','tissueChart','material']:p.pop(key,None)
 label=source[pid]['name'];kind='root'if pid in marker['rootPartIds']else'valve'if pid in marker['valvePartIds']else'papillary';color=[126,108,176]if p['system']=='venous'else[177,72,73]if p['system']=='arterial'else[211,171,153]if kind=='valve'else[161,77,70]if kind=='papillary'else[185,87,88];material='HEART:female:'+kind+':'+p['system'];a['materials'][material]=dict(color=color,roughness=.67,metalness=0,pairedNormalSigns=True,geometryOrigin=origin.tolist());p.update(id=mapping[pid],conceptId=mapping[pid],sourceId='Heart:transferred:'+pid,positions=put((v-origin).astype('<f4')),normals=put(words),indices=put(f.astype('<u4')),normalSigns=put(np.ones(len(v),'i1')),vertexCount=len(v),indexCount=f.size,triangles=len(f),bounds=mesh.bounds.tolist(),geometryOrigin=origin.tolist(),material=material,pairedTissue=kind!='root',sectionSolid=kind!='root'and proof['topology']['closed_oriented_manifold'],tissueWallAssembly='heart:female',tissueWallLayer=kind,defaultHidden=False,cardiacTemplateSourcePartId=pid,cardiacTemplateSourceCatalogueSha256=digest_file(mp));p.pop('suppressed',None);p.pop('replacementParts',None);aux.append(p)
bank=out/'apparatus.bin';atomic_bytes(bank,bytes(buffer));ci=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=len(buffer)))
for p in aux:p['chunk']=ci
a['parts'].extend(aux);wi=len(a['chunks']);wb=Path(wall['bank']);a['chunks'].append(dict(url='/local-models/'+str(wb.relative_to(ROOT/'.local-models')),bytes=wb.stat().st_size))
for name,packed in wall['parts'].items():
 oldpid='TA98HEART:male:'+name;p=copy.deepcopy(source[oldpid]);p.update(packed);p.update(id=mapping[oldpid],conceptId=mapping[oldpid],chunk=wi,tissueWallAssembly='heart:female');a['parts'].append(p)
for key,spec in male['materials'].items():
 if key.startswith('HEART:'):a['materials'][key]=dict(copy.deepcopy(spec),geometryOrigin=origin.tolist())
def rename(value):
 if isinstance(value,str):return mapping.get(value,value.replace('TA98HEART:male:','TA98HEART:female:').replace('heart:male','heart:female'))
 if isinstance(value,list):return[rename(v)for v in value]
 if isinstance(value,dict):return{k:rename(v)for k,v in value.items()}
 return value
newmarker=rename(copy.deepcopy(marker));newmarker.update(sex='female',sourceCatalogueSha256=digest_file(fp),templateCatalogueSha256=digest_file(mp),wallBankSha256=wall['bank_sha256'],transfer=fit,recipientDistalVascularConnections='Native aortic arch retained. Matching of transferred venous/pulmonary root ends to recipient distal trees remains provisional.');a['finiteCardiacIntegration']=newmarker
selected=set(marker['wallPartIds']+auxids);conceptids={source[pid]['conceptId']for pid in selected}|{g['conceptId']for g in marker['termGroups']};newconcepts=[rename(copy.deepcopy(c))for c in male['concepts']if c['id']in conceptids];newids={c['id']for c in newconcepts};a['concepts']=[c for c in a['concepts']if c['id']not in newids]+newconcepts
replacement={
 'LOCAL:female:r_atrium:1512':[mapping['TA98HEART:male:myocardium:right-atrium']],
 'LOCAL:female:l_atrium:1507':[mapping['TA98HEART:male:myocardium:left-atrium']],
 'LOCAL:female:ventricles:1506':[mapping[pid]for pid in marker['wallPartIds']if any(x in pid for x in ['right-ventricle','left-ventricle','interventricular'])],
 'LOCAL:female:superior_vena_cava:1510':[mapping['ZA:Superior vena cava']],
 'LOCAL:female:inferior_vena_cava:1511':[mapping['ZA:Inferior vena cava (thoracic part)']],
 'LOCAL:female:pulmonary_trunk:1509':[mapping['ZA:Pulmonary trunk']],
}
for oldid,key in [('LOCAL:female:aortic_semilunar_valve:1504','aortic-valve'),('LOCAL:female:pulmonary_semilunar_valve:1505','pulmonary-valve'),('LOCAL:female:tricuspid_atrioventricular_valve:1503','tricuspid-valve'),('LOCAL:female:l_bicuspid_atrioventricular_valve:1502','mitral-valve')]:replacement[oldid]=next(g['partIds']for g in newmarker['termGroups']if g['key']==key)
for pid,ids in replacement.items():old[pid].update(suppressed=True,defaultHidden=True,replacementParts=ids)
for c in a['concepts']:c['elements']=list(dict.fromkeys(pid for oldid in c.get('elements',[])for pid in replacement.get(oldid,[oldid])))
active=set(newmarker['wallPartIds']+newmarker['valvePartIds']+newmarker['papillaryPartIds']+newmarker['rootPartIds']);a['defaultHidden']=sorted((set(a.get('defaultHidden',[]))|set(replacement))-active);a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts']);atomic_json(out/'atlas.json',a);atomic_json(out/'transfer.json',dict(passed=all(p['passed']for p in proofs),source_catalogue_sha256=digest_file(mp),target_catalogue_sha256=digest_file(fp),wall_bank_sha256=wall['bank_sha256'],apparatus_bank_sha256=digest_file(bank),parts=proofs,fit=fit,all_tissues_use_one_similarity=True,recipient_skeleton_lungs_diaphragm_fixed=True,retired_native_parts=list(replacement),remaining_vascular_limits=newmarker['recipientDistalVascularConnections'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False));print('Female heart catalogue prepared',len(active),'parts',flush=True)
