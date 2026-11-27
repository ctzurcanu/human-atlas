"""Place the certified heart banks once; never requantize their tissue vertices."""
from pathlib import Path
import argparse,copy,json,re,sys,time
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import audit
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='female-instance-v112');args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);mp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';fp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';male=json.loads(mp.read_text());a=json.loads(fp.read_text());source={p['id']:p for p in male['parts']};old={p['id']:p for p in a['parts']};marker=male['finiteCardiacIntegration'];fit=json.loads((b/'female-fit-v68/fit.json').read_text());quaternion=Rotation.from_matrix(fit['rotation']).as_quat().tolist();instance=dict(rotation=quaternion,scale=fit['scale']);r=Rotation.from_quat(quaternion).as_matrix();scale=fit['scale'];translation=np.asarray(fit['translation_m']);reader=AtlasReader();auxids=marker['valvePartIds']+marker['papillaryPartIds']+[pid for pid in marker['rootPartIds']if pid!='ZA:Ascending aorta'];mapping={pid:pid.replace('TA98HEART:male:','TA98HEART:female:')for pid in marker['wallPartIds']};mapping.update({pid:'LOCAL:female:heart-template:'+re.sub('[^a-z0-9]+','-',pid.lower()).strip('-')for pid in auxids});mapping.update({'ZA:Right atrium':'LOCAL:female:r_atrium','ZA:Left atrium':'LOCAL:female:l_atrium','ZA:Right ventricle':'TA98HEART:female:right-ventricle','ZA:Left ventricle':'TA98HEART:female:left-ventricle','ZA:Ascending aorta':'LOCAL:female:aortic_arch:1513'});chunks={};parts=[];proofs=[]
def chunk(ci):
 if ci not in chunks:chunks[ci]=len(a['chunks']);a['chunks'].append(copy.deepcopy(male['chunks'][ci]))
 return chunks[ci]
def rename(value):
 if isinstance(value,str):return mapping.get(value,value.replace('TA98HEART:male:','TA98HEART:female:').replace('heart:male','heart:female'))
 if isinstance(value,list):return[rename(v)for v in value]
 if isinstance(value,dict):return{k:rename(v)for k,v in value.items()}
 return value
for pid in marker['wallPartIds']+auxids:
 original=source[pid];p=rename(copy.deepcopy(original));p.update(id=mapping[pid],conceptId=mapping[pid],chunk=chunk(original['chunk']),geometryInstance=copy.deepcopy(instance),geometryOrigin=(scale*(r@np.asarray(original.get('geometryOrigin',[0,0,0])))+translation).tolist(),defaultHidden=False,cardiacTemplateSourcePartId=pid,cardiacTemplateSourceCatalogueSha256=digest_file(mp));p.pop('suppressed',None);p.pop('replacementParts',None)
 if pid in marker['wallPartIds']:p['tissueWallAssembly']='heart:female'
 selected_component=None
 if 'pulmonary vein'in original['name'].lower():
  m=reader.read_part(male,original,frame='model:local-ta98:assembled-metres',catalogue_path=mp);v,iv=np.unique(m.vertices,axis=0,return_inverse=True);f=iv[m.faces];edges=np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]];count,label=connected_components(coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(v),len(v))).tocsr(),directed=False);candidates=[]
  for i in range(count):
   ids=np.flatnonzero(label[f[:,0]]==i);distance=float(np.linalg.norm(m.vertices[np.unique(m.faces[ids])]-[.018,1.405,-.012],axis=1).min());candidates.append((distance,i,ids))
  _,selected_component,ids=min(candidates,key=lambda q:q[0]);assert np.array_equal(ids,np.arange(ids[0],ids[-1]+1)),'A noncontiguous source component needs an explicit literal index-only bank';p['indices']=original['indices']+int(ids[0])*12;p['indexCount']=len(ids)*3;p['triangles']=len(ids);p['name']='Root of '+original['name'].lower();p['coverageLimitation']='Existing connected proximal pulmonary venous source component. Distal tributaries and recipient vascular joining remain provisional.'
 material='HEART:female:instance:'+str(len(parts));spec=copy.deepcopy(male.get('materials',{}).get(original.get('material'),{}));spec.setdefault('color',[120,112,175]if original['system']=='venous'else[182,73,74]if original['system']=='arterial'else[208,167,151]if pid in marker['valvePartIds']else[160,77,70]);spec.update(geometryOrigin=p['geometryOrigin'],geometryInstance=copy.deepcopy(instance));a['materials'][material]=spec;p['material']=material
 m=reader.read_part(a,p,frame='model:local-female-ta98:assembled-metres');p['bounds']=m.bounds.tolist();parts.append(p)
 if pid not in marker['wallPartIds']:
  proof=audit(m,require_closed=False,max_pairs=1300000,max_tests=14000000);atomic_json(out/(str(len(proofs))+'-audit.json'),proof);assert proof['passed'],(pid,proof['failures']);proofs.append(dict(sourcePartId=pid,partId=p['id'],passed=True,geometrySha256=m.geometry_sha,component=selected_component,geometryInstance=proof.get('geometry_instance')));print('Shared apparatus',len(proofs),len(auxids),pid,'passed',flush=True)
a['parts'].extend(parts);newmarker=rename(copy.deepcopy(marker));newmarker.update(sex='female',sourceCatalogueSha256=digest_file(fp),templateCatalogueSha256=digest_file(mp),transfer=fit,geometryInstance=instance,recipientDistalVascularConnections='Distal recipient vessel joins remain provisional. Native aortic arch retained.');newmarker['remainingLimits']+=' Distal recipient vessel joins remain provisional.';a['finiteCardiacIntegration']=newmarker;selected=set(marker['wallPartIds']+auxids);conceptids={source[pid]['conceptId']for pid in selected}|{g['conceptId']for g in marker['termGroups']};newconcepts=[rename(copy.deepcopy(c))for c in male['concepts']if c['id']in conceptids];newids={c['id']for c in newconcepts};a['concepts']=[c for c in a['concepts']if c['id']not in newids]+newconcepts
for c in a['concepts']:
 if c['id']==newmarker['wholeConceptId']:c['coverageLimitation']=newmarker['remainingLimits']
for p in parts:
 if p['name'].startswith('Root of '):
  c=next((c for c in a['concepts']if c['id']==p['conceptId']),None)
  if c:c.update(name=p['name'],ta98Kind='parent',coverageLimitation=p['coverageLimitation'])
replacement={
 'LOCAL:female:r_atrium:1512':[mapping['TA98HEART:male:myocardium:right-atrium']],
 'LOCAL:female:l_atrium:1507':[mapping['TA98HEART:male:myocardium:left-atrium']],
 'LOCAL:female:ventricles:1506':[mapping[pid]for pid in marker['wallPartIds']if any(q in pid for q in ['right-ventricle','left-ventricle','interventricular'])],
 'LOCAL:female:superior_vena_cava:1510':[mapping['ZA:Superior vena cava']],
 'LOCAL:female:inferior_vena_cava:1511':[mapping['ZA:Inferior vena cava (thoracic part)']],
 'LOCAL:female:pulmonary_trunk:1509':[mapping['ZA:Pulmonary trunk']],
}
for oldid,key in [('LOCAL:female:aortic_semilunar_valve:1504','aortic-valve'),('LOCAL:female:pulmonary_semilunar_valve:1505','pulmonary-valve'),('LOCAL:female:tricuspid_atrioventricular_valve:1503','tricuspid-valve'),('LOCAL:female:l_bicuspid_atrioventricular_valve:1502','mitral-valve')]:replacement[oldid]=next(g['partIds']for g in newmarker['termGroups']if g['key']==key)
for pid,ids in replacement.items():old[pid].update(suppressed=True,defaultHidden=True,replacementParts=ids)
for c in a['concepts']:c['elements']=list(dict.fromkeys(pid for oldid in c.get('elements',[])for pid in replacement.get(oldid,[oldid])))
active=set(newmarker['wallPartIds']+newmarker['valvePartIds']+newmarker['papillaryPartIds']+newmarker['rootPartIds']);a['defaultHidden']=sorted((set(a.get('defaultHidden',[]))|set(replacement))-active);a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts']);catalogue=out/'atlas.json';atomic_json(catalogue,a);e=json.loads((b/'regional-bank-v102/emission.json').read_text());e.update(frame='model:local-female-ta98:assembled-metres',catalogue=str(catalogue),origin=next(p['geometryOrigin']for p in parts if p['id']==newmarker['wallPartIds'][0]),part_ids={name:mapping[pid]for name,pid in e['part_ids'].items()},parts={name:next(p for p in parts if p['id']==mapping['TA98HEART:male:'+name])for name in e['parts']},source_catalogue_sha256=digest_file(fp),geometryInstance=instance);atomic_json(out/'emission.json',e);atomic_json(out/'transfer.json',dict(passed=True,source_catalogue_sha256=digest_file(mp),target_catalogue_sha256=digest_file(fp),parts=proofs,all_tissues_use_one_similarity=True,no_vertex_or_normal_requantization=True,new_geometry_banks=0,recipient_skeleton_lungs_diaphragm_fixed=True,retired_native_parts=list(replacement),remaining_limits=newmarker['remainingLimits'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False));print('Female shared-bank heart prepared',len(active),'parts',flush=True)
