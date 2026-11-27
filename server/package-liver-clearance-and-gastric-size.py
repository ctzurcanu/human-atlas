"""Apply only passed common fields and rebind literal runtime evidence."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.packed_field_transport import transport
from geometry_tools.liver_impression_field import posterior_impression
from geometry_tools.gastric_body_field import expand_body
from geometry_tools.liver_completion import validate_liver_completion
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);args=ap.parse_args();sex=args.sex
b=ROOT/'.local-models/geometry-workbench/liver-completion';source=b/'final-ready-v265'/(sex+'.json');out=b/'final-clearance-v270'/sex;out.mkdir(parents=True,exist_ok=False);start=time.monotonic();a=json.loads(source.read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader()
probe=json.loads((b/'final-ready-v265'/(sex+'-port-safe-v267')/'fit.json').read_text());assert probe['passed'];fit=probe['chosen'];S=a['liverTransferIntegration']['scale']if sex=='female'else 1.;T=np.array(a['liverTransferIntegration']['translation_m'])if sex=='female'else np.zeros(3)
def field(p):return posterior_impression(p,fit['amplitude_m'],down=fit.get('down_m',0.))
ids=set(a['liverIntegration']['partIds']+a['liverCompletionIntegration']['partIds']+a['hepaticConnectionIntegration']['wholePartIds']+a['hepaticConnectionIntegration']['regionPartIds']);ids={pid for pid in ids if not parts[pid].get('suppressed')}
rewrites,records=transport(a,ids,field,out/'hepatic',source_scale=S,source_translation=T)
# Closed render classifiers share some position blocks with real tissue and
# some computational-only closures. Transport all of them from the same source.
for marker in [a['liverSectionIntegration'],a['liverCompletionIntegration']['functionalSections']]:
 oldurl=marker['bankUrl'];entry=rewrites[oldurl];old=r.asset(oldurl).read_bytes();path=Path(entry['path']);data=bytearray(path.read_bytes())
 for desc in marker.get('volumes',[])+marker['exclusions']:
  n=desc['vertexCount'];scale=float(desc.get('geometryScale',1.));origin=np.array(desc['origin']);raw=r.view(old,desc['positions'],n,'<f4',3);src=(raw.astype(float)*scale+origin-T)/S;v,J=field(src);local=raw.copy();changed=np.any(v!=src,axis=1);local[changed]=((v[changed]*S+T-origin)/scale).astype('<f4');data[desc['positions']:desc['positions']+local.nbytes]=local.tobytes();f=r.view(old,desc['indices'],desc['indexCount'],'<u4').reshape(-1,3);desc['geometrySha256']=Mesh(local.astype(float)*scale+origin,f,'model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres').geometry_sha
 atomic_bytes(path,bytes(data),overwrite=True);entry['sha256']=digest_file(path)
for marker in [a['liverSectionIntegration'],a['liverCompletionIntegration']['functionalSections']]:
 entry=rewrites[marker['bankUrl']];marker.update(bankUrl=entry['url'],bankSha256=entry['sha256'])
host=parts[a['liverFeatureIntegration']['hostPartId']];hosturl=a['chunks'][host['chunk']]['url'];hostsha=digest_file(r.asset(hosturl));a['liverIntegration'].update(bankUrl=hosturl,bankSha256=hostsha);a['liverFeatureIntegration']['bankSha256']=hostsha;a['liverCompletionIntegration']['hostBankSha256']=hostsha;a['liverCompletionIntegration']['functionalSections']['sourceBankSha256']=hostsha
for role,bank in a['hepaticConnectionIntegration']['banks'].items():
 if bank['url']in rewrites:bank.update(url=rewrites[bank['url']]['url'],sha256=rewrites[bank['url']]['sha256'])
cache={}
def rebind(value):
 if isinstance(value,dict):
  if {'partId','face','barycentric','point'}<=value.keys():
   pid=value['partId']
   if pid not in cache:cache[pid]=r.read_part(a,parts[pid],frame='rebinding:assembled-metres')
   value['point']=(np.array(value['barycentric'])@cache[pid].triangles[value['face']]).tolist()
  for child in value.values():rebind(child)
 elif isinstance(value,list):
  for child in value:rebind(child)
for marker in [a['liverFeatureIntegration'],a['liverCompletionIntegration']]:rebind(marker['termGroups'])
for g in a['liverFeatureIntegration']['termGroups']:g['surface']['sourceBankSha256']=hostsha
# Independently preserve actual receiving mouth positions and signed normals.
before=json.loads(source.read_text());bp={p['id']:p for p in before['parts']};seams=[]
for interface in a['liverFinalIntegration']['interfaces']:
 for side in ['first','second']:
  pid=interface[side];m=r.read_part(a,parts[pid],frame='interface:assembled-metres');old=r.read_part(before,bp[pid],frame=m.frame);ff=interface[side+'_faces'];assert np.array_equal(m.triangles[ff],old.triangles[ff]),('Native mouth moved',pid)
 seams.append(dict(first=interface['first'],second=interface['second'],fixed_literal=True))
if sex=='male':
 gp=json.loads((b/'male-gastric-volume-v268/fit.json').read_text());assert gp['passed'];strength=gp['chosen']['strength'];gids=a['gastricCompleteWall']['partIds'];gmap,grecords=transport(a,gids,lambda p:expand_body(p,strength),out/'gastric')
 for pid in gids:
  if parts[pid]['regionTerm']in ['A05.5.01.007','A05.5.01.014']:
   old=r.read_part(before,bp[pid],frame='gastric:assembled-metres');new=r.read_part(a,parts[pid],frame=old.frame);assert np.array_equal(old.vertices,new.vertices)and old.attributes_sha==new.attributes_sha,('Gastric port region changed',pid)
 a['gastricCompleteWall']['sizeCorrection']=dict(schema='gastric-body-size-field/v1',strength=strength,field='geometry_tools/gastric_body_field.py',outerVolumeBeforeMl=gp['trials'][0]['volume_ml'],outerVolumeAfterMl=gp['chosen']['volume_ml'],cardiaAndPylorusFixed=True,physiologicalCalibration=False)
 if len(gmap)==1:a['gastricCompleteWall']['payloadSha256']=next(iter(gmap.values()))['sha256']
 records+=grecords
a['liverFinalIntegration']['clearanceField']=dict(source=str(source),probe=str(b/'final-ready-v265'/(sex+'-port-safe-v267')/'fit.json'),parameters=fit,positiveCommonField=True,nativeMouthsFixed=True,oesophagusFixed=True,diaphragmFixed=True)
witnesses=validate_liver_completion(a,reader=AtlasReader(),frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres');atomic_json(out/'atlas.json',a);atomic_json(out/'field-evidence.json',dict(passed=True,records=records,reference_groups=len(witnesses),fixed_annuli=seams,elapsed_seconds=time.monotonic()-start));print(sex,'common field packed',len(ids),'hepatic parts;',len(witnesses),'reference groups replay',flush=True)
