"""Prepare a local heart update using one male bank and a positive female instance."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True)
ap.add_argument('--bank',required=True);ap.add_argument('--stage',required=True)
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart'
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/args.sex/'atlas.json'
a=json.loads(source.read_text());parts={p['id']:p for p in a['parts']};before=copy.deepcopy(parts)
marker=a['finiteCardiacIntegration'];e=json.loads((h/args.bank/'emission.json').read_text())
out=h/args.stage;out.mkdir(exist_ok=False);bank=Path(e['bank']);assert digest_file(bank)==e['bank_sha256']
root_map={p.get('cardiacTemplateSourcePartId'):p['id']for p in a['parts']if p.get('cardiacTemplateSourcePartId')}
mapping={}
for name,pid in e['part_ids'].items():
    target=pid
    if args.sex=='female':
        target=pid.replace('TA98HEART:male:','TA98HEART:female:')if pid.startswith('TA98HEART:male:')else root_map.get(pid)
    assert target in parts,(pid,target)
    mapping[pid]=target
reference=parts[marker['wallPartIds'][0]]
origin=list(reference['geometryOrigin']);pose=copy.deepcopy(reference.get('geometryInstance'))
ci=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size))
layout={};part_ids={}
for name,record in e['parts'].items():
    male_id=e['part_ids'][name];pid=mapping[male_id];p=parts[pid]
    for key in ['uvs','geometryScale','sourceGeometryOffset','geometryInstance','vertexColors']:
        p.pop(key,None)
    row=copy.deepcopy(record);row['id']=pid;row['geometryOrigin']=origin
    row['tissueDrawGroups']=[dict(g,**({'neighbor':mapping[g['neighbor']]}if g.get('neighbor')else{}))for g in record['tissueDrawGroups']]
    if pose:row['geometryInstance']=pose
    p.update(row,chunk=ci,sectionSolid=True,pairedTissue=True,
             tissueWallAssembly='heart:'+args.sex,defaultHidden=False)
    if name.startswith('root:'):
        old=a.get('materials',{}).get(p.get('material'),{})
        color=old.get('color',[118,83,114]);material='HEART:'+args.sex+':pulmonary-vein-wall'
        appearance=dict(color=[v/255 for v in color]if max(color)>1 else color,
                        pattern='solid',scale=.0005,relief=.000001)
        a['materials'][material]=dict(color=color,roughness=.55,metalness=0,
            pairedNormalSigns=True,tissueAppearance=appearance,geometryOrigin=origin,
            **({'geometryInstance':pose}if pose else{}))
        p.update(material=material,tissueWallLayer='pulmonary-vein-wall',
            tissueAppearance=appearance,
            cardiacRootWall=dict(schema='finite-pv-root-wall/v1',
                bankSha256=e['bank_sha256'],sourceMalePartId=male_id,
                distalRecipientJoinVerified=False,innerRadialFraction=.8,
                histologicalThicknessVerified=False))
    layout[name]=p;part_ids[name]=pid
candidate=out/'atlas.json';atomic_json(candidate,a)
reader=AtlasReader();frame='model:'+('local-ta98'if args.sex=='male'else'local-female-ta98')+':assembled-metres'
for name,p in layout.items():
    mesh=reader.read_part(a,p,frame=frame,catalogue_path=candidate)
    p['bounds']=mesh.bounds.tolist()
marker.update(wallBankSha256=e['bank_sha256'],
    rootWallPartIds=[part_ids[name]for name in part_ids if name.startswith('root:')],
    rootWallEvidence=e['root_interface_evidence'])
a['triangles']=sum(p['indexCount']//3 for p in a['parts']);atomic_json(candidate,a,overwrite=True)
emission=dict(e,frame=frame,origin=origin,parts=layout,part_ids=part_ids,
              catalogue=str(candidate),sourceCatalogueSha256=digest_file(source),
              installed=False)
atomic_json(out/'emission.json',emission)
if pose:atomic_json(out/'transfer.json',dict(positive_similarity=pose,origin=origin,
    source_bank_sha256=e['bank_sha256'],source_normal_words_literal=True,
    source_positions_literal=True,new_geometry_banks=0))
atomic_json(out/'construction.json',dict(prepared=True,sex=args.sex,
    source_catalogue_sha256=digest_file(source),bank_sha256=e['bank_sha256'],
    changed_parts=[pid for pid,p in parts.items()if p!=before[pid]],
    new_parts=0,new_textures=0,shared_bank_bytes=bank.stat().st_size,
    elapsed_seconds=time.monotonic()-start,installed=False))
print(args.sex,'catalogue prepared;',len(layout),'parts;',bank.stat().st_size,'shared bytes',flush=True)
