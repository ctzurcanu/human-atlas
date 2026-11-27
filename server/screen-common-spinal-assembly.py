"""Check the actual fitted cord exterior, pia and outer sac against fixed bone."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts

ap=argparse.ArgumentParser();ap.add_argument('stage');ap.add_argument('--cord');ap.add_argument('--carrier');ap.add_argument('--packed',action='store_true');ap.add_argument('--encoded-stage',action='store_true',help='Screen a source-bound pending bank without granting it installation status');args=ap.parse_args();stage=Path(args.stage).resolve();start=time.monotonic();frame='model:local-female-ta98:assembled-metres'
cord=Path(args.cord).resolve()if args.cord else stage/'cord';carrier=Path(args.carrier).resolve()if args.carrier else stage/'carrier'
if args.packed:
    d=json.loads((cord/('pending.json'if args.encoded_stage else'tissue.json')).read_text());raw=(cord/'cells.bin').read_bytes();selected=[];position_key=None
    assert d['bank_sha256']==digest_file(cord/'cells.bin')
    for name,p in d['parts'].items():
        if name.split('@')[0]!='pia':continue
        key=(p['positions'],p['vertexCount'])
        if key!=position_key:
            positions=np.ndarray((p['vertexCount'],3),'<f4',buffer=raw,offset=p['positions']).astype(float)+[0.,1.35,-.05];position_key=key
        faces=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices'])
        for g in p['tissueDrawGroups']:
            if g.get('neighbor','').startswith('TA98CORD:female:pia@'):continue
            selected.append(faces[g['start']//3:(g['start']+g['count'])//3])
    assert selected,'The packed bank has no pial charts'
    pia=Mesh(positions,np.concatenate(selected),frame);cord_hash=digest_file(cord/'cells.bin')
else:
    d=json.loads((cord/'cells.json').read_text());z=np.load(cord/'cells.npz');uses=d['cells']['pia'];pia_faces=np.concatenate([z[u['surface']]if u['side']==1 else z[u['surface']][:,[0,2,1]]for u in uses]);pia=Mesh(z['positions'],pia_faces,frame);cord_hash=digest_file(cord/'cells.npz')
q=np.load(carrier/'surface.npz');sac=Mesh(q['vertices'],q['faces'],frame)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';atlas=json.loads(path.read_text());reader=AtlasReader();checks=[]
for p in atlas['parts']:
    if p.get('suppressed')or not(p['id'].startswith('LOCAL:female:')and p['system']=='skeletal'and('vertebra'in p['id']or ':sacrum:'in p['id'])):continue
    bone=reader.read_part(atlas,p,frame=frame,catalogue_path=path)
    for name,m in [('pia',pia),('outer-dura',sac)]:
        r,s=pair_contacts(m,bone,0.,max_pairs=200000,max_tests=2000000);assert not s['truncated']
        wrong=np.isin(r['relation'],['crossing','coplanar-overlap']);pairs=r['pairs'][wrong]
        checks.append(dict(tissue=name,fixed_part=p['id'],fixed_geometry_sha256=bone.geometry_sha,crossings=int(wrong.sum()),examples=pairs[:12].tolist(),complete=True))
        if wrong.any():print(json.dumps(checks[-1]),flush=True)
r,s=pair_contacts(sac,pia,0.,max_pairs=200000,max_tests=2000000);assert not s['truncated'];bad=np.isin(r['relation'],['crossing','coplanar-overlap'])
checks.append(dict(tissue='outer-dura',fixed_part='fitted-pia',crossings=int(bad.sum()),examples=r['pairs'][bad][:12].tolist(),complete=True))
proof=dict(source_cell_sha256=cord_hash,actual_packed_float32=bool(args.packed),source_carrier_sha256=digest_file(carrier/'surface.npz'),checks=checks,
    passed=not any(x['crossings']for x in checks),skeleton_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start)
atomic_json(stage/'fixed-context-screen.json',proof);print(json.dumps(dict(passed=proof['passed'],checks=len(checks),crossings=sum(x['crossings']for x in checks),elapsed_seconds=time.monotonic()-start)),flush=True)
assert proof['passed'],'Uninstalled coupled candidate retains actual tissue crossings'
