"""Audit emitted covering cells against literal fixed bone and pial tissue."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts

ap=argparse.ArgumentParser();ap.add_argument('bank');args=ap.parse_args();start=time.monotonic();bank=Path(args.bank).resolve()
data=json.loads((bank/'tissue.json').read_text());assert data['passed'];raw=(bank/'cells.bin').read_bytes();origin=np.array([0.,1.35,-.05]);frame=data['frame']
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';atlas=json.loads(path.read_text());reader=AtlasReader();meshes={}
for name,p in data['parts'].items():
    v=np.ndarray((p['vertexCount'],3),'<f4',buffer=raw,offset=p['positions']).astype(float)+origin
    f=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices'])
    meshes[name]=Mesh(v,f,frame)
rows=[]
for p in atlas['parts']:
    pid=p['id']
    if p.get('suppressed')or not (pid.startswith('LOCAL:female:')and p['system']=='skeletal'and('vertebra'in pid or pid.startswith('LOCAL:female:sacrum:'))or pid=='TA98CORD:female:pia'):continue
    fixed=reader.read_part(atlas,p,frame=frame,catalogue_path=path)
    for name,mesh in meshes.items():
        result,summary=pair_contacts(mesh,fixed,0.,max_pairs=500000,max_tests=3000000)
        assert not summary['truncated'],'Fixed-context screen exceeded its declared pair budget'
        bad=np.isin(result['relation'],['crossing','coplanar-overlap']);pairs=result['pairs'][bad]
        row=dict(covering=name,part_id=pid,fixed_geometry_sha256=fixed.geometry_sha,crossings=int(bad.sum()),examples=pairs[:12].tolist(),summary=summary)
        rows.append(row);print(json.dumps(dict(covering=name,fixed=p['name'],crossings=int(bad.sum()),elapsed_seconds=time.monotonic()-start)),flush=True)
        atomic_json(bank/'fixed-context-screen.json',dict(bank_sha256=digest_file(bank/'cells.bin'),checks=rows,complete=False,skeleton_changed=False,anatomical_acceptance=False),overwrite=True)
proof=dict(bank_sha256=digest_file(bank/'cells.bin'),checks=rows,complete=True,passed=not any(r['crossings']for r in rows),
    scope='Native female skeleton and installed pia; unregistered HRA donor bones are separate sources, not the placement frame',
    excluded_donor_sources=[dict(part_id=p['id'],status=p.get('reviewStatus'))for p in atlas['parts']if p['id'].startswith('HRA:')and('vertebra'in p['id']or 'sacrum'in p['id'])],
    skeleton_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start)
atomic_json(bank/'fixed-context-screen.json',proof,overwrite=True);assert proof['passed'],'Covering crosses a fixed source bone or the finite pia'
