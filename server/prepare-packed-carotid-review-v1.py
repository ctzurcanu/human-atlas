"""Render the actual exported cell bank, including a separately selectable fork."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1/packed-left-v2'
D=B/'actual-bank-render';D.mkdir(exist_ok=False)
audit=json.loads((B/'independent-catalogue-audit.json').read_text())
assert audit['passed']and audit['complete']and not audit['occupied_pair_stage_pending']
ready=json.loads((B/'readiness.json').read_text());byid={p['id']:p for p in ready['parts'].values()}
assert hashlib.sha256((B/'cells.bin').read_bytes()).hexdigest()==audit['bank_sha256']
entries=[];lo=np.array([.012,1.597,-.020]);hi=np.array([.064,1.628,.029])
for pid,proof in audit['parts'].items():
    if proof['role']not in ['CCA','bifurcation','ICA','ECA']:continue
    m=Mesh.load(proof['decoded_path']);assert m.geometry_sha==proof['geometry_sha256']
    assert m.attributes_sha==proof['attributes_sha256']
    tri=m.triangles;keep=np.all(tri.max(1)>=lo,axis=1)&np.all(tri.min(1)<=hi,axis=1)
    groups=[];faces=[];start=0
    for g in byid[pid]['tissueDrawGroups']:
        ff=m.faces[g['start']//3:(g['start']+g['count'])//3]
        kk=keep[g['start']//3:(g['start']+g['count'])//3]
        ff=ff[kk]
        if not len(ff):continue
        groups.append(dict(start=start,count=int(ff.size),**({'neighbor':g['neighbor']}if g.get('neighbor')else{})))
        faces.extend(ff);start+=ff.size
    cropped=Mesh(m.vertices,np.asarray(faces),m.frame,m.attributes,dict(m.source,diagnostic_uncapped_crop=True))
    p=D/(proof['role']+'-actual-bank-crop.npz');cropped.save(p)
    entries.append(dict(mesh=str(p),name='Actual exported '+proof['role'],part_id=pid,tissue_draw_groups=groups,color=[.8,.16,.12,1]))
atomic_json(D/'assembled-spec.json',dict(frame='model:local-ta98:assembled-metres',assembled_tissues=True,meshes=entries,resolution=1100,wall_seconds=40,
    views=[dict(name='actual-bank-assembled-fork',direction=[1,.2,-1]),dict(name='actual-bank-communicating-lumen-cut',direction=[0,1,0],section=dict(normal=[0,1,0],point=[.033,1.6098,0],clear_outer=True))]))
fork=next(p for p in audit['parts'].values()if p['role']=='bifurcation')
atomic_json(D/'bifurcation-spec.json',dict(frame='model:local-ta98:assembled-metres',meshes=[dict(mesh=fork['decoded_path'],name='Actual exported A12.2.04.009 tissue region',color=[.8,.16,.12,1])],resolution=1000,wall_seconds=40,
    views=[dict(name='actual-selectable-bifurcation',direction=[1,.4,-1]),dict(name='actual-bifurcation-open-outlets',direction=[0,1,.15])]))
print(D)
