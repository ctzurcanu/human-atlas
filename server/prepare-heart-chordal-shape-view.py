"""Source-preserving diagnostic colours for intrinsic leaflet/strand shape."""
from pathlib import Path
import copy,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
h=ROOT/'.local-models/geometry-workbench/heart';src=h/'apparatus-source-v203';fields=h/'chordal-shape-fields-v204';out=h/'chordal-shape-view-v205';out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};info=json.loads((src/'source.json').read_text());reader=AtlasReader();data=bytearray();ci=len(a['chunks']);rows=[]
def put(values):
    while len(data)%4:data.append(0)
    offset=len(data);data.extend(np.ascontiguousarray(values).tobytes());return offset
for index in [1,5,8,9,10]:
    original=Mesh.load(src/f'part-{index}.npz');pid=info['parts'][index]['part_id'];p=parts[pid];record=copy.deepcopy(p);bank=reader._bank(reader.asset(a['chunks'][p['chunk']]['url']));used=np.array(original.source['source_vertex_ids'],int);n=p['vertexCount'];z=np.load(fields/f'part-{index}-fields.npz');ratio=z['transverse_to_longitudinal_ratio'][1][z['source_vertex_map']];known=np.isfinite(ratio);t=np.clip((.35-np.nan_to_num(ratio,nan=.35))/.28,0,1);blue=np.array([.20,.52,.98]);orange=np.array([1.,.43,.10]);colors=(1-t[:,None])*blue+t[:,None]*orange;colors[~known]=[.85,.2,.85]
    p['positions']=put(reader.view(bank,record['positions'],n,'<f4',3)[used]);p['indices']=put(original.faces.astype('<u4'))
    for name,dtype,width in [('normals','<i2',3),('normalSigns','i1',1),('uvs','<f4',2)]:
        if name in record:p[name]=put(reader.view(bank,record[name],n,dtype,width)[used])
    p['colors']=put(np.rint(colors*255).astype('u1'));p.update(chunk=ci,vertexCount=len(used),material='HEART:intrinsic-shape-diagnostic')
    rows.append(dict(part_id=pid,source_geometry_sha256=original.geometry_sha,unmeasured_controls=int((~known).sum()),source_exterior_unchanged=True))
bank_path=out/'shape-view.bin';atomic_bytes(bank_path,bytes(data));a['chunks'].append(dict(url='/local-models/'+str(bank_path.relative_to(ROOT/'.local-models')),bytes=len(data)));a['materials']['HEART:intrinsic-shape-diagnostic']=dict(color=[255,255,255],vertexColors=True,roughness=.8,metalness=0)
candidate=out/'atlas.json';atomic_json(candidate,a)
for row in rows:
    actual=reader.read_part(a,parts[row['part_id']],frame='model:local-ta98:assembled-metres',catalogue_path=candidate);source=Mesh.load(src/f"part-{next(r['i']for r in info['parts']if r['part_id']==row['part_id'])}.npz");assert actual.geometry_sha==source.geometry_sha
    for key in ['normals','normalSigns','uvs']:
        if key in source.attributes:assert np.array_equal(actual.attributes[key],source.attributes[key]),(row['part_id'],key)
atomic_json(out/'diagnostic.json',dict(source_catalogue_sha256=digest_file(path),geometry_modified=False,parts=rows,new_bytes=len(data),anatomical_partition_established=False,installed=False))
print('Source-identical shape diagnostic:',len(rows),'AV sources;',len(data),'bytes',flush=True)
