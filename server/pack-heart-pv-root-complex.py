"""Pack the source-normal cardiac/root charts once with reciprocal normal signs."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,atomic_bytes,digest_file
from geometry_tools.tissue import normal_words
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--interface-stage',required=True)
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart'
source=h/args.source;out=h/args.stage;out.mkdir(exist_ok=False)
meta=json.loads((source/'arrangement.json').read_text());z=np.load(source/'arrangement.npz')
interface_stage=h/args.interface_stage
assert digest_file(interface_stage/'arrangement.npz')==meta['source_arrangement_sha256']
old=h/'regional-bank-v102';e=json.loads((old/'emission.json').read_text());base=np.load(old/'arrangement.npz')
o=np.array(meta['origin']);p=z['positions'];oldp=(base['positions']-o).astype('<f4').astype(float)+o
rv,rf=z['root_positions'],z['root_faces'];root_words=normal_words(rv,{'root':rf})['root']
surfaces={};words={};arrays={'positions':p};surface_arrays={};word_arrays={};normal_evidence=[]
for row in meta['surfaces']:
    sid=row['id'];faces=z[row['array']];parent=z[row['parents']]
    if row.get('parentSurface')=='occupied-outer':
        controls=rv[rf[parent]];sw=root_words[parent].astype(float)
        sign=row.get('parentNormalSign',1)
    else:
        controls=oldp[base[e['surfaces'][sid]][parent]]
        sw=base[e['corner_words'][sid]][parent].astype(float);sign=1
    a,b,c=controls[:,0],controls[:,1],controls[:,2]
    u,v=b-a,c-a;nn=np.cross(u,v);drop=np.argmax(abs(nn),axis=1)
    axes=np.array([[1,2],[0,2],[0,1]])[drop];ii=np.arange(len(faces))[:,None]
    uv=np.take_along_axis(u,axes,axis=1);vv=np.take_along_axis(v,axes,axis=1)
    delta=p[faces]-a[:,None];xy=np.take_along_axis(delta,np.repeat(axes[:,None],3,axis=1),axis=2)
    det=uv[:,0]*vv[:,1]-uv[:,1]*vv[:,0];assert np.all(det!=0),sid
    s=(xy[:,:,0]*vv[:,None,1]-xy[:,:,1]*vv[:,None,0])/det[:,None]
    t=(uv[:,None,0]*xy[:,:,1]-uv[:,None,1]*xy[:,:,0])/det[:,None]
    weights=np.stack([1-s-t,s,t],axis=2)
    field=np.einsum('tck,tkj->tcj',weights,sw)
    field/=np.linalg.norm(field,axis=2)[:,:,None]
    field=np.rint(field*32767).astype('<i2');literal=0
    for corner in range(3):
        for k in range(3):
            exact=np.all(p[faces[:,corner]]==controls[:,k],axis=1)
            field[exact,corner]=sw[exact,k].astype('<i2');literal+=int(exact.sum())
    field=(field.astype('int32')*sign).astype('<i2')
    geo=np.cross(p[faces[:,1]]-p[faces[:,0]],p[faces[:,2]]-p[faces[:,0]])
    length=np.linalg.norm(geo,axis=1);assert np.all(length>0),('Zero stored face',sid)
    assert np.all(np.einsum('ij,ij->i',geo,nn*sign)>0),('Stored facet reversed from its literal source parent',sid)
    geo/=length[:,None]
    n=field.astype(float);n/=np.linalg.norm(n,axis=2)[:,:,None]
    dot=np.einsum('tki,ti->tk',n,geo);bad=dot<=0
    if np.any(bad):
        # A newly cut, very narrow facet can expose a hard tissue boundary.
        # Use its actual face normal only at those nonoutward new corners.
        fallback=np.broadcast_to(np.rint(geo*32767).astype('<i2')[:,None],field.shape)
        field[bad]=fallback[bad]
    surfaces[sid]=faces;words[sid]=field
    key=f'surface_{len(surface_arrays)}';wk=f'words_{len(word_arrays)}'
    arrays[key],arrays[wk]=faces,field;surface_arrays[sid]=key;word_arrays[sid]=wk
    normal_evidence.append(dict(surface=sid,literal_source_corners=literal,
        cut_face_hard_corners=int(bad.sum()),reciprocal_chart_stored_once=True))
part_ids=dict(e['part_ids'])
part_ids.update(dict(zip(meta['root_names'],[
    'ZA:Left inferior pulmonary vein','ZA:Left superior pulmonary vein',
    'ZA:Right inferior pulmonary vein','ZA:Right superior pulmonary vein'])))
packed=pack_finite_cells(p,surfaces,words,meta['cells'],part_ids,o)
atomic_bytes(out/'cells.bin',packed.data);atomic_npz(out/'arrangement.npz',**arrays)
atomic_json(out/'normal-charts.json',normal_evidence)
atomic_json(out/'emission.json',dict(schema='cardiac-shared-bank-emission/v1',
    frame='model:local-ta98:assembled-metres',origin=o.tolist(),
    bank=str(out/'cells.bin'),bank_sha256=digest_file(out/'cells.bin'),
    parts=packed.cells,part_ids=part_ids,cells=meta['cells'],surfaces=surface_arrays,
    corner_words=word_arrays,positions=len(packed.positions),source_positions=len(p),
    surface_triangles=sum(len(f)for f in surfaces.values()),bytes=len(packed.data),
    source_arrangement_sha256=digest_file(source/'arrangement.npz'),
    root_interface_evidence=str(interface_stage/'partition-evidence.json'),
    root_interface_evidence_sha256=digest_file(interface_stage/'partition-evidence.json'),
    source_root_stage=meta['root_source_stage'],
    full_emitted_audits_pending=True,elapsed_seconds=time.monotonic()-start,
    installed=False,anatomical_acceptance=False))
print('Packed',len(part_ids),'cells;',len(packed.data),'bytes',round(time.monotonic()-start,2),'s',flush=True)
