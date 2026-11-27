"""One literal Float32 coordinate bank for reciprocal atrial/root charts."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True)
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart'
source=h/args.source;out=h/args.stage;out.mkdir(exist_ok=False)
meta=json.loads((source/'arrangement.json').read_text());z=np.load(source/'arrangement.npz');o=np.array(meta['origin'])
p,first,remap=np.unique((z['positions']-o).astype('<f4'),axis=0,return_index=True,return_inverse=True)
p=p.astype(float)+o;arrays={'positions':p};surfaces={};rows=[];encoding=[]
for row in meta['surfaces']:
    original=z[row['array']];f=remap[original];parents=z[row['parents']]
    q=p[f];zero=np.all(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])==0,axis=1)
    f,parents=f[~zero],parents[~zero]
    shift=np.argmin(f,axis=1);cyclic=np.take_along_axis(f,(shift[:,None]+np.arange(3))%3,axis=1)
    sign=np.where(cyclic[:,1]<cyclic[:,2],1,-1)
    canonical,inv=np.unique(np.sort(f,axis=1),axis=0,return_inverse=True)
    coefficient=np.bincount(inv,weights=sign,minlength=len(canonical)).astype('int64')
    assert np.max(abs(coefficient),initial=0)<=1,('Repeated tissue coverage',row['id'])
    selected=[]
    for i in np.flatnonzero(coefficient):
        selected.append(int(np.flatnonzero((inv==i)&(sign==coefficient[i]))[0]))
    f,parents=f[selected],parents[selected]
    encoding.append(dict(surface=row['id'],source_triangles=len(original),
        zero_area_removed=int(zero.sum()),opposed_removed=len(original)-int(zero.sum())-len(f),
        stored_triangles=len(f)))
    if not len(f):continue
    surfaces[row['id']]=f;key=f'surface_{len(rows)}';pk=f'parent_{len(rows)}'
    arrays[key],arrays[pk]=f,parents
    rows.append(dict(row,array=key,parents=pk))
uses={name:[u for u in uu if u['surface']in surfaces]for name,uu in meta['cells'].items()}
checks=[]
for name,uu in uses.items():
    assert uu,name
    ff=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uu])
    mesh=Mesh(p,ff,'model:local-ta98:assembled-metres')
    proof=audit(mesh,require_closed=True,self_intersections=False)
    atomic_json(out/(name.replace(':','-')+'-topology.json'),proof)
    row=dict(name=name,passed=proof['passed'],failures=proof['failures'],
        open_edges=proof['topology']['open_edges'],nonmanifold=proof['topology']['nonmanifold_edges'],
        zero=proof['topology']['zero_area_faces'])
    checks.append(row);print(json.dumps(row),flush=True)
arrays['root_positions']=z['root_positions'];arrays['root_faces']=z['root_faces']
atomic_npz(out/'arrangement.npz',**arrays)
atomic_json(out/'arrangement.json',dict(meta,surfaces=rows,cells=uses,
    source_arrangement_sha256=digest_file(source/'arrangement.npz')))
atomic_json(out/'construction.json',dict(passed=all(c['passed']for c in checks),
    encoding=encoding,cells=checks,elapsed_seconds=time.monotonic()-start,
    literal_float32_aliases_only=True,positive_area_caps_added=0,installed=False))
