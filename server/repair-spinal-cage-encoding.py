"""Repair diagnosed Float32 pairs in independent, 100 nm bounded patches."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.bounded_float32 import untangle

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--prechecks',required=True);ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic()
source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);z=np.load(source/'cells.npz');d=json.loads((source/'cells.json').read_text());origin=np.array([0.,1.35,-.05])
relative=(z['positions']-origin).astype('<f4');encoded,first,mapping=np.unique(relative.astype(float),axis=0,return_index=True,return_inverse=True);encoded+=origin;expected=z['positions'][first]
surfaces={};removed=0
for sid in d['surfaces']:
    f=mapping[z[sid]];q=encoded[f];keep=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0;removed+=int((~keep).sum());surfaces[sid]=f[keep]
checks={}
for path in sorted(Path(args.prechecks).parent.glob(Path(args.prechecks).name+'-precheck-*.json')):
    record=json.loads(path.read_text());assert record['source_sha256']==digest_file(source/'cells.npz');checks.update(record['checks']['cells'])
assert len(checks)==len(d['cells']),'Every emitted tissue needs a complete precheck before rounding repair'
pairs=[];failures=[]
for name,check in checks.items():
    if check['passed']:continue
    assert check['failures']==['self_intersection_check'];summary=check['self_intersections'];assert not summary['truncated']and summary['defect_pairs']==len(summary['examples'])
    f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in d['cells'][name]])
    assert Mesh(encoded,f,'model:local-female-ta98:assembled-metres').geometry_sha==check['geometry_sha256'],'Witnesses must bind the actual emitted faces'
    pairs.extend(f[np.asarray(summary['examples'],int)]);failures.append(name)
assert pairs;pairs=np.asarray(pairs,int);ff=np.concatenate(list(surfaces.values()));q=expected[ff];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];words=np.repeat(normal[:,None,:],3,axis=1)
groups=[]
for index,pair in enumerate(pairs):
    controls=set(map(int,np.unique(ff[np.isin(ff,np.unique(pair)).any(1)])))
    touched=set(map(int,np.flatnonzero(np.isin(ff,list(controls)).any(1))))
    group=dict(pairs={index},controls=controls,faces=touched)
    again=True
    while again:
        again=False
        for previous in list(groups):
            if group['faces']&previous['faces']:
                for key in group:group[key]|=previous[key]
                groups.remove(previous);again=True
    groups.append(group)
positions=encoded.copy();records=[]
atomic_json(out/'patch-plan.json',dict(groups=[dict(pair_ids=sorted(g['pairs']),controls=sorted(g['controls']),incident_faces=sorted(g['faces']))for g in groups],
    source_sha256=digest_file(source/'cells.npz'),maximum_coordinate_move_m=1e-7,per_patch_seconds=10,installed=False))
for index,group in enumerate(groups):
    if time.monotonic()-start>45:raise TimeoutError('Saved patch repair stage reached its native deadline')
    before=positions.copy()
    positions,patch=untangle(expected,positions,ff,pairs[sorted(group['pairs'])],words,origin=origin,mutable_controls=sorted(group['controls']),new_steps=256,
        maximum_displacement_m=1e-7,time_limit=10.,maximum_rounds=8,feasibility_only=True)
    moved=np.flatnonzero(np.any(positions!=before,axis=1));atomic_npz(out/('patch-'+str(index)+'.npz'),controls=moved,positions=positions[moved])
    atomic_json(out/('patch-'+str(index)+'.json'),patch);records.append(patch)
    print(json.dumps(dict(patch=index,groups=len(groups),pairs=len(group['pairs']),moved_controls=len(moved),elapsed_seconds=time.monotonic()-start)),flush=True)
record=dict(patches=records,declared_coordinate_displacement_bound_m=1e-7,maximum_control_displacement_m=float(np.linalg.norm(positions-encoded,axis=1).max()),
    remaining_diagnosed_crossings=0,full_mesh_acceptance=False)
atomic_npz(out/'cells.npz',positions=positions,**surfaces)
record.update(source_cell_bank_sha256=digest_file(source/'cells.npz'),diagnosed_cells=failures,
    diagnosed_pairs=len(pairs),collapsed_encoded_source_faces=removed,complete_assembly_audit_pending=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'rounding-repair.json',record)
meta=dict(d,surfaces={k:len(f)for k,f in surfaces.items()},construction={**d['construction'],'shared_float32_rounding_repair':str(out/'rounding-repair.json'),'rounding_bound_m':1e-7})
atomic_json(out/'cells.json',meta)
print(json.dumps(dict(passed=True,diagnosed_cells=failures,pairs=len(pairs),elapsed_seconds=time.monotonic()-start,result=record)),flush=True)
