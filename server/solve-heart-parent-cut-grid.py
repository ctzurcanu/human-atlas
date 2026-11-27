"""Solve the coupled shared cut against its unrounded source and parent fields."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.bounded_float32 import untangle
from geometry_tools.distance import triangle_pairs
ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--orientation-only',action='store_true');ap.add_argument('--save-orientation-candidate',action='store_true');args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.source;out=h/args.stage;out.mkdir(exist_ok=False)
m=json.loads((src/'arrangement.json').read_text());z=np.load(src/'arrangement.npz');p=z['positions'];o=np.array(m['origin']);rows=sorted(m['surfaces'],key=lambda r:r['id'])
raw=np.load(h/'pv-root-interfaces-v178/arrangement.npz');word,first=np.unique((raw['positions']-o).astype('<f4'),axis=0,return_index=True);target=raw['positions'][first];anchor=np.load(h/'pv-root-scan-v180/arrangement.npz')['positions'];assert np.array_equal(word.astype(float)+o,anchor)
old=json.loads((h/'regional-bank-v102/emission.json').read_text());oz=np.load(h/'regional-bank-v102/arrangement.npz');op=(oz['positions']-o).astype('<f4').astype(float)+o;rv,rf=z['root_positions'],z['root_faces'];ff=[];fields=[];group=[];offset=0
for row in rows:
    f=z[row['array']];pa=z[row['parents']];q=rv[rf[pa]]if row.get('parentSurface')=='occupied-outer'else op[oz[old['surfaces'][row['id']]][pa]]
    normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])*row.get('parentNormalSign',1);normal/=np.linalg.norm(normal,axis=1)[:,None]
    ff.append(f);fields.append(np.repeat(normal[:,None],3,axis=1))
    if row['id']=='interface:8:-1':group=np.unique(f[pa==4002]).tolist()
f=np.concatenate(ff);fields=np.concatenate(fields);diag=json.loads((src/'diagnostics.json').read_text());assert diag['geometry_sha256']==Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha
pairs=np.array(diag['remaining_global_defect_pairs'],int).reshape(-1,2);guard=LocalCrossingGuard(p,f,pairs,maximum_displacement_m=5e-7,maximum_pair_tests=8000000)
all_original=np.r_[oz['positions'],*[Mesh.load(h/m['root_source_stage']/f'root-{i}.npz').vertices for i in range(4)]];dtype=np.dtype((np.void,12));fixed=np.isin((anchor-o).astype('<f4').view(dtype).ravel(),(all_original-o).astype('<f4').view(dtype).ravel());assert len(group)<=32
support=[];sat=pairs.copy();evidence=[];passed=False
for attempt in range(6):
    if time.monotonic()-start>35:break
    incident=np.flatnonzero(np.isin(f,group).any(1));ids=np.unique(np.r_[incident,sat.ravel()]);used,local=np.unique(f[ids],return_inverse=True);lf=local.reshape(-1,3);lookup=np.full(len(p),-1,int);lookup[used]=np.arange(len(used));lg=lookup[group]
    halfspaces=[]
    for v,g in zip(group,lg):
        bound=5e-8 if fixed[v] else 2.5e-7
        for axis in range(3):
            for sign in [-1,1]:halfspaces.append(dict(control=int(g),normal=(np.eye(3)[axis]*sign).tolist(),point=anchor[v].tolist(),minimum=-bound))
    halfspaces.extend(dict(control=int(lookup[v]),normal=n,point=anchor[v].tolist(),minimum=-bound)for v,n,bound in support)
    try:
        candidate,proof=untangle(target[used],p[used],lf,None if args.orientation_only and attempt==0 else lookup[f[sat]],fields[ids],origin=o,
            mutable_controls=lg,original_controls=lookup[np.array(group)[fixed[group]]],original_steps=16,new_steps=64,
            additional_controls=lg,feasibility_only=False,
            maximum_displacement_m=2.5e-7,time_limit=min(8.,35-(time.monotonic()-start)),maximum_rounds=8,point_halfspaces=halfspaces,separation_axis_limit=1)
        values=candidate[lg];delta=values-anchor[group];length=np.linalg.norm(delta,axis=1);bounds=np.where(fixed[group],5e-8,2.5e-7)
        if np.any(length>bounds):
            support.extend((v,(-d/size).tolist(),float(bound))for v,d,size,bound in zip(group,delta,length,bounds)if size>bound);continue
        proposal=guard.check_move(group,values)
        if proposal is None:
            new=np.array(guard.last_new_pairs,int).reshape(-1,2);sat=np.unique(np.sort(np.r_[sat,new],axis=1),axis=0);evidence.append(dict(attempt=attempt,new_pair_constraints=len(new)))
            if args.save_orientation_candidate:
                trial=p.copy();trial[group]=values;near=guard.pairs_for(group);tri=trial[f[near]]
                bad=np.isin(triangle_pairs(tri[:,0],tri[:,1])['relation'],['crossing','coplanar-overlap'])
                remaining=(set(guard.remaining)-set(map(tuple,near)))|set(map(tuple,near[bad]));p=trial
                guard=LocalCrossingGuard(p,f,sorted(remaining),maximum_displacement_m=5e-7,maximum_pair_tests=8000000)
                evidence.append(dict(unaccepted_candidate_saved=True,complete_affected_pair_cover=True,remaining_crossings=len(remaining)));break
            continue
        guard.commit(proposal);p=guard.positions.copy();evidence.append(proof);passed=not guard.remaining;break
    except (ValueError,TimeoutError,AssertionError)as error:
        evidence.append(dict(attempt=attempt,error=str(error)));break
arrays={k:z[k]for k in z.files};arrays['positions']=p;atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',m)
q=p[f];nn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);bad=np.flatnonzero(np.einsum('ij,ij->i',nn,fields[:,0])<=0)
report=dict(passed=passed and not len(bad),remaining_global_crossings=len(guard.remaining),remaining_parent_reversals=len(bad),maximum_source_adjustment_m=float(np.linalg.norm(p[fixed]-anchor[fixed],axis=1).max()),maximum_total_adjustment_m=float(np.linalg.norm(p-anchor,axis=1).max()),attempts=evidence,elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'diagnostics.json',dict(geometry_sha256=Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha,remaining_global_defect_pairs=sorted(guard.remaining),complete_guarded_worklist=True));atomic_json(out/'construction.json',report);print(json.dumps({k:v for k,v in report.items()if k!='attempts'}),flush=True)
