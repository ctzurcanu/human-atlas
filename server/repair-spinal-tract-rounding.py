"""Saved, bounded storage repairs for a shared tract partition.

Only derived controls may move, within 100 nm per coordinate. Native tissue
controls stay fixed. Every completed patch is separately reproducible; none
is an installable assembly without the subsequent complete bank audit.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.bounded_float32 import untangle

ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['plan','solve','assemble']);ap.add_argument('source')
ap.add_argument('--precheck');ap.add_argument('--plan');ap.add_argument('--output',required=True)
ap.add_argument('--origin',type=float,nargs=3);ap.add_argument('--group',type=int)
ap.add_argument('--patches',nargs='*');ap.add_argument('--axis-limit',type=int,default=4)
ap.add_argument('--fixed-domain',help='Declared source control bank to keep immutable, rather than the old native core')
args=ap.parse_args();started=time.monotonic()
source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
if args.stage=='plan':
    d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');origin=np.array(args.origin)
    check=json.loads(Path(args.precheck).read_text());assert check['source_sha256']==digest_file(source/'cells.npz')
    q,first,mapping=np.unique((z['positions']-origin).astype('<f4').astype(float),axis=0,return_index=True,return_inverse=True)
    q+=origin;p=z['positions'][first];surfaces={};offset={};offset_end=0
    for sid in d['surfaces']:
        f=mapping[z[sid]];t=q[f];f=f[np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)>0]
        surfaces[sid]=f;offset[sid]=offset_end;offset_end+=len(f)
    pairs=[]
    for name,a in check['checks']['cells'].items():
        if a['passed']:continue
        assert a['failures']==['self_intersection_check'];uses=d['cells'][name]
        f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
        ids=np.concatenate([offset[u['surface']]+np.arange(len(surfaces[u['surface']]))for u in uses])
        mesh=Mesh(q,f,'model:local-female-ta98:assembled-metres');assert mesh.geometry_sha==a['geometry_sha256']
        summary=a['self_intersections']
        if summary['defect_pairs']!=len(summary['examples']):summary=self_contact_summary(mesh,max_pairs=500000,max_tests=2000000,example_limit=2000)
        assert not summary['truncated']and summary['defect_pairs']==len(summary['examples'])
        pairs.extend(ids[np.asarray(summary['examples'],int)])
    pairs=np.unique(np.sort(np.asarray(pairs,int),axis=1),axis=0);f=np.concatenate(list(surfaces.values()));groups=[]
    for number,pair in enumerate(pairs):
        controls=set(map(int,np.unique(f[np.isin(f,np.unique(f[pair])).any(1)])))
        incident=set(map(int,np.flatnonzero(np.isin(f,list(controls)).any(1))))
        group=dict(pairs={number},controls=controls,faces=incident)
        again=True
        while again:
            again=False
            for old in list(groups):
                if old['faces']&group['faces']:
                    for key in group:group[key]|=old[key]
                    groups.remove(old);again=True
        groups.append(group)
    fixed_path=Path(args.fixed_domain).resolve()if args.fixed_domain else ROOT/'.local-models/geometry-workbench/spinal-cord-completion/female-31-segments-v2/cells.npz'
    native=np.load(fixed_path)['positions']
    distance,_=cKDTree(native).query(p);locked=set(map(int,np.flatnonzero(distance<2.1e-12)))
    tq=p[f];normal=np.cross(tq[:,1]-tq[:,0],tq[:,2]-tq[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    atomic_npz(out/'plan.npz',source=p,encoded=q,faces=f,pairs=f[pairs],corner_words=np.repeat(normal[:,None],3,axis=1),**surfaces)
    atomic_json(out/'plan.json',dict(source_sha256=digest_file(source/'cells.npz'),precheck_sha256=digest_file(args.precheck),origin=origin.tolist(),
        groups=[dict(pairs=sorted(g['pairs']),controls=sorted(g['controls']-locked),fixed_controls=len(g['controls']&locked))for g in groups],
        surfaces=list(surfaces),maximum_coordinate_move_m=1e-7,patch_seconds=15.,source_native_controls_fixed=True,
        fixed_control_source=str(fixed_path),fixed_control_source_sha256=digest_file(fixed_path),
        anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started))
    print(json.dumps(dict(groups=len(groups),pairs=len(pairs),elapsed_seconds=time.monotonic()-started)),flush=True)
elif args.stage=='solve':
    plan=Path(args.plan).resolve();d=json.loads((plan/'plan.json').read_text());assert d['source_sha256']==digest_file(source/'cells.npz')
    z=np.load(plan/'plan.npz');g=d['groups'][args.group]
    used=np.unique(z['faces']);lower,upper=z['source'][used,1].min(),z['source'][used,1].max();halfspaces=[]
    for control in g['controls']:
        if min(abs(z['source'][control,1]-lower),abs(z['source'][control,1]-upper))<2.1e-12:
            for sign in [-1.,1.]:halfspaces.append(dict(control=control,normal=[0.,sign,0.],point=z['encoded'][control].tolist(),minimum=0.))
    result,proof=untangle(z['source'],z['encoded'],z['faces'],z['pairs'][g['pairs']],z['corner_words'],origin=d['origin'],
        mutable_controls=g['controls'],new_steps=512,maximum_displacement_m=d['maximum_coordinate_move_m'],
        time_limit=d['patch_seconds'],maximum_rounds=8,feasibility_only=True,separation_axis_limit=args.axis_limit,point_halfspaces=halfspaces)
    moved=np.flatnonzero(np.any(result!=z['encoded'],axis=1));atomic_npz(out/'patch.npz',controls=moved,positions=result[moved])
    atomic_json(out/'patch.json',dict(proof,plan_sha256=digest_file(plan/'plan.npz'),group=args.group,
        full_geometry_audits_pending=True,installed=False,total_elapsed_seconds=time.monotonic()-started))
    print(json.dumps(dict(group=args.group,moved=len(moved),seconds=time.monotonic()-started,solved=True)),flush=True)
else:
    plan=Path(args.plan).resolve();d=json.loads((plan/'plan.json').read_text());z=np.load(plan/'plan.npz');q=z['encoded'].copy();seen=set()
    for path in args.patches:
        patch=Path(path);proof=json.loads((patch/'patch.json').read_text());assert proof['plan_sha256']==digest_file(plan/'plan.npz')
        group=proof['group'];assert group not in seen;seen.add(group);p=np.load(patch/'patch.npz');q[p['controls']]=p['positions']
    assert seen==set(range(len(d['groups']))),'Every independent patch must complete before bank assembly'
    metadata=json.loads((source/'cells.json').read_text());atomic_npz(out/'cells.npz',positions=q,**{sid:z[sid]for sid in d['surfaces']})
    metadata['construction'].update(storage_rounding_plan=str(plan/'plan.json'),maximum_derived_coordinate_move_m=d['maximum_coordinate_move_m'],
        original_native_controls_fixed=True,complete_emitted_bank_checks_pending=True,storage_origin=d['origin'])
    atomic_json(out/'cells.json',metadata)
    atomic_json(out/'rounding.json',dict(completed_groups=sorted(seen),plan_sha256=digest_file(plan/'plan.npz'),
        maximum_coordinate_displacement_m=float(abs(q-z['encoded']).max()),anatomical_acceptance=False,installed=False))
    print(json.dumps(dict(groups=len(seen),assembly_ready_for_audit=True,elapsed_seconds=time.monotonic()-started)),flush=True)
