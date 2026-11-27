"""Preserve source control identities while selecting bounded Float32 words.

Distinct source controls must not be welded merely because nearest rounding
lands them on the same GPU coordinate. Anatomical boundary controls are fixed;
diagnosed interior/derived controls can move, within 100 nm per coordinate.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.bounded_float32 import untangle

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--domain',required=True);ap.add_argument('--output',required=True)
ap.add_argument('--origin',type=float,nargs=3,default=[0.,1.35,-.05]);args=ap.parse_args();started=time.monotonic()
source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');p,mapping=np.unique(z['positions'],axis=0,return_inverse=True)
surfaces={sid:mapping[z[sid]]for sid in d['surfaces']};f=np.concatenate(list(surfaces.values()));origin=np.array(args.origin)
q=(p-origin).astype('<f4').astype(float)+origin;tri=p[f];n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
stored=q[f];normal=np.cross(stored[:,1]-stored[:,0],stored[:,2]-stored[:,0]);bad=np.einsum('ij,ij->i',n,normal)<=0
active=np.unique(f);_,inverse,count=np.unique(q[active],axis=0,return_inverse=True,return_counts=True)
duplicates=[active[inverse==i]for i in np.flatnonzero(count>1)]
controls=np.unique(np.r_[f[bad].ravel(),*[g for g in duplicates]])
if len(controls):controls=np.unique(f[np.isin(f,controls).any(1)])
domain=np.load(args.domain);native=domain['positions'];boundary=native[domain['interior_material']<0]
distance,_=cKDTree(boundary).query(p);mutable=np.flatnonzero(distance>2.1e-12)
if len(controls):
    halfspaces=[];low,high=p[active,1].min(),p[active,1].max()
    for control in controls:
        if min(abs(p[control,1]-low),abs(p[control,1]-high))<2.1e-12:
            for sign in [-1.,1.]:halfspaces.append(dict(control=int(control),normal=[0.,sign,0.],point=q[control].tolist(),minimum=0.))
    result,proof=untangle(p,q,f,None,np.repeat(n[:,None],3,axis=1),origin=origin,mutable_controls=mutable,
        additional_controls=controls,new_steps=512,maximum_displacement_m=1e-7,time_limit=10.,
        maximum_rounds=8,feasibility_only=False,point_halfspaces=halfspaces)
else:result=q;proof=dict(changed_controls=[],maximum_coordinate_displacement_m=0.)
atomic_npz(out/'cells.npz',positions=result,**surfaces)
atomic_json(out/'cells.json',dict(d,construction=dict(d['construction'],source_controls_not_merged_at_encoding=True,
    structured_anatomical_boundary_controls_fixed=True,storage_repair=str(out/'encoding.json'),storage_origin=origin.tolist())))
atomic_json(out/'encoding.json',dict(proof,source_sha256=digest_file(source/'cells.npz'),domain_sha256=digest_file(args.domain),
    diagnosed_orientation_facets=int(bad.sum()),diagnosed_coincident_control_groups=len(duplicates),
    full_emitted_bank_audit_pending=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(diagnosed_facets=int(bad.sum()),coincident_groups=len(duplicates),moved_controls=len(proof['changed_controls']),
    maximum_coordinate_move_m=proof['maximum_coordinate_displacement_m'],seconds=time.monotonic()-started)),flush=True)
