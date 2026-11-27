"""Condition shared tissue incidences before coincident Float32 controls weld.

Only existing endpoints may survive a short-edge contraction. Every affected
material must pass the manifold link condition; physical/macro boundary controls
are locked. The independent packed-cell audits remain mandatory afterwards.
"""
from pathlib import Path
import argparse, json, time
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.algorithms.shared_edge_collapse import collapse_short_edges
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

ap=argparse.ArgumentParser()
ap.add_argument('source'); ap.add_argument('--domain',required=True)
ap.add_argument('--output',required=True)
ap.add_argument('--origin',type=float,nargs=3,default=[0.,1.35,-.05])
ap.add_argument('--limit',type=float,default=1e-7)
ap.add_argument('--seconds',type=float,default=5.)
args=ap.parse_args(); started=time.monotonic()
source=Path(args.source).resolve(); domain=Path(args.domain).resolve()
out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False)
m=json.loads((source/'cells.json').read_text()); z=np.load(source/'cells.npz')
assert m['passed'], 'A failed source is not a storage candidate'
p,inv=np.unique(z['positions'],axis=0,return_inverse=True)
origin=np.array(args.origin); q=(p-origin).astype('<f4').astype(float)+origin
surfaces={sid:inv[z[sid]] for sid in m['surfaces']}
words={}
for sid,f in surfaces.items():
    t=p[f]; n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    n/=np.linalg.norm(n,axis=1)[:,None]
    words[sid]=np.repeat(n[:,None],3,axis=1)
d=np.load(domain/'domain.npz')
protected=d['positions'][d['interior_material']<0]
locked=np.flatnonzero(cKDTree(protected).query(p)[0]<2.1e-12)
result=collapse_short_edges(q,surfaces,m['cells'],words,
    desired_positions=p,maximum_displacement_m=args.limit,
    maximum_steps=64,wall_seconds=args.seconds,candidate_face_rings=1,
    preserve_material_membership=False,locked=locked)
atomic_json(out/'contractions.json',result['evidence'])
assert result['evidence']['passed'], result['evidence']['budget_or_stop_reason']
surfaces={sid:f for sid,f in result['surfaces'].items() if len(f)}
cells={name:[u for u in uses if u['surface']in surfaces]
       for name,uses in m['cells'].items()}
assert all(cells.values()), 'No anatomical tissue may disappear at storage'
checks={}
for name,uses in cells.items():
    f=np.concatenate([surfaces[u['surface']] if u['side']==1
                      else surfaces[u['surface']][:,[0,2,1]] for u in uses])
    used,inv=np.unique(f,return_inverse=True)
    check=topology(Mesh(q[used],inv.reshape(-1,3),'spinal-storage-candidate'))
    checks[name]={k:check[k]for k in ['closed_oriented_manifold','open_edges',
        'nonmanifold_edges','zero_area_faces','components','oriented_triangle_integral_m3']}
assert all(c['closed_oriented_manifold']for c in checks.values()), checks
atomic_npz(out/'cells.npz',positions=q,**surfaces)
atomic_json(out/'cells.json',dict(m,cells=cells,
    surfaces={sid:len(f)for sid,f in surfaces.items()},topology=checks,
    construction={**m.get('construction',{}),
        'pre_storage_source_sha256':digest_file(source/'cells.npz'),
        'locked_domain_sha256':digest_file(domain/'domain.npz'),
        'link_checked_storage_contractions':result['evidence'],
        'full_emitted_checks_pending':True},
    elapsed_seconds=time.monotonic()-started,installed=False))
print(json.dumps(dict(passed=True,cells=len(cells),
    contractions=len(result['evidence']['contractions']),
    maximum_source_displacement_m=result['evidence']['maximum_accepted_cluster_displacement_m'],
    elapsed_seconds=time.monotonic()-started)),flush=True)
