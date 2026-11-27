"""Re-triangulate numerically planar slivers before coordinate storage.

One physical chart is changed once for both owners. No control moves and no
boundary edge changes. Plane residual is limited to two picometres; this is
numerical conditioning, not a free anatomical remodelling tolerance.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.planar_sliver_flip import flip_planar_slivers
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);args=ap.parse_args()
started=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');v=z['positions'];surfaces={};proof={}
for sid in d['surfaces']:
    if time.monotonic()-started>25:raise TimeoutError('Planar facet stage reached its saved 25-second budget')
    f,r=flip_planar_slivers(v,z[sid],minimum_altitude_m=5e-8,plane_tolerance_m=2e-12,maximum_flips=100,seconds=2.)
    surfaces[sid]=f
    if r['flips']:proof[sid]=r
for name,uses in d['cells'].items():
    f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
    t=topology(Mesh(v,f,'model:local-female-ta98:assembled-metres'));assert t['closed_oriented_manifold'],name
atomic_npz(out/'cells.npz',positions=v,**surfaces)
d['construction'].update(planar_diagonal_conditioning=str(out/'facet-conditioning.json'),source_controls_literal=True)
atomic_json(out/'cells.json',d)
atomic_json(out/'facet-conditioning.json',dict(charts=proof,source_sha256=digest_file(source/'cells.npz'),
    controls_literal=True,shared_boundary_edges_literal=True,full_emitted_audit_pending=True,anatomical_acceptance=False,
    installed=False,elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(charts=len(proof),flips=sum(len(r['flips'])for r in proof.values()),
    controls_moved=0,elapsed_seconds=time.monotonic()-started)),flush=True)
