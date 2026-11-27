"""Bounded actual Float32 serialization checks, without altering source controls."""
from pathlib import Path
import json,sys,time,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json
from geometry_tools.model import Mesh
from geometry_tools.contact import self_contact_summary
args=argparse.ArgumentParser();args.add_argument('--version',default='v5');args=args.parse_args()
T=time.monotonic();D=ROOT/('.local-models/geometry-workbench/circulatory-completion/male-left-lingual-finite-network-'+args.version);layout=json.loads((D/'owned-bank-v1/tissue.json').read_text());z=np.load(layout['arrangement']);v=z['positions'];faces=[]
for ref in layout['cells']['lingual']:
    f=z[layout['surfaces'][ref['surface']]];faces.append(f if ref['side']==1 else f[:,::-1])
f=np.concatenate(faces);u=np.unique(f);center=(v[u].min(0)+v[u].max(0))/2;origins=[np.zeros(3),center,np.array([.025,1.625,.035]),np.array([.020,1.624,.040]),np.array([.016594554721975487,1.623810152874436,.04152941869862918])];rows=[]
for origin in origins:
    p=(v-origin).astype('<f4').astype(float)+origin;m=Mesh(p,f,layout['frame']);summary,_=self_contact_summary(m,max_pairs=250000,max_tests=2000000);rows.append(dict(origin_m=origin.tolist(),maximum_position_roundoff_m=float(np.linalg.norm(p[u]-v[u],axis=1).max()),self_contact=summary))
passed=[r for r in rows if r['self_contact']['complete_within_pair_budget']and r['self_contact']['defect_pairs']==0];best=min(passed,key=lambda r:r['maximum_position_roundoff_m'])if passed else None
atomic_json(D/'position-encoding-v1.json',dict(schema='source-control-float32-serialization/v1',shape_control_change=False,tested=rows,selected=best,elapsed_seconds=time.monotonic()-T));print(json.dumps(dict(selected=best,tested=[dict(origin=r['origin_m'],error=r['maximum_position_roundoff_m'],defects=r['self_contact']['defect_pairs'])for r in rows],elapsed_seconds=time.monotonic()-T)))
