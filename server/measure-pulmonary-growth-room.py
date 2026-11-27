"""Conservative lung-growth proposal from actual closest thoracic surfaces.

Per-sample inward clearance halfspaces avoid assigning a false inside/outside
sign to open muscular source patches. A full triangle check must follow.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.open3d_backend import library
from geometry_tools.common import atomic_json,atomic_npz,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-male-ta98:assembled-metres';out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
def read(pid):return reader.read_part(a,parts[pid],frame=frame,catalogue_path=path)
def combine(ids):
 ms=[read(i)for i in ids];off=np.r_[0,np.cumsum([len(m.vertices)for m in ms])];return Mesh(np.concatenate([m.vertices for m in ms]),np.concatenate([m.faces+off[i]for i,m in enumerate(ms)]),frame)
lung=combine(a['pulmonarySegmentIntegration']['partIds']);hidden=set(a.get('defaultHidden',[]));lo,hi=lung.bounds+np.array([[-.03]*3,[.03]*3])
fixed=[p['id']for p in a['parts']if p['system']in ['skeletal','muscular']and p['id']not in hidden and not p.get('suppressed')and not p.get('defaultHidden')and np.all(np.array(p['bounds'])[1]>=lo)and np.all(np.array(p['bounds'])[0]<=hi)]
wall=combine(fixed);o3d=library();origin=wall.vertices.mean(0);scale=np.ptp(wall.vertices,axis=0).max();s=o3d.t.geometry.RaycastingScene(nthreads=2);s.add_triangles(o3d.core.Tensor(((wall.vertices-origin)/scale).astype(np.float32)),o3d.core.Tensor(wall.faces.astype(np.uint32)))
# Every distinct lung vertex and every triangle centre is tested. Opposed
# shared charts repeat constraints without changing their inward direction.
q=np.unique(np.r_[lung.vertices,lung.triangles.mean(1)],axis=0);pivot=(lung.bounds[0]+lung.bounds[1])*.5
c=s.compute_closest_points(o3d.core.Tensor(((q-origin)/scale).astype(np.float32)),nthreads=2)['points'].numpy().astype(float)*scale+origin
d=np.linalg.norm(q-c,axis=1);assert d.min()>1e-9;n=(q-c)/d[:,None]
# Linear variables are extra scale and translation in mm. Clearances already
# below 0.5mm are preserved rather than silently called repaired.
design=np.c_[np.einsum('ij,ij->i',(q-pivot)*1000,n),n]
result=linprog([-1,0,0,0],A_ub=-design,b_ub=d*1000-np.minimum(d*1000,.5),bounds=[(0,.2),(-10,10),(-10,10),(-10,10)],method='highs')
assert result.success,result.message
x=result.x;factor=1+x[0];translation=x[1:]/1000;mat=np.eye(4);mat[:3,:3]*=factor;mat[:3,3]=pivot*(1-factor)+translation
gap=d*1000+design@x;report=dict(source_catalogue_sha256=digest_file(path),lung_matrix=mat.tolist(),relative_scale=float(factor),translation_m=translation.tolist(),sample_count=len(q),minimum_source_clearance_m=float(d.min()),minimum_proposed_plane_clearance_m=float(gap.min()/1000),maximum_linear_clearance_growth=True,source_near_contacts_preserved_not_accepted=True,fixed_ids=fixed,proposal_only=True,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'proposal.json',report);atomic_npz(out/'clearance.npz',points=q,nearest=c,distance_m=d,inward_normals=n,proposed_gap_m=gap/1000);print(json.dumps({k:v for k,v in report.items()if k!='fixed_ids'},indent=2),flush=True)
