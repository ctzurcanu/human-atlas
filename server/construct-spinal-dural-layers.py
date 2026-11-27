"""Finite dura and arachnoid with a literal shared interface and open cranial rim."""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import linprog
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.tissue import layer_arrangement
from geometry_tools.algorithms.tube_continuity import finite_wall_cells
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

ap=argparse.ArgumentParser();ap.add_argument('carrier');ap.add_argument('--output',required=True);ap.add_argument('--outer-carrier',action='store_true');args=ap.parse_args()
start=time.monotonic();source=Path(args.carrier).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
data=np.load(source/'surface.npz');v0=data['vertices'];f0=data['faces'];used,inv=np.unique(f0,return_inverse=True);v=v0[used];f=inv.reshape(-1,3)
assert json.loads((source/'surface.json').read_text())['audit']['passed'],'Reject an invalid or incomplete source carrier before constructing coats'
q=v[f];fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);norm=fn/np.linalg.norm(fn,axis=1)[:,None];summed=np.zeros_like(v)
for k in range(3):np.add.at(summed,f[:,k],fn)
directions=summed/np.linalg.norm(summed,axis=1)[:,None]
adj=[[]for _ in v]
for i,tri in enumerate(f):
    for j in tri:adj[int(j)].append(i)
repaired=[]
for j,ids in enumerate(adj):
    nn=norm[ids]
    # Polyhedral miter offsets keep incident planes at nearly equal depth.
    # Unit-normal offsets warp the many tiny cuff triangles differently and
    # can fold an otherwise thin, finite layer.
    miter=np.linalg.lstsq(nn,np.ones(len(nn)),rcond=.01)[0]
    if np.min(nn@miter)>.1 and np.linalg.norm(miter)<4:
        directions[j]=miter;continue
    result=linprog([0,0,0,-1],A_ub=np.c_[-nn,np.ones(len(nn))],b_ub=np.zeros(len(nn)),bounds=[(-1,1)]*3+[(0,2)],method='highs')
    assert result.success and result.x[3]>1e-8,'Source vertex has no common outward offset cone'
    direction=result.x[:3];direction/=np.min(nn@direction);assert np.linalg.norm(direction)<8
    directions[j]=direction
    repaired.append(dict(vertex=j,minimum_face_dot=float(np.min(nn@direction))))
# Regional illustrative dimensions in the range of published human spinal
# tissue measurements. The arachnoid remains 40 um and apposed to the dura;
# no fabricated "subdural layer" occupies the potential interface.
dura=np.interp(v[:,1],[.987,1.178,1.575],[.00012,.00020,.00028])
if args.outer_carrier:dura=np.interp(v[:,1],[.987,1.178,1.575],[.00010,.00012,.00014])
branch=np.clip((np.abs(v[:,0])-.007)/.008,0,1);dura=dura*(1-branch)+.00010*branch
arachnoid=np.full(len(v),.00004)
history=[]
for iteration in range(5):
    # The native source is a zero-thickness dural carrier, not an established
    # outer histological boundary. Retain it at the dura/arachnoid interface:
    # add dura externally and the measured thin arachnoid internally. This
    # avoids assigning the entire two-coat thickness inside narrow cuffs.
    fields=np.stack([v,v-directions*dura[:,None],v-directions*(dura+arachnoid)[:,None]])if args.outer_carrier else np.stack([v+directions*dura[:,None],v,v-directions*arachnoid[:,None]])
    proof=finite_wall_cells(fields,f)
    history.append(dict(iteration=iteration,invalid=proof['invalid_counts_by_interval']))
    if proof['finite_wall_cells_pass']:break
    bad=np.unique(f[np.unique(np.concatenate([np.asarray(x,int)for x in proof['invalid_faces_by_interval']]))])
    proposed=np.maximum(.00008,dura[bad]*.7)
    if np.array_equal(proposed,dura[bad]):break
    dura[bad]=proposed
atomic_npz(out/'fields.npz',fields=fields,faces=f,frame=np.array('model:local-female-ta98:assembled-metres'))
atomic_json(out/'field-check.json',dict(prisms=proof,history=history,normal_cone_repairs=repaired,
    dura_thickness_m=[float(dura.min()),float(dura.max())],arachnoid_thickness_m=.00004,
    references=['https://pubmed.ncbi.nlm.nih.gov/8815466/','https://pubmed.ncbi.nlm.nih.gov/21033455/','https://pubmed.ncbi.nlm.nih.gov/26806431/'],
    source_sha256=digest_file(source/'surface.npz'),anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
assert proof['finite_wall_cells_pass'],proof['invalid_counts_by_interval']
positions,surfaces,cells,_=layer_arrangement(fields,f,['dura','arachnoid'])
checks={}
for name,uu in cells.items():
    ff=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uu])
    checks[name]=topology(Mesh(positions,ff,'model:local-female-ta98:assembled-metres'))
assert all(t['closed_oriented_manifold']for t in checks.values())
atomic_npz(out/'cells.npz',positions=positions,**surfaces)
atomic_json(out/'cells.json',dict(cells=cells,surfaces={k:len(fs)for k,fs in surfaces.items()},attributes={n:dict(native_region=n,axial_region='whole-covering')for n in cells},
    construction=dict(finite_wall_cells_pass=True,source_sha256=digest_file(source/'surface.npz'),cranial_rim_is_annular=True,
        source_carrier_role='Reference-guided outer dura with inward finite tissue coats'if args.outer_carrier else'Reference-guided dura/arachnoid interface; original source did not establish a histological surface',
        dura_arachnoid_interface_shared=True,subarachnoid_cavity_not_plugged=True),
    topology=checks,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(cells=list(cells),triangles_per_interface=len(f),normal_cone_repairs=len(repaired),
    dura_thickness_mm=[float(dura.min()*1000),float(dura.max()*1000)],elapsed_seconds=time.monotonic()-start)),flush=True)
