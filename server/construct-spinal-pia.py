"""A finite source-conforming pial layer, sharing the cord's outer charts.

The selected carrier is witnessed against the native exterior, excluding the
central canal and the artificial cranial truncation. Thickness varies by
region and is illustrative; root sleeves and pial microfenestrations are not
inferred by this construction.
"""
from pathlib import Path
import argparse,copy,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.audit import audit,topology
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim,finite_wall_cells

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');v=z['positions'].copy()
surfaces={k:z[k].copy()for k in d['surfaces']};cells=copy.deepcopy(d['cells']);pia=[];carriers=[];native_uses=[]
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
native=np.load(base/'female-surface-v1/surface.npz');index=SurfaceIndex(Mesh(native['positions'],native['faces'],'model:local-female-ta98:assembled-metres'))
high=float(native['positions'][:,1].max());selected=[]
nq=native['positions'][native['faces']];nfn=np.cross(nq[:,1]-nq[:,0],nq[:,2]-nq[:,0]);nfn/=np.linalg.norm(nfn,axis=1)[:,None]
source_end_cut=((nq[:,:,1].mean(1)>high-.004)&(nfn[:,1]>.2))|((nq[:,:,1].mean(1)<native['positions'][:,1].min()+.0005)&(nfn[:,1]<-.3))
for cell,uu in list(cells.items()):
    for use in list(uu):
        if use.get('neighbor'):continue
        sid=use['surface'];f=surfaces[sid];q=v[f];centres=q.mean(1);query,_=index.closest_points(centres)
        artificial_end_cut=source_end_cut[query['face_ids']]
        take=(query['distance_m']<.0000001)&~artificial_end_cut
        if not take.any():continue
        key='pia-inner:'+sid;surfaces[key]=f[take];carriers.append(f[take]if use['side']==1 else f[take][:,[0,2,1]])
        new_use=dict(use,surface=key,neighbor='pia');uu.append(new_use)
        pia.append(dict(surface=key,side=-use['side'],neighbor=cell));native_uses.append((key,use['side'],cell))
        if (~take).any():surfaces[sid]=f[~take]
        else:uu.remove(use);del surfaces[sid]
carrier=np.concatenate(carriers);outer_ids=np.unique(carrier);remap=np.full(len(v),-1,int);remap[outer_ids]=np.arange(len(outer_ids))
local=v[outer_ids];f=remap[carrier]
# Interpolate the original envelope's field. Averaging newly cut slivers
# independently introduces an artificial kink at each axial partition.
nv,inverse=np.unique(native['positions'],axis=0,return_inverse=True);nf=inverse[native['faces']]
q=nv[nf];fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);fn[source_end_cut]=0;nn=np.zeros_like(nv)
for k in range(3):np.add.at(nn,nf[:,k],fn)
nn/=np.maximum(np.linalg.norm(nn,axis=1)[:,None],1e-100)
native_thickness=np.interp(nv[:,1],[1.147,1.158,1.20,1.44,1.565],[.000085,.000085,.00015,.00015,.00009])
field=nn*native_thickness[:,None]
query,_=index.closest_points(local)
delta=np.einsum('ni,nij->nj',query['barycentric'],field[nf[query['face_ids']]])
thickness=np.linalg.norm(delta,axis=1);outer=local+delta
prisms=finite_wall_cells(np.stack([outer,local]),f)
atomic_json(out/'prism-check.json',prisms)
atomic_npz(out/'pial-fields.npz',inner=local,outer=outer,faces=f)
assert prisms['finite_wall_cells_pass'],prisms['invalid_counts_by_interval']
offset=len(v);v=np.r_[v,outer];mapping=np.full(offset,-1,int);mapping[outer_ids]=offset+np.arange(len(outer_ids))
for key,side,_ in native_uses:
    outer_key=key.replace('pia-inner:','pia-outer:');surfaces[outer_key]=mapping[surfaces[key]]
    pia.append(dict(surface=outer_key,side=side,neighbor=None))
for ordinal,ring in enumerate(boundary_loops(carrier)):
    key='pia-rim:'+str(ordinal)
    # Inner faces are reversed relative to the native carrier; rim ownership
    # follows the outer pial carrier without spanning either opening.
    surfaces[key]=annular_rim(mapping[ring],ring)
    pia.append(dict(surface=key,side=1,neighbor=None))
cells['pia']=pia
faces=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in pia])
proof=audit(Mesh(v,faces,'model:local-female-ta98:assembled-metres'),require_closed=True,max_tests=1000000)
atomic_npz(out/'cells.npz',positions=v,**surfaces)
attributes=d.get('attributes',{});attributes['pia']=dict(native_region='pia',axial_region='whole-covering',thickness_m=[float(thickness.min()),float(thickness.max())])
atomic_json(out/'cells.json',dict(cells=cells,surfaces={k:len(f)for k,f in surfaces.items()},attributes=attributes,bands=d.get('bands'),
    construction=dict(source_sha256=digest_file(source/'cells.npz'),source_envelope_preserved=True,
        pial_inner_interface_shared=True,cranial_continuation_open=True,central_canal_not_coated=True,prism_check=prisms),
    pia_audit=proof,remaining_limits=['Root sleeves, filum and denticulate attachments are separate pending interfaces.',
        'Thickness is a reference-guided regional illustration; donor-specific thickness and microscopic fenestrations are not established.'],
    references=['https://pubmed.ncbi.nlm.nih.gov/9558933/','https://pubmed.ncbi.nlm.nih.gov/8886963/'],
    passed=proof['passed'],anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
assert proof['passed'],proof['failures']
print(json.dumps(dict(parts=len(cells),pia_faces=len(faces),thickness_mm=[float(thickness.min()*1000),float(thickness.max()*1000)],
    elapsed_seconds=time.monotonic()-start)),flush=True)
