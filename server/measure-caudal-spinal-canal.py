"""Fixed-bone ray witnesses for the caudal dural sac; no bone deformation."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.open3d_backend import SurfaceIndex

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
catalogue=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
a=json.loads(catalogue.read_text());r=AtlasReader();positions=[];faces=[];sources=[]
for p in a['parts']:
    if p.get('suppressed')or not p['id'].startswith('LOCAL:female:'):continue
    if 'lumbar_vertebra_'not in p['id']and not p['id'].startswith('LOCAL:female:sacrum:'):continue
    m=r.read_part(a,p,frame='model:local-female-ta98:assembled-metres',catalogue_path=catalogue)
    faces.append(m.faces+sum(len(v)for v in positions));positions.append(m.vertices)
    sources.append(dict(id=p['id'],geometry_sha256=m.geometry_sha,bounds_m=m.bounds.tolist()))
v=np.concatenate(positions);f=np.concatenate(faces);q=v[f]
nonzero=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0
mesh=Mesh(v,f[nonzero],'model:local-female-ta98:assembled-metres');index=SurfaceIndex(mesh)
ys=np.linspace(.94,1.18,121);origins=np.c_[np.zeros(len(ys)),ys,np.full(len(ys),.15)]
direction=np.tile([0.,0.,-1.],(len(ys),1));hits=[]
for step in range(8):
    ray,proof=index.cast_rays(origins,direction);hits.append(ray)
    origins[ray['hit']]=ray['points'][ray['hit']]+direction[ray['hit']]*.000005
    origins[~ray['hit'],2]=-.2
points=np.stack([h['points']for h in hits],axis=1)
atomic_npz(out/'canal-rays.npz',y=ys,hit_points=points,
    face_ids=np.stack([h['face_ids']for h in hits],axis=1))
rows=[]
for y,pp in zip(ys,points):
    zz=pp[np.isfinite(pp).all(1),2];zz=zz[np.r_[True,np.abs(np.diff(zz))>.00002]]if len(zz)else zz
    # Four ordered surface hits witness body/front-lamina separation. They
    # are retained as witnesses, not an occupancy claim for an open source.
    if len(zz)>=4 and zz[1]-zz[2]>.003:
        rows.append(dict(y_m=float(y),anterior_bone_z_m=float(zz[1]),posterior_bone_z_m=float(zz[2]),
            canal_midpoint_z_m=float((zz[1]+zz[2])/2),ap_clearance_m=float(zz[1]-zz[2])))
atomic_json(out/'canal-witnesses.json',dict(source_catalogue_sha256=digest_file(catalogue),sources=sources,
    geometry_sha256=mesh.geometry_sha,rows=rows,excluded_zero_area_source_faces=int((~nonzero).sum()),
    approximate_rays=True,anatomical_acceptance=False,skeleton_changed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(rows=len(rows),examples=rows[::max(1,len(rows)//12)],sources=sources,elapsed_seconds=time.monotonic()-start)),flush=True)
