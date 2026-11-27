"""Find a smooth caudal sac course in the fixed bones' measured free corridor.

Distance fields propose the course; actual emitted surfaces still require the
independent triangle crossing screen. An open sacral source is not certified
as a closed bone volume by these queries.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.ndimage import gaussian_filter1d
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.open3d_backend import SurfaceIndex,library
from geometry_tools.audit import topology

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';a_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
a=json.loads(a_path.read_text());reader=AtlasReader();indices=[];records=[]
for p in a['parts']:
    if p.get('suppressed')or not p['id'].startswith('LOCAL:female:'):continue
    if 'lumbar_vertebra_'not in p['id']and not p['id'].startswith('LOCAL:female:sacrum:'):continue
    mesh=reader.read_part(a,p,frame='model:local-female-ta98:assembled-metres',catalogue_path=a_path)
    indices.append(SurfaceIndex(mesh));records.append(dict(part_id=p['id'],geometry_sha256=mesh.geometry_sha,closed=topology(mesh)['closed_oriented_manifold']))
def distances(points):
    values=[]
    for index in indices:
        local=((points-index.origin)/index.scale).astype('float32')
        sd=index.scene.compute_signed_distance(library().core.Tensor(local),nthreads=2,nsamples=3).numpy().astype(float)*index.scale
        values.append(sd)
    return np.min(values,axis=0)
ys=np.linspace(.987,1.178,241);zs=np.linspace(-.113,-.017,385)
yy,zz=np.meshgrid(ys,zs,indexing='ij');points=np.c_[np.zeros(yy.size),yy.ravel(),zz.ravel()]
field=distances(points).reshape(yy.shape)
guide=np.interp(ys,[.987,1.005,1.025,1.045,1.065,1.085,1.105,1.125,1.145,1.178],[-.080,-.067,-.054,-.045,-.044,-.039,-.038,-.039,-.042,-.050])
cost=((zs[None]-guide[:,None])/.009)**2-.6*np.minimum(field/.005,1.2)
cost[field<.001]=1e12
# The dynamic program can only progress through a continuous corridor. A
# slope penalty discourages following individual vertebral surface facets.
parent=np.full(cost.shape,-1,int);score=cost[0].copy()
for i in range(1,len(ys)):
    candidates=[]
    for shift in range(-7,8):
        row=np.full(len(zs),1e12)
        target=np.arange(len(zs));source=target-shift;keep=(source>=0)&(source<len(zs))
        row[keep]=score[source[keep]]+.025*shift*shift
        candidates.append(row)
    stack=np.asarray(candidates);choice=stack.argmin(0);score=cost[i]+stack[choice,np.arange(len(zs))]
    parent[i]=np.arange(len(zs))-(choice-7)
j=int(np.argmin(score));assert score[j]<1e10,'No continuous free corridor in the declared sagittal grid'
path=[]
for i in range(len(ys)-1,-1,-1):path.append(zs[j]);j=int(parent[i,j])if i else j
path=np.asarray(path[::-1]);smooth=gaussian_filter1d(path,3.,mode='nearest')
clearance=distances(np.c_[np.zeros(len(ys)),ys,smooth])
assert clearance.min()>.0008,'Smoothed course leaves its measured corridor'
# The circumferential construction starts conservatively inside the measured
# distance to every bone. Later root cuffs may extend to their native exits.
radius=np.minimum(.007,.82*clearance)
atomic_npz(out/'course.npz',y=ys,z=smooth,radius=radius,grid_z=zs,signed_distance=field,unsmoothed_z=path)
atomic_json(out/'course.json',dict(schema='fixed-spinal-canal-course/v1',source_catalogue_sha256=digest_file(a_path),bones=records,
    minimum_center_clearance_m=float(clearance.min()),minimum_proposed_radius_m=float(radius.min()),
    maximum_proposed_radius_m=float(radius.max()),skeleton_changed=False,
    field_scope='Native Open3D distance proposal; open sacral-source sign is not a volume certificate; final exact source triangle crossing checks are mandatory',
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(stations=len(ys),minimum_clearance_mm=float(clearance.min()*1000),radius_mm=[float(radius.min()*1000),float(radius.max()*1000)],elapsed_seconds=time.monotonic()-start)),flush=True)
