"""Fit the pial envelope inside the fixed, bone-clear meningeal carrier.

Solve the exact interval of permissible cord eccentricity for every candidate
area-preserving section shape. Do not smooth a result after its constraints.
This produces a candidate field; actual mapped triangles must still be audited.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.interpolate import PchipInterpolator
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file

ap=argparse.ArgumentParser();ap.add_argument('course');ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic();source=Path(args.course).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
c=dict(np.load(source/'course.npz'));base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';pia=np.load(base/'female-pia-v4/pial-fields.npz');tri=pia['outer'][pia['faces']]
ys=c['y'];shapes=np.linspace(.85,1.15,61);beta=np.zeros((len(ys),len(shapes)));valid=np.ones_like(beta,dtype=bool);lower=np.full_like(beta,-.0015);upper=np.full_like(beta,.0015)
for i,y in enumerate(ys):
    if y<tri[:,:,1].min()or y>tri[:,:,1].max():continue
    cuts=[]
    for a,b in [(0,1),(1,2),(2,0)]:
        selected=(tri[:,a,1]-y)*(tri[:,b,1]-y)<0;pa,pb=tri[selected,a],tri[selected,b]
        if len(pa):cuts.append(pa+(pb-pa)*(y-pa[:,1])[:,None]/(pb[:,1]-pa[:,1])[:,None])
    if not cuts:continue
    p=np.concatenate(cuts);dx=p[:,None,0]-c['center_x'][i];dz=p[:,None,2]-c['native_center_z'][i]
    rx,rz=c['rx'][i]-.00021,c['rz'][i]-.00021
    limit=.985-(dx/shapes[None,:]/rx)**2;valid[i]=np.all(limit>0,axis=0)
    reach=rz*np.sqrt(np.maximum(limit,0));lo=np.max(dz*shapes[None,:]-reach,axis=0);hi=np.min(dz*shapes[None,:]+reach,axis=0)
    lower[i]=np.maximum(lower[i],lo);upper[i]=np.minimum(upper[i],hi);valid[i]&=lower[i]<=upper[i]
    beta[i]=np.clip(c['cord_eccentricity_z'][i],lower[i],upper[i])
bad=np.flatnonzero(~valid.any(1));atomic_npz(out/'feasible-intervals.npz',y=ys,shapes=shapes,lower=lower,upper=upper,valid=valid)
assert not len(bad),dict(no_feasible_section_y=ys[bad].tolist())
node=((shapes[None,:]-c['section_ap_scale'][:,None])/.02)**2+((beta-c['cord_eccentricity_z'][:,None])/.0003)**2*.1
cost=np.where(valid[0],node[0],1e20);parents=np.zeros_like(valid,dtype=int)
for i in range(1,len(ys)):
    values=cost[:,None]+((shapes[:,None]-shapes[None,:])/.015)**2+((beta[i-1,:,None]-beta[i,None,:])/.0002)**2
    parents[i]=values.argmin(0);cost=np.where(valid[i],values.min(0)+node[i],1e20)
selected=np.zeros(len(ys),int);selected[-1]=int(cost.argmin())
for i in range(len(ys)-1,0,-1):selected[i-1]=parents[i,selected[i]]
aa=shapes[selected];bb=beta[np.arange(len(ys)),selected]
c.update(section_ap_scale=aa,cord_eccentricity_z=bb,cord_shift_z=c['center_z']-c['native_center_z']-bb)
atomic_npz(out/'course.npz',**c)
record=json.loads((source/'course.json').read_text());record.update(source_course_sha256=digest_file(source/'course.npz'),
    meningeal_course_and_envelope_literal=True,cord_inside_meninges_constraints=True,post_fit_smoothing=False,
    ap_scale=[float(aa.min()),float(aa.max())],eccentricity_m=[float(bb.min()),float(bb.max())],
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'course.json',record)
print(json.dumps(dict(passed=True,ap_scale=record['ap_scale'],maximum_eccentricity_mm=float(abs(bb).max()*1000),elapsed_seconds=time.monotonic()-start)),flush=True)
