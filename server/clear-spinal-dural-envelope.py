"""Smooth radial clearance of the dural carrier around the literal cord.

The map is a positive radial expansion with fixed caudal sac and distant
root exits. It moves no skeleton, cord, pia or peripheral nerve control.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic()
source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);base=source.parent
z=np.load(source/'surface.npz');v=z['vertices'];f=z['faces'];cord=np.load(base/'female-surface-v1/surface.npz');tri=cord['positions'][cord['faces']]
ys=np.linspace(tri[:,:,1].min()+1e-7,tri[:,:,1].max()-1e-7,850);centers=[]
for y in ys:
    pieces=[]
    for i,j in [(0,1),(1,2),(2,0)]:
        d=(tri[:,i,1]-y)*(tri[:,j,1]-y);sel=d<0;a,b=tri[sel,i],tri[sel,j]
        if len(a):pieces.append(a+(b-a)*(y-a[:,1])[:,None]/(b[:,1]-a[:,1])[:,None])
    p=np.concatenate(pieces);centers.append((p[:,[0,2]].min(0)+p[:,[0,2]].max(0))/2)
centers=gaussian_filter1d(np.asarray(centers),1.,axis=0);center=PchipInterpolator(ys,centers)(np.clip(v[:,1],ys[0],ys[-1]))
delta=v[:,[0,2]]-center;r=np.linalg.norm(delta,axis=1)
def smooth(x):
    x=np.clip(x,0,1);return x*x*(3-2*x)
radial_gate=1-smooth((r-.008)/.012)
axial_gate=smooth((v[:,1]-1.171)/.008)*(1-smooth((v[:,1]-1.561)/.014))
amount=.00085*radial_gate*axial_gate
factor=1+amount/np.sqrt(r*r+.002*.002);vv=v.copy();vv[:,[0,2]]=center+delta*factor[:,None]
mesh=Mesh(vv,f,'model:local-female-ta98:assembled-metres');proof=audit(mesh,require_closed=False,max_pairs=250000,max_tests=2000000)
atomic_npz(out/'surface.npz',vertices=vv,faces=f,frame=np.array(mesh.frame))
record=json.loads((source/'surface.json').read_text());record.update(source_sha256=digest_file(source/'surface.npz'),
    original_upper_positions_literal=False,clearance_map=dict(kind='smooth radially monotone expansion about literal cord sections',
        maximum_parameter_m=.00085,caudal_stationary_below_y_m=1.171,distant_root_exits_stationary_above_radius_m=.020,
        cord_and_pia_literal=True,skeleton_literal=True,maximum_carrier_displacement_m=float(np.linalg.norm(vv-v,axis=1).max())),
    audit=proof,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'surface.json',record);assert proof['passed'],proof['failures']
print(json.dumps(dict(passed=True,self_crossings=proof['self_intersections']['defect_pairs'],
    maximum_carrier_displacement_mm=record['clearance_map']['maximum_carrier_displacement_m']*1000,elapsed_seconds=time.monotonic()-start)),flush=True)
