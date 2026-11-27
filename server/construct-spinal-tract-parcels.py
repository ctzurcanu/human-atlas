"""Disjoint probabilistic white-matter parcels in one native segment.

The shared calculation grid is warped using actual fitted cord sections.
It never replaces the cord wall. The resulting closed masks must next be
intersected with native white matter and audited; they are not anatomy yet.
"""
from pathlib import Path
from itertools import permutations
import argparse,json,time
import numpy as np
from scipy.interpolate import interp1d
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.labelled_tetrahedra import partition_labelled_tetrahedra
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

ap=argparse.ArgumentParser();ap.add_argument('--fields',required=True);ap.add_argument('--segment',required=True)
ap.add_argument('--output',required=True);ap.add_argument('--seconds',type=float,default=25.)
args=ap.parse_args();started=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';field_path=Path(args.fields).resolve()
metadata=json.loads((field_path/'fields.json').read_text());plan=json.loads((field_path/'plan.json').read_text())['parts'];field=np.load(field_path/'fields.npz')
registration=json.loads((base/'female-root-regions-v2/root-regions.json').read_text())['segments']
names=[r['name']for r in registration];index=names.index(args.segment);low,high=registration[index]['y_range_m']
assert index<26,'No distal template samples: do not extrapolate S2-Co1'
source=base/'female-cage-rounded-v3';d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');positions=z['positions']
# Pia/cord shared charts bound the actual fitted core, not the larger dural sac.
boundary=set()
for name,uses in d['cells'].items():
    if name.split('@')[0]!='pia':continue
    for use in uses:
        if use.get('neighbor')and not use['neighbor'].startswith('pia'):boundary.add(use['surface'])
tri=positions[np.concatenate([z[s]for s in sorted(boundary)])]
ylow,yhigh=tri[:,:,1].min(),tri[:,:,1].max()
def section(y):
    yy=np.clip(y,ylow+1e-5,yhigh-1e-5);distance=tri[:,:,1]-yy;hits=[tri[distance==0]]
    for a,b in [(0,1),(1,2),(2,0)]:
        keep=distance[:,a]*distance[:,b]<0
        weight=distance[keep,a]/(distance[keep,a]-distance[keep,b])
        hits.append(tri[keep,a]+weight[:,None]*(tri[keep,b]-tri[keep,a]))
    q=np.concatenate(hits);lo,hi=q[:,[0,2]].min(0),q[:,[0,2]].max(0)
    return (lo+hi)/2,(hi-lo)/2

# A slight longitudinal extension makes these classification solids cross,
# rather than terminate numerically just short of, the native segment caps.
cut_low=float(np.float32(low-1.35))+1.35;cut_high=float(np.float32(high-1.35))+1.35
# Include the literal native cut planes as calculation stations. Adjacent
# segments then evaluate the identical field and transverse grid at their
# shared boundary, rather than two slightly different interpolated traces.
ys=np.r_[cut_low-2e-5,np.linspace(cut_low,cut_high,max(3,int(np.ceil((high-low)/.006))+1)),cut_high+2e-5]
uv=field['normalized_controls'];axis=field['axis'];n=len(axis);points=[];scores=[]
middle=field['segment_y'][:30].mean(1);order=np.argsort(middle)
interp=interp1d(middle[order],field['scores'][order],axis=0,bounds_error=False,
    fill_value=(field['scores'][order[0]],field['scores'][order[-1]]))
eligible=np.array([names.index(row['first_level'])<=index<=names.index(row['last_level'])for row in plan])
for y in ys:
    center,radii=section(y);xy=center+uv*radii;points.append(np.c_[xy[:,0],np.full(len(xy),y),xy[:,1]])
    prob=interp(y);prob[:,~eligible]=-1.
    scores.append(np.c_[np.full(len(prob),metadata['probability_floor']),prob])
p=np.concatenate(points);s=np.concatenate(scores)
shape=[len(ys),n,n];base_nodes=np.arange(len(p)).reshape(shape)[:-1,:-1,:-1].ravel()
strides=np.array([n*n,n,1]);tet=[]
for order in permutations(range(3)):
    tet.append(base_nodes[:,None]+np.r_[0,np.cumsum(strides[list(order)])])
t=np.concatenate(tet);det=np.linalg.det(p[t[:,1:]]-p[t[:,:1]]);assert np.all(det!=0)
bad=det<0;t[bad]=t[bad][:,[0,2,1,3]]
active=[0,*[i+1 for i,e in enumerate(eligible)if e]]
result=partition_labelled_tetrahedra(p,t,s[:,active],seconds=args.seconds,maximum_tetrahedra=25000,
    identical_score_priority=list(range(len(active))))
arrays={'positions':result['positions']};cells={};checks={}
for local,f in result['cells'].items():
    global_label=active[local]
    if global_label==0 or not len(f):continue
    name=plan[global_label-1]['name'];sid='parcel:'+name;arrays[sid]=f
    cells[name]=[dict(surface=sid,side=1,neighbor=None)]
    q=topology(Mesh(result['positions'],f,'native-spinal-tract-classifier'))
    checks[name]={k:q[k]for k in ['closed_oriented_manifold','open_edges','nonmanifold_edges','invalid_vertex_links','zero_area_faces','oriented_triangle_integral_m3']}
    assert q['closed_oriented_manifold'],(name,checks[name])
atomic_npz(out/'masks.npz',**arrays)
atomic_json(out/'masks.json',dict(cells=cells,segment=args.segment,y_range_m=[low,high],fields_sha256=digest_file(field_path/'fields.npz'),
    source_native_tissue_sha256=digest_file(source/'cells.npz'),topology=checks,construction=result['evidence'],
    literal_native_segment_cap_stations=[cut_low,cut_high],
    anatomical_acceptance=False,construction_only=True,not_native_intersection_yet=True,installed=False,
    elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(segment=args.segment,finite_parcels=len(cells),triangles=sum(len(a)for k,a in arrays.items()if k!='positions'),
    elapsed_seconds=time.monotonic()-started,construction_only=True)),flush=True)
