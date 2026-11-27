"""Literal, uncapped plane contours along named pulmonary-vein entrances."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.algorithms.simple_compartment_section import polygon_contains

def contours(v,f,c,n):
    s=(v-c)@n;ss=s[f];take=np.flatnonzero((ss.min(1)<0)&(ss.max(1)>0));tri=f[take]
    edge=np.stack([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]],axis=1);hit=s[edge[:,:,0]]*s[edge[:,:,1]]<0
    if np.any(hit.sum(1)!=2):raise ValueError('A section passes through a source control; select a separate physical station')
    edge=np.sort(edge[hit],axis=1);unique,iv=np.unique(edge,axis=0,return_inverse=True);segments=iv.reshape(-1,2);fraction=s[unique[:,0]]/(s[unique[:,0]]-s[unique[:,1]]);p=v[unique[:,0]]+fraction[:,None]*(v[unique[:,1]]-v[unique[:,0]])
    adjacency=[[]for _ in p]
    for a,b in segments:adjacency[a].append(int(b));adjacency[b].append(int(a))
    if any(len(x)!=2 for x in adjacency):raise ValueError('Actual section is not a collection of closed curves')
    remaining=set(range(len(p)));loops=[]
    while remaining:
        first=min(remaining);q=first;previous=-1;ids=[]
        while q in remaining:
            ids.append(q);remaining.remove(q);other=adjacency[q];previous,q=q,other[0]if other[0]!=previous else other[1]
        if q!=first:raise ValueError('A measured contour does not close')
        loops.append(p[ids])
    return loops

ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--all-roots',action='store_true');args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);m=Mesh.load(b/'blood-spaces-v126/installed-boundary.npz');gates=json.loads((b/'blood-spaces-v126/blood-spaces.json').read_text())['gates'];seed=np.array([.018,1.405,-.012]);allrows=[]
for g in gates:
 if not args.all_roots and'pulmonary vein'not in g['part_id']:continue
 c0=np.array(g['center_m']);direction=seed-c0
 if'pulmonary vein'not in g['part_id']:
  direction=np.array(g['normal']);target=np.array([.025,1.36,.06]if'Pulmonary trunk'in g['part_id']else[.05,1.36,.025]if'Ascending aorta'in g['part_id']else[-.022,1.39,.022]);direction*=1 if np.dot(direction,target-c0)>0 else-1
 direction/=np.linalg.norm(direction);axis=np.eye(3)[np.argmin(abs(direction))];u=np.cross(direction,axis);u/=np.linalg.norm(u);v=np.cross(direction,u);basis=np.c_[u,v];rows=[]
 for t in [-.008,-.004,-.002,0,.001,.002,.003,.004,.005,.006,.008,.010,.012,.014,.016,.018,.020]:
  c=c0+t*direction;loops=contours(m.vertices,m.faces,c,direction);enclosing=[]
  for p in loops:
   xy=(p-c)@basis
   if not polygon_contains(xy,np.zeros(2)):continue
   area=abs(np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1)))/2;enclosing.append(dict(area_mm2=float(area*1e6),maximum_radius_mm=float(np.linalg.norm(xy,axis=1).max()*1000),points_m=p.tolist()))
  enclosing.sort(key=lambda q:q['area_mm2']);row=dict(distance_from_source_mouth_mm=t*1000,point_m=c.tolist(),normal=direction.tolist(),enclosing=enclosing,other_curve_count=len(loops)-len(enclosing));rows.append(row)
 root=dict(part_id=g['part_id'],stations=rows);allrows.append(root);atomic_json(out/f'root-{len(allrows)-1}.json',root);print(g['part_id'],[(r['distance_from_source_mouth_mm'],[(round(x['area_mm2'],1),round(x['maximum_radius_mm'],1))for x in r['enclosing']])for r in rows],flush=True)
atomic_json(out/'necks.json',dict(measured=True,whole_wall_sha256=m.geometry_sha,roots=allrows,elapsed_seconds=time.monotonic()-start,no_caps=True,geometry_changed=False,anatomical_acceptance=False))
