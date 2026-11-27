"""Measure actual source cross-sections from papillary ends toward leaflet sails."""
from pathlib import Path
import json,sys,time
import numpy as np
import igl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';src=h/'apparatus-source-v203';out=h/'chordal-sections-v207';out.mkdir(exist_ok=False)
def contours(p,f,values,level,origin,n):
    d=values-level;edges=np.stack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]],axis=1);cross=d[edges[:,:,0]]*d[edges[:,:,1]]<0;active=np.flatnonzero(cross.sum(1)==2)
    cut=edges[active][cross[active]];keys,inv=np.unique(np.sort(cut,axis=1),axis=0,return_inverse=True);a,b=keys.T;t=d[a]/(d[a]-d[b]);q=p[a]+t[:,None]*(p[b]-p[a]);segments=inv.reshape(-1,2);adj={}
    for a,b in segments:adj.setdefault(int(a),[]).append(int(b));adj.setdefault(int(b),[]).append(int(a))
    assert all(len(v)==2 for v in adj.values());todo=set(adj);loops=[]
    while todo:
        first=min(todo);loop=[];previous=None;at=first
        while at in todo:
            todo.remove(at);loop.append(at);choices=adj[at];nxt=choices[0]if choices[0]!=previous else choices[1];previous,at=at,nxt
        assert at==first;loops.append(q[loop])
    u=np.cross(n,np.eye(3)[np.argmin(abs(n))]);u/=np.linalg.norm(u);v=np.cross(n,u);rows=[]
    for poly in loops:
        xy=(poly-origin)@np.c_[u,v];area=abs(np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1)))/2;perimeter=np.linalg.norm(poly-np.roll(poly,1,axis=0),axis=1).sum();cov=np.linalg.eigvalsh(np.cov(xy.T));rows.append(dict(points_m=poly.tolist(),area_m2=float(area),perimeter_m=float(perimeter),roundness=float(4*np.pi*area/perimeter**2),transverse_ratio=float(cov[0]/cov[1]),diameter_bound_m=float(np.linalg.norm(np.ptp(xy,axis=0)))))
    return rows
rows=[]
for i,j in [(1,13),(5,12),(8,14),(9,15),(10,11)]:
    if time.monotonic()-start>25:raise TimeoutError('Native chordal section measurement deadline')
    m=Mesh.load(src/f'part-{i}.npz');p,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];pm=Mesh.load(src/f'part-{j}.npz');dd,_,_=igl.point_mesh_squared_distance(p,pm.vertices,pm.faces);distance=np.sqrt(dd);origin=p[distance<distance.min()+.0004].mean(0);body=p[distance>=np.quantile(distance,.9)].mean(0);normal=body-origin;normal/=np.linalg.norm(normal);value=(p-origin)@normal;low,high=np.quantile(value,[.05,.95]);profiles=[]
    for coordinate in np.linspace(low,high,80):
        if np.min(abs(value-coordinate))<2e-8:coordinate+=5e-8
        if np.min(abs(value-coordinate))<1e-10:continue
        loops=contours(p,f,value,coordinate,origin,normal)
        if not loops:continue
        profiles.append(dict(coordinate_m=float(coordinate),loops=len(loops),minimum_roundness=min(r['roundness']for r in loops),maximum_diameter_m=max(r['diameter_bound_m']for r in loops),sections=loops))
    eligible=[r for r in profiles if r['loops']>=3 and r['minimum_roundness']>.25 and r['maximum_diameter_m']<.004]
    candidate=max(eligible,key=lambda r:r['coordinate_m'])if eligible else None
    record=dict(source=i,papillary_source=j,origin_m=origin.tolist(),normal=normal.tolist(),profiles=profiles,measured_neck_candidate=candidate,anatomical_partition_established=False)
    atomic_json(out/f'part-{i}.json',record);rows.append(dict(source=i,profiles=len(profiles),candidate_coordinate_m=candidate['coordinate_m']if candidate else None,candidate_loops=candidate['loops']if candidate else None));print(json.dumps(rows[-1]),flush=True)
atomic_json(out/'measurements.json',dict(parts=rows,geometry_modified=False,elapsed_seconds=time.monotonic()-start))
