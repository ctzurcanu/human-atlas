"""Complete the remaining full-budget field audits and septal bundle course."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
import igl
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json
from geometry_tools.geometry_instance import instance_pose
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--resume-stage');ap.add_argument('--conus-stage');args=ap.parse_args();base=ROOT/'.local-models/geometry-workbench/heart';src=base/'typed-economical-v661';out=base/args.stage;out.mkdir(exist_ok=False);started=time.monotonic();meta=copy.deepcopy(json.loads((src/'construction.json').read_text()));previous=json.loads((src/'float32-audit.json').read_text())['geometries'];rows={};resumed=json.loads((base/args.resume_stage/'float32-audit.json').read_text())['geometries']if args.resume_stage else {}
def right_bundle():
    path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};bundle=json.loads(path.with_name('cardiac-label-sites.json').read_text());sites={g['term']:g for g in bundle['integration']['termGroups']};reader=AtlasReader();pid='TA98HEART:male:myocardium:interventricular-muscular-septum';p=parts[pid];ep=parts['TA98HEART:male:epicardium'];rot,scale=instance_pose(ep);origin=np.array(ep['geometryOrigin']);local=lambda q:(np.asarray(q)-origin)@rot/scale
    m=reader.read_part(a,p,frame='assembled-metres',catalogue_path=path);rv=reader.read_part(a,parts['TA98HEART:male:myocardium:right-ventricle'],frame='assembled-metres',catalogue_path=path);spans=np.concatenate([np.arange(d['start']//3,(d['start']+d['count'])//3)for d in p['tissueDrawGroups']if d.get('neighbor')=='TA98HEART:male:endocardium']);tri=m.triangles[spans];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);faces=m.faces[spans[np.einsum('ij,ij->i',normal,rv.vertices.mean(0)-tri.mean(1))>0]];v,remap=np.unique(m.vertices,axis=0,return_inverse=True);f=remap[faces]
    landmarks=json.loads((base/'cardiac-label-landmarks-v604/report.json').read_text());st=next(x for x in landmarks['landmarks']if x['cid']==20044);ends=np.asarray(st['centroid_m'])+np.asarray(st['principal_extent_90_m'])[:,0,None]*np.asarray(st['principal_axes'])[0];target=ends[np.argmin(ends[:,1])];start=np.array(sites['A12.1.06.005']['anchor']['point']);pp,fi,projected=igl.point_mesh_squared_distance(np.array([start,target]),v,f);ids=[int(face[np.argmin(np.linalg.norm(v[face]-point,axis=1))])for face,point in zip(f[fi],projected)];edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0);w=np.linalg.norm(v[edges[:,0]]-v[edges[:,1]],axis=1);graph=coo_matrix((np.r_[w,w],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),)*2).tocsr();distance,pre=dijkstra(graph,indices=ids[0],return_predecessors=True);route=[ids[1]]
    if np.isfinite(distance[ids[1]]):
        while route[-1]!=ids[0]:route.append(int(pre[route[-1]]))
        pts=v[route[::-1]]
    else:
        # The source chart is incomplete. Bound a continuous illustrative
        # route by its native projected endpoints; do not claim a geodesic.
        pts=projected
    length=np.r_[0,np.cumsum(np.linalg.norm(np.diff(pts,axis=0),axis=1))];controls=np.array([np.interp(np.linspace(0,length[-1],9),length,pts[:,k])for k in range(3)]).T;_,face,projected=igl.point_mesh_squared_distance(controls,v,f);n=np.cross(v[f[face,1]]-v[f[face,0]],v[f[face,2]]-v[f[face,0]]);n/=np.linalg.norm(n,axis=1)[:,None];controls=local(projected-n*.00025);keep=np.r_[True,np.linalg.norm(np.diff(controls,axis=0),axis=1)>.00001];controls=controls[keep];t=np.r_[0,np.cumsum(np.linalg.norm(np.diff(controls,axis=0),axis=1))];curve=CubicSpline(t,controls,bc_type='natural')(np.linspace(0,t[-1],64));direction=np.gradient(curve,axis=0);direction/=np.linalg.norm(direction,axis=1)[:,None];verts=[];old=None;theta=np.arange(12)*np.pi/6
    for point,d in zip(curve,direction):
        x=np.cross(d,[1,0,0]if abs(d[0])<.8 else[0,1,0])if old is None else old-d*(old@d);x/=np.linalg.norm(x);y=np.cross(d,x);old=x;verts.extend(point+.00018/scale*(np.cos(theta)[:,None]*x+np.sin(theta)[:,None]*y))
    ff=[]
    for i in range(63):
        for j in range(12):
            aa=i*12+j;bb=i*12+(j+1)%12;cc=(i+1)*12+(j+1)%12;dd=(i+1)*12+j;ff.extend([[aa,bb,cc],[aa,cc,dd]])
    for j in range(1,11):ff.extend([[0,j+1,j],[63*12,63*12+j,63*12+j+1]])
    vv=np.asarray(verts).astype('float32').astype(float);ff=np.array(ff);volume=np.einsum('ij,ij->i',vv[ff[:,0]],np.cross(vv[ff[:,1]],vv[ff[:,2]])).sum()
    if volume<0:ff=ff[:,::-1]
    return Mesh(vv,ff,'heart-native-local-metres'),pid
for term,g in meta['geometries'].items():
    m=Mesh.load(src/g['file']);evidence=previous[term]
    if term in resumed and resumed[term]['passed']:
        m=Mesh.load(base/args.resume_stage/g['file']);evidence=resumed[term];assert m.geometry_sha==evidence['geometry_sha256']
    if term=='A12.1.02.008'and args.conus_stage:
        source=base/args.conus_stage;ev=json.loads((source/'float32-audit.json').read_text())['geometries'][term];assert ev['passed'];m=Mesh.load(source/g['file']);assert m.geometry_sha==ev['geometry_sha256'];evidence=ev
    if not evidence['passed']:
        if term=='A12.1.06.006':
            m,pid=right_bundle();g.update(carrierPartId=pid,method='finite continuous right-bundle reference bounded by native RV-facing septal landmarks; source chart gaps and microscopic containment remain provisional',nativeSurfaceRoute=False)
        else:m=Mesh.load(base/'typed-field-audited-v660'/g['file'])
        evidence=audit(m,require_closed=True,max_pairs=1000000,max_tests=4000000)
    m.save(out/g['file']);g.update(vertices=len(m.vertices),triangles=len(m.faces));rows[term]=evidence;atomic_json(out/'float32-audit.json',dict(complete=False,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True);print(term,evidence['passed'],evidence['failures'],len(m.faces),'triangles',round(time.monotonic()-started,2),'s',flush=True)
atomic_json(out/'construction.json',meta);atomic_json(out/'float32-audit.json',dict(complete=True,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True);failed=[t for t,e in rows.items()if not e['passed']];print('Complete',len(rows),'failed',failed,flush=True)
if failed:raise SystemExit(2)
