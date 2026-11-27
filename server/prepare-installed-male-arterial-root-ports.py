"""Extract actual unchanged arterial mouth evidence in the installed heart pose."""
from pathlib import Path
import json,sys,time
from collections import defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.audit import topology

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-installed-arterial-root-ports-v1';D.mkdir(exist_ok=False)
cp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cp.read_text());assert digest_file(cp)=='a0d7b46e033d145617a8937e1aff489db90b7a2808aa76ec90192f155ad95518'
reader=AtlasReader();rows={p['id']:p for p in a['parts']};physical={x['part']for x in json.loads((ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/preview-proof.json').read_text())['physical_tissues']};triangles=[];owner=[];sourceface=[]
for pid in sorted(physical):
 p=rows[pid];groups=[np.arange(g['start']//3,(g['start']+g['count'])//3)for g in p.get('tissueDrawGroups',[])if g.get('neighbor')not in physical]
 if not groups:continue
 keep=np.concatenate(groups);m=reader.read_part(a,p,frame='model:local-ta98:assembled-metres');triangles.append(m.triangles[keep]);owner.extend([pid]*len(keep));sourceface.extend(keep)
tri=np.concatenate(triangles);v,iv=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);f=iv.reshape(-1,3);m=Mesh(v,f,'model:local-ta98:assembled-metres');m.save(D/'actual-host-boundary.npz');np.savez_compressed(D/'actual-host-provenance.npz',owner=np.asarray(owner),source_face=np.asarray(sourceface));owner=np.asarray(owner);sourceface=np.asarray(sourceface)
edges=defaultdict(list)
for fi,face in enumerate(f):
 for x,y in zip(face,np.roll(face,-1)):edges[tuple(sorted((int(x),int(y))))].append(fi)
chosen={e:ii for e,ii in edges.items()if len(ii)==2 and len(set(owner[ii]))==2 and 'TA98HEART:male:endocardium'in owner[ii]};graph=defaultdict(set)
for x,y in chosen:graph[x].add(y);graph[y].add(x)
unseen=set(graph);curves=[]
while unseen:
 seed=min(unseen);todo=[seed];seen={seed};unseen.remove(seed)
 while todo:
  at=todo.pop()
  for nxt in graph[at]:
   if nxt in unseen:unseen.remove(nxt);seen.add(nxt);todo.append(nxt)
 closed=all(len(graph[k])==2 for k in seen);ordered=[]
 if closed:
  at=seed;prev=None
  while not ordered or at!=seed:
   ordered.append(at);nxt=min(graph[at]-({prev}if prev is not None else set()));prev,at=at,nxt
 else:ordered=sorted(seen)
 pts=v[ordered];ef=[(e,ii)for e,ii in chosen.items()if e[0]in seen]
 curves.append(dict(closed=closed,controls=len(ordered),vertices=ordered,points_m=pts.tolist(),center_m=pts.mean(0).tolist(),span_mm=(np.ptp(pts,axis=0)*1000).tolist(),parts=sorted({str(owner[k])for _,ii in ef for k in ii}),source_edges=[dict(vertices=list(e),owners=[str(owner[k])for k in ii],source_faces=sourceface[ii].tolist())for e,ii in ef]))
roots=[]
for pid,cusps in [('ZA:Ascending aorta',['ZA:Right coronary leaflet','ZA:Left coronary leaflet','ZA:Non-coronary leaflet']),('ZA:Pulmonary trunk',['ZA:Right semilunar leaflet of pulmonary valve','ZA:Left semilunar leaflet of pulmonary valve','ZA:Anterior semilunar leaflet of pulmonary valve'])]:
 vessel=reader.read_part(a,rows[pid],frame=m.frame);vv,iv=np.unique(vessel.vertices,axis=0,return_inverse=True);loops=boundary_loops(iv[vessel.faces]);loops.sort(key=lambda r:vv[r,1].mean());assert len(loops)==2;prox=vv[loops[0]];center=prox.mean(0);cuspcenter=np.mean([reader.read_part(a,rows[k],frame=m.frame).vertices.mean(0)for k in cusps],axis=0)
 ranked=sorted(range(len(curves)),key=lambda i:np.linalg.norm(np.asarray(curves[i]['points_m'])-center,axis=1).min())[:8]
 roots.append(dict(part=pid,proximal_native_ring=prox.tolist(),proximal_native_vertices=loops[0].tolist(),center_m=center.tolist(),cusp_seed_m=cuspcenter.tolist(),near_curves=[dict(index=i,minimum_distance_mm=float(np.linalg.norm(np.asarray(curves[i]['points_m'])-center,axis=1).min()*1000),**curves[i])for i in ranked]))
report=dict(schema='installed-male-arterial-root-mouths/v1',source_catalogue_sha256=digest_file(cp),heart_pose_preserved=True,vessels_fixed=True,host_topology=topology(m),curves=curves,roots=roots,elapsed_seconds=time.monotonic()-T);atomic_json(D/'ports.json',report)
print(json.dumps(dict(source_sha256=report['source_catalogue_sha256'],elapsed_seconds=report['elapsed_seconds'],curves=len(curves),roots=[dict(part=r['part'],center=r['center_m'],curves=[{k:c[k]for k in ['index','minimum_distance_mm','closed','controls','center_m','span_mm','parts']}for c in r['near_curves']])for r in roots])))
