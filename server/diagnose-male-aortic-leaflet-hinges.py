"""Read-only actual cusp boundary/host/coaptation evidence; no inferred annulus."""
from pathlib import Path
import sys,time,json
from collections import defaultdict
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-native-hinge-diagnostic-v1';D.mkdir(exist_ok=False);cp=ROOT/'.local-models/geometry-workbench/jobs/b9c8595514cdc1b8671263ec425dbf7bbad8f059c484638ae4b24211c0b2d9dc/prepared-2/proposal.json';a=json.loads(cp.read_text());pmap={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';ids=['ZA:Right coronary leaflet','ZA:Left coronary leaflet','ZA:Non-coronary leaflet'];meshes={pid:reader.read_part(a,pmap[pid],frame=frame)for pid in ids};hosts={pid:reader.read_part(a,pmap[pid],frame=frame)for pid in ['TA98HEART:male:endocardium','TA98HEART:male:epicardium','TA98HEART:male:myocardium:left-ventricle']};records=[];review=[];all_nearest={}
for j,pid in enumerate(ids):
 m=meshes[pid];m.save(D/(str(j)+'.npz'));top=topology(m);tri=m.triangles;fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);fn/=np.linalg.norm(fn,axis=1)[:,None];edges=defaultdict(list)
 for fi,face in enumerate(m.faces):
  for aa,bb in zip(face,np.roll(face,-1)):edges[tuple(sorted([int(aa),int(bb)]))].append(fi)
 ridge=[];angles=[]
 for edge,faces in edges.items():
  if len(faces)!=2:continue
  angle=float(np.degrees(np.arccos(np.clip(fn[faces[0]]@fn[faces[1]],-1,1))));angles.append(angle)
  if angle>=35:ridge.append(edge)
 graph=defaultdict(set)
 for aa,bb in ridge:graph[aa].add(bb);graph[bb].add(aa)
 unseen=set(graph);components=[]
 while unseen:
  seed=min(unseen);stack=[seed];seen={seed};unseen.remove(seed)
  while stack:
   at=stack.pop()
   for nxt in graph[at]:
    if nxt in unseen:unseen.remove(nxt);seen.add(nxt);stack.append(nxt)
  components.append(dict(controls=len(seen),edge_count=sum(len(graph[k])for k in seen)//2,closed_degree_two=all(len(graph[k])==2 for k in seen),branch_degrees=sorted(set(len(graph[k])for k in seen))))
 distances=[];values={}
 for target,q in {**hosts,**{k:v for k,v in meshes.items()if k!=pid}}.items():
  d,face,witness=igl.point_mesh_squared_distance(m.vertices,q.vertices,q.faces);mm=np.sqrt(d)*1000;_,contacts=pair_contacts(m,q,0.,max_pairs=100000,max_tests=1000000);assert contacts['complete_within_pair_budget'];distances.append(dict(target=target,minimum_control_distance_mm=float(mm.min()),median_control_distance_mm=float(np.median(mm)),maximum_control_distance_mm=float(mm.max()),controls_within_0_25mm=int((mm<=.25).sum()),controls_within_0_01mm=int((mm<=.01).sum()),actual_crossing_face_pairs=contacts['crossing_face_pairs']));values[target.replace(':','_')]=mm
 atomic_npz(D/(str(j)+'-source-edge-data.npz'),vertices=m.vertices,faces=m.faces,ridge_edges=np.asarray(ridge,int).reshape(-1,2),frame=np.asarray(frame),**values);record=dict(part=pid,source_geometry_sha256=m.geometry_sha,source_bank_sha256=m.source['bank_sha256'],source_vertices=len(m.vertices),source_faces=len(m.faces),topology=top,declared_hinge_curve=False,declared_coaptation_edge=False,declared_parent_annular_interface=False,actual_open_boundary_edges=top['open_edges'],sharp_ridge_threshold_degrees=35,maximum_native_dihedral_degrees=max(angles),unclassified_sharp_ridge_components=components,distances=distances);records.append(record);review.append(dict(mesh=str(D/(str(j)+'.npz')),name=pid+' actual closed source; hinge/free edge unclassified',color=[[.2,.75,.38,1],[.24,.52,.9,1],[.92,.68,.2,1]][j]))
report=dict(schema='male-aortic-actual-cusp-hinge-evidence/v1',source_catalogue_sha256=digest_file(cp),leaflets=records,all_source_faces_vertices_and_normal_words_literal=True,hinge_annulus_constructed=False,convex_hull_used=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(D/'evidence.json',report);center=np.mean(np.concatenate([m.vertices for m in meshes.values()]),axis=0);atomic_json(D/'review-spec.json',dict(frame=frame,meshes=review,resolution=1100,wall_seconds=45,views=[dict(name='actual three native closed leaflets, oblique',direction=[.4,.7,1]),dict(name='actual leaflet side profile',direction=[-1,.1,.4]),dict(name='uncapped leaflet middle section',direction=[0,0,1],section=dict(normal=[0,0,1],point=center.tolist(),clear_outer=True))]));print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'],hinge_annulus_constructed=False,leaflets=[dict(part=r['part'],open_edges=r['actual_open_boundary_edges'],components=len(r['topology']['components']),ridge_components=r['unclassified_sharp_ridge_components'],contacts=r['distances'])for r in records])))
