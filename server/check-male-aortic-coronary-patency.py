"""Independent source-guide and retained far-course checks on the ostia candidate."""
from pathlib import Path
import sys,json,time
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh
from geometry_tools.common import ROOT,atomic_json,atomic_npz
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.audit import audit

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-coronary-ostia-v2';wall=Mesh.load(D/'finite-wall.npz');blood=Mesh.load(D/'closed-blood-classifier.npz');plan=json.loads((D/'construction-plan.json').read_text());sourcepath=ROOT/'.local-models/geometry-workbench/heart/male-native-aortic-source-defects-v2';parent=Mesh.load(ROOT/'.local-models/geometry-workbench/heart/male-aortic-section-transport-v1/candidate-outer.npz');left=Mesh.load(sourcepath/'Left-coronary-artery.npz');right=Mesh.load(sourcepath/'Right-coronary-artery.npz');rcabranch=max(recover_station_rings(right.vertices,right.faces),key=lambda b:max(right.vertices[b['rings'][0],1].mean(),right.vertices[b['rings'][-1],1].mean()));selected=np.concatenate(rcabranch['rings'][:13]);right=Mesh(right.vertices,right.faces[np.all(np.isin(right.faces,selected),axis=1)],right.frame);sources=[parent,left,right];guides=[];sourceproof=[];outer=Mesh.load(D/'outer-sheet.npz');literal={tuple(sorted(map(tuple,t)))for t in outer.triangles}
for i,m in enumerate(sources):
 branch=recover_station_rings(m.vertices,m.faces)[0];rings=list(branch['rings'])
 if i in [1,2]and m.vertices[rings[0],1].mean()<m.vertices[rings[-1],1].mean():rings=rings[::-1]
 centers=np.array([m.vertices[r].mean(0)for r in rings]);guides.append(centers);station=np.full(len(m.vertices),-1,int)
 for j,r in enumerate(rings):station[r]=j
 retained=np.array([tuple(sorted(map(tuple,t)))in literal for t in m.triangles]);modified_intervals=sorted(set(station[m.faces[~retained]].min(1).tolist()));sourceproof.append(dict(component=i,native_source_faces=len(m.faces),literal_retained_source_faces=int(retained.sum()),modified_source_intervals=modified_intervals,furthest_modified_source_station=max(modified_intervals,default=-1)))
segments=[]
for i,c in enumerate(guides):segments.extend((i,j,p,q)for j,(p,q)in enumerate(zip(c[:-1],c[1:])))
for sleeve in plan['root_sleeves']:
 i=sleeve['child'];target=np.array(sleeve['buried_classifier_center']);root=np.array(sleeve['native_source_root']);nearest=guides[0][np.argmin(np.linalg.norm(guides[0]-target,axis=1))];segments.extend([(i,-2,nearest,target),(i,-1,target,root)])
midpoints=np.array([(p+q)/2 for _,_,p,q in segments]);winding=igl.winding_number(blood.vertices,blood.faces,midpoints);assert np.min(winding)>.99;tri=wall.triangles;lo=tri.min(1);hi=tri.max(1);hits=[]
for i,j,p,q in segments:
 direction=q-p;sel=np.flatnonzero(np.all(hi>=np.minimum(p,q),axis=1)&np.all(lo<=np.maximum(p,q),axis=1));t=tri[sel];e1=t[:,1]-t[:,0];e2=t[:,2]-t[:,0];h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h);good=abs(det)>1e-24;sel=sel[good];t=t[good];e1=e1[good];e2=e2[good];h=h[good];det=det[good];s=p-t[:,0];u=np.einsum('ij,ij->i',s,h)/det;k=np.cross(s,e1);v=np.einsum('j,ij->i',direction,k)/det;mu=np.einsum('ij,ij->i',e2,k)/det;wrong=np.flatnonzero((u>=0)&(v>=0)&(u+v<=1)&(mu>=0)&(mu<=1));hits.extend(dict(component=i,source_interval=j,wall_face=int(sel[x]),fraction=float(mu[x]))for x in wrong)
assert not hits,hits;bloodcheck=audit(blood,require_closed=True,expected_components=1,expected_genus=[0]);assert bloodcheck['passed'];atomic_json(D/'blood-classifier-audit.json',bloodcheck);report=dict(schema='male-aortic-coronary-independent-patency/v1',source_bound_construction=plan['source_catalogue_sha256'],whole_native_and_root_connector_segments=len(segments),actual_wall_segment_crossings=len(hits),midpoints_inside_one_blood_volume=True,minimum_blood_winding=float(winding.min()),blood_volume_components=1,blood_volume_genus=0,native_outer_retention=sourceproof,heart_and_arch_joins_verified=False,wall_thickness_clinically_verified=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(D/'independent-patency.json',report);atomic_npz(D/'native-and-connector-guides.npz',segments=np.array([[p,q]for _,_,p,q in segments]),component=np.array([i for i,_,_,_ in segments]),source_interval=np.array([j for _,j,_,_ in segments]));print(json.dumps(report))
