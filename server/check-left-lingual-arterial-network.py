"""Actual parent aperture, guide patency and oral/tissue-entry checks."""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import audit

args=argparse.ArgumentParser();args.add_argument('--version',default='v3');args.add_argument('--encoded',action='store_true');args.add_argument('--side',choices=['left','right'],default='left');args.add_argument('--route-version');args.add_argument('--layout',default='owned-bank-v2');args=args.parse_args()
T=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/('male-'+args.side+'-lingual-finite-network-'+args.version);O=D/('independent-check-emitted-v1'if args.encoded else'independent-check-v1');O.mkdir(exist_ok=False);P=B/('male-'+args.side+'-lingual-routes-'+(args.route_version or args.version));plan=json.loads((P/'routes.json').read_text());inputs=json.loads((B/'male-lingual-source-plan-v1/inputs.json').read_text());sources={r['role']:Mesh.load(ROOT/r['mesh'])for r in inputs['fixed_sources']};wall=Mesh.load(D/'finite-wall.npz');blood=Mesh.load(D/'closed-blood-classifier.npz');outer=Mesh.load(D/'outer-sheet.npz');data=dict(np.load(D/'surfaces-and-ancestry.npz'));portdata=json.loads((D/'ports.json').read_text());component_names=portdata['component_names'];component_terms=portdata['component_terms']
if args.encoded:
    layout=json.loads((D/args.layout/'tissue.json').read_text());assert layout['passed'];origin=np.array(layout['parts']['lingual']['geometryOrigin'])
    def encoded(mesh):return Mesh((mesh.vertices-origin).astype('<f4').astype(float)+origin,mesh.faces,mesh.frame)
    wall,blood,outer=map(encoded,[wall,blood,outer]);data['vertices']=(data['vertices']-origin).astype('<f4').astype(float)+origin
newfaces=np.r_[data['outer_faces'][data['outer_component']>0],data['inner_faces'][data['inner_component']>0],*[data[f'annulus_{i}_faces']for i in range(2,7)]];new=Mesh(data['vertices'],newfaces,wall.frame);new.save(O/'actual-new-arterial-wall.npz')
paths=[];segments=[]
for component,name in enumerate(component_names[1:],1):
    z=np.load(D/(name+'-fields.npz'));points=z['guide_points'];paths.append(points);segments.extend((component,i,a,b)for i,(a,b)in enumerate(zip(points[:-1],points[1:])))
parent_centers=np.array(plan['fixed_native_carotid']['source_station_components'][0]['centres_m']);segments.extend((0,i,a,b)for i,(a,b)in enumerate(zip(parent_centers[:-1],parent_centers[1:])))
mid=np.array([(a+b)/2 for _,_,a,b in segments]);winding=igl.winding_number(blood.vertices,blood.faces,mid);tri=wall.triangles;lo=tri.min(1);hi=tri.max(1);hits=[]
for component,interval,a,b in segments:
    direction=b-a;sel=np.flatnonzero(np.all(hi>=np.minimum(a,b),axis=1)&np.all(lo<=np.maximum(a,b),axis=1));t=tri[sel];e1=t[:,1]-t[:,0];e2=t[:,2]-t[:,0];h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h);good=abs(det)>1e-24;sel=sel[good];t=t[good];e1=e1[good];e2=e2[good];h=h[good];det=det[good];s=a-t[:,0];u=np.einsum('ij,ij->i',s,h)/det;k=np.cross(s,e1);v=np.einsum('j,ij->i',direction,k)/det;mu=np.einsum('ij,ij->i',e2,k)/det;wrong=np.flatnonzero((u>=0)&(v>=0)&(u+v<=1)&(mu>=0)&(mu<=1));hits.extend(dict(component=component,interval=interval,wall_face=int(sel[i]),fraction=float(mu[i]))for i in wrong)
atomic_npz(O/'actual-guide-segments.npz',segments=np.array([[a,b]for _,_,a,b in segments]),component=np.array([c for c,_,_,_ in segments]),source_interval=np.array([i for _,i,_,_ in segments]))
targets={'suprahyoid':sources['hyoglossus-'+args.side],'dorsal-lingual-0':Mesh.load(P/'inner-tongue-fixed-reference.npz'),'dorsal-lingual-1':Mesh.load(P/'inner-tongue-fixed-reference.npz'),'sublingual':sources['sublingual-gland-'+args.side],'deep-lingual':Mesh.load(P/'inner-tongue-fixed-reference.npz')};terminals=[]
for port in portdata['ports'][2:]:
    component=port['component'];name=component_names[component];target=targets[name];points=wall.vertices[port['outer_indices']];center=points.mean(0);center_distance=float(np.sqrt(igl.point_mesh_squared_distance(center[None],target.vertices,target.faces)[0][0]));radius=float(np.linalg.norm(points-center,axis=1).max());inside=abs(igl.winding_number(target.vertices,target.faces,points));ball_clearance=center_distance-radius;terminals.append(dict(name=name,entire_annular_terminal_inside_recipient=bool(inside.min()>.99 and ball_clearance>0),minimum_source_vertex_winding=float(inside.min()),whole_annulus_clearance_lower_m=ball_clearance,source_recipient_geometry_sha256=target.geometry_sha))
contexts=[]
for role in ['tongue-mucosa','greater-horn-'+args.side,'hyoglossus-'+args.side,'sublingual-gland-'+args.side]:
    mesh=sources[role];contacts,summary=pair_contacts(new,mesh,0.,max_pairs=100000,max_tests=1000000);assert summary['complete_within_pair_budget'];atomic_npz(O/(role+'-contacts.npz'),**contacts);contexts.append(dict(role=role,geometry_sha256=mesh.geometry_sha,summary=summary,interpretation='Muscle/gland entry is intentional only for assigned branches; exposed mucosal/bony crossings are forbidden'))
whole=Mesh.load(P/'whole-tongue-fixed-reference.npz');used=np.unique(new.faces);points=new.vertices[used];inside=abs(igl.winding_number(whole.vertices,whole.faces,points))>.99
def smooth(x):
    t=np.clip(x,0.,1.);return t*t*(3-2*t)
reflection=1.6265-.0025*smooth((points[:,2]-.055)/.017)-.002*smooth(abs(points[:,0])/.016);oral_box=(abs(points[:,0])<.0181)&(points[:,2]>=.0234)&(points[:,2]<=.0838);oral_above=points[:,1]>reflection+.0001;outside=oral_box&oral_above&~inside
atomic_npz(O/'oral-reference-classification.npz',source_vertex_ids=used,positions=points,inside_closed_whole_tongue=inside,within_oral_reference_domain=oral_box,above_reflection=oral_above,outside_tissue_in_oral_reference_domain=outside)
parent=sources['external-carotid-'+args.side];literal={tuple(sorted(map(tuple,t)))for t in outer.triangles};retained=np.array([tuple(sorted(map(tuple,t)))in literal for t in parent.triangles]);patency=bool(not hits and winding.min()>.99);mucosa=next(r for r in contexts if r['role']=='tongue-mucosa')['summary'];bone=next(r for r in contexts if r['role']=='greater-horn-'+args.side)['summary'];passed=bool(patency and all(r['entire_annular_terminal_inside_recipient']for r in terminals)and not mucosa['crossing_face_pairs']and not bone['crossing_face_pairs']and not outside.any())
report=dict(schema='left-lingual-source-tissue-patency-check/v1',whole_guide_segments=len(segments),wall_segment_hits=hits,minimum_blood_winding=float(winding.min()),patency_pass=patency,terminal_tissue_interfaces=terminals,context=contexts,
    native_ECA_outer_faces=len(parent.faces),literal_native_ECA_outer_faces=int(retained.sum()),actual_parent_ostium_source_faces_clipped=np.flatnonzero(~retained).tolist(),
    oral_reference_outside_tissue_vertices=int(outside.sum()),oral_domain_scope='Native closed tongue plus pre-existing approximate FOM reflection. This is not a full independently reconstructed oral-cavity volume.',
    numerical_route_gate_pass=passed,anatomical_acceptance=False,installable=False,elapsed_seconds=time.monotonic()-T)
atomic_json(O/'evidence.json',report)
if args.encoded:
    print(json.dumps(dict(actual_float32_wall_and_ports_checked=True,segments=len(segments),patency=patency,mucosal_crossings=mucosa['crossing_face_pairs'],hyoid_crossings=bone['crossing_face_pairs'],oral_outside_tissue_vertices=int(outside.sum()),passed=passed,elapsed_seconds=report['elapsed_seconds'])));sys.exit(0)
entries=[dict(mesh=str(D/'finite-wall.npz'),name='Left five-term lingual arterial network; native ECA position fixed',color=[.88,.15,.12,1])]
for role,color in [('tongue-mucosa',[.55,.64,.84,.22]),('hyoglossus-'+args.side,[.8,.62,.24,.3]),('greater-horn-'+args.side,[.85,.83,.72,1]),('sublingual-gland-'+args.side,[.77,.47,.65,.35])]:entries.append(dict(mesh=str(ROOT/next(r['mesh']for r in inputs['fixed_sources']if r['role']==role)),name=role+' fixed source',color=color))
root=np.array(plan['fixed_native_carotid']['source_station_components'][0]['centres_m'][plan['source_origin_station']]);atomic_json(D/'review-spec.json',dict(frame=wall.frame,resolution=1100,wall_seconds=60,meshes=entries,views=[dict(name='Left native ECA to lingual family and fixed tissue recipients',direction=[0,.1,1]),dict(name='Medial hyoglossus route and deep/sublingual branches',direction=[1,.1,0]),dict(name='Uncapped true ECA branch opening',direction=[0,1,0],section=dict(normal=[0,1,0],point=root.tolist(),clear_outer=True))]));print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'],segments=len(segments),patency=patency,terminal_interfaces=terminals,mucosal_crossings=mucosa['crossing_face_pairs'],hyoid_crossings=bone['crossing_face_pairs'],oral_outside_tissue_vertices=int(outside.sum()),passed=passed)))
