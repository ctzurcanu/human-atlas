"""One bounded, uninstalled three-sinus root and communicating coronary wall.

The basal cardiac mouth and leaflet hinge ownership stay pending. Fixed heart,
distal native AA profiles and all native coronary positions are never fitted.
"""
from pathlib import Path
import sys,time,json,hashlib
import numpy as np
import manifold3d as md
from scipy.spatial.transform import Rotation,Slerp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,CommonPositionBank
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim,continuity_audit,loft_faces,finite_wall_cells
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.audit import audit
from geometry_tools.contact import self_contact_summary

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-three-sinus-root-prototype-v1b';D.mkdir(exist_ok=False)
cp=ROOT/'.local-models/geometry-workbench/jobs/b9c8595514cdc1b8671263ec425dbf7bbad8f059c484638ae4b24211c0b2d9dc/prepared-2/proposal.json';source_sha=digest_file(cp);assert source_sha=='9cf8b30456620f7b7bd6f52311ba818fc94f935beee24e6e07cb953d14106f64';a=json.loads(cp.read_text());rows={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
aa=reader.read_part(a,rows['ZA:Ascending aorta'],frame=frame);lca=reader.read_part(a,rows['ZA:Left coronary artery'],frame=frame);rca=reader.read_part(a,rows['ZA:Right coronary artery'],frame=frame)
branch=recover_station_rings(aa.vertices,aa.faces)[0];native_rings,_=align_station_correspondence(aa.vertices,aa.faces,branch['rings']);native_centers=np.array([aa.vertices[r].mean(0)for r in native_rings]);cut=14
native_used=np.concatenate(native_rings[cut:]);native_faces=aa.faces[np.all(np.isin(aa.faces,native_used),axis=1)];cut_loop=next(loop for loop in boundary_loops(native_faces)if set(loop)==set(native_rings[cut]));mouth=cut_loop;end_center=aa.vertices[mouth].mean(0);q=aa.vertices[mouth]-end_center
normal=np.cross(q,np.roll(q,-1,axis=0)).sum(0);normal/=np.linalg.norm(normal);forward=native_centers[15]-native_centers[14]
assert normal@forward>0
end_u=q[0]-normal*(q[0]@normal);end_u/=np.linalg.norm(end_u);end_v=np.cross(normal,end_u);end_frame=np.c_[end_u,end_v,normal];beta=np.arctan2(q@end_v,q@end_u)
hinge=json.loads((ROOT/'.local-models/geometry-workbench/heart/male-aortic-provisional-source-hinges-v1/evidence.json').read_text());root_center=np.array(hinge['root_center_m']);root_axis=np.array(hinge['root_axis']);start_u=end_u-root_axis*(end_u@root_axis);start_u/=np.linalg.norm(start_u);start_frame=np.c_[start_u,np.cross(root_axis,start_u),root_axis]
cusp_ids=['ZA:Right coronary leaflet','ZA:Left coronary leaflet','ZA:Non-coronary leaflet'];cusps=[reader.read_part(a,rows[pid],frame=frame)for pid in cusp_ids];cusp_centers=np.array([m.vertices.mean(0)for m in cusps]);basal_height=min(float(((m.vertices-root_center)@root_axis).min())for m in cusps)-.001;start_center=root_center+root_axis*basal_height
K=18;parameters=np.linspace(0,1,K+1);frames=Slerp([0,1],Rotation.from_matrix(np.array([start_frame,end_frame])))(parameters).as_matrix();centers=(1-parameters[:,None])*start_center+parameters[:,None]*end_center
reference_radius=float(np.linalg.norm(q,axis=1).mean());local_angles=[];amplitudes=[]
for cusp in cusps:
    delta=cusp.vertices.mean(0)-root_center;local_angles.append(float(np.arctan2(delta@start_frame[:,1],delta@start_frame[:,0])))
    # Source-derived directional extent sets a sinus bulge, not a copied mesh.
    radius=np.linalg.norm(cusp.vertices-root_center-((cusp.vertices-root_center)@root_axis)[:,None]*root_axis,axis=1)
    amplitudes.append(max(.001,float(radius.max()+.001-reference_radius)))
plan=dict(schema='male-three-sinus-root-prototype-plan/v1',source_catalogue_sha256=source_sha,authoritative_inventory='TA98-CIRCULATION.md A12.2.03 descendants',
    parts=[dict(term='A12.2.03.001',representation='one paired finite root/ascending wall'),dict(term='A12.2.03.002',representation='three sinus blood spaces, not solids'),dict(term='A12.2.03.003',representation='supravalvular boundary curve'),dict(term='A12.2.03.004',representation='aortic bulb wall region'),dict(term='A12.2.03.101',representation='existing fixed RCA root mouth'),dict(term='A12.2.03.201',representation='existing fixed LCA root mouth')],
    interleaflet_support='three finite-wall angular sectors sharing coordinates and chart boundaries; separate cusp endpoints are not falsely welded',
    source_cusp_attachment_status='three independent source-qualified peripheral references; anatomical U attachment ownership pending',
    retained_native_aortic_stations=list(range(14,25)),native_distal_and_coronary_positions_literal=True,heart_and_skeleton_literal=True,
    basal_cardiac_host_join='pending actual matched opening',prototype_thickness_m=.0009,thickness_status='illustrative, not established histology',
    sinus_amplitude_m=amplitudes,root_course='one new straight center course with proper gradual section rotation; no installed object translated',installable=False,anatomical_acceptance=False)
atomic_json(D/'construction-plan.json',plan)

ring_points=[];region_tags=[]
for j,(t,c,basis)in enumerate(zip(parameters,centers,frames)):
    shape=np.zeros(12);tags=[]
    for k,angle in enumerate(beta):
        delta=np.angle(np.exp(1j*(angle-np.asarray(local_angles))));weights=np.exp(-.5*(delta/.46)**2);shape[k]=sum(weights*np.array(amplitudes));tags.append(int(np.argmax(weights)))
    gain=np.sin(np.pi*t)**2;radius=np.linalg.norm(q,axis=1)+gain*shape
    radial=np.cos(beta)[:,None]*basis[:,0]+np.sin(beta)[:,None]*basis[:,1];ring_points.append(c+radius[:,None]*radial);region_tags.append(tags)
ring_points[-1]=aa.vertices[mouth].copy()
v=np.r_[aa.vertices,np.concatenate(ring_points[:-1])];root_rings=[np.arange(len(aa.vertices)+12*j,len(aa.vertices)+12*(j+1))for j in range(K)]+[mouth]
root_faces=np.array([[root_rings[int(index)//12][int(index)%12]for index in face]for face in loft_faces(K+1,12)],int)
parent=Mesh(v,np.r_[root_faces,native_faces],frame);parent.save(D/'root-and-native-distal-outer.npz')
all_rings=root_rings+[r for r in native_rings[15:]];V=[v];F=[parent.faces];R=[all_rings];C=[np.array([v[r].mean(0)for r in all_rings])]
for source in [lca,rca]:
    branch=max(recover_station_rings(source.vertices,source.faces),key=lambda b:max(source.vertices[b['rings'][0],1].mean(),source.vertices[b['rings'][-1],1].mean()));rings=list(branch['rings'])
    if source.vertices[rings[0],1].mean()<source.vertices[rings[-1],1].mean():rings=rings[::-1]
    if source is rca:rings=rings[:13]
    used=np.concatenate(rings);faces=source.faces[np.all(np.isin(source.faces,used),axis=1)];V.append(source.vertices.copy());F.append(faces);R.append(rings);C.append(np.array([source.vertices[r].mean(0)for r in rings]))
thickness=[.0009,.0003,.0003];ORIGIN=np.array([0.,1.4,.03125]);external=[(0,0),(0,len(R[0])-1),(1,len(R[1])-1),(2,len(R[2])-1)];targets={}
for i in [1,2]:
    station=int(np.argmin(np.linalg.norm(C[0]-C[i][0],axis=1)));targets[i]=C[0][station]+.60*(C[i][0]-C[0][station])
inner=v.copy()
for ring,c in zip(R[0],C[0]):
    delta=v[ring]-c;inner[ring]=v[ring]-.0009*delta/np.linalg.norm(delta,axis=1)[:,None]
fields=np.array([v,inner]);atomic_npz(D/'root-fields.npz',fields=fields,faces=parent.faces,frame=np.array(frame),root_faces=root_faces,native_distal_faces=native_faces,
    rings=np.array(all_rings),region_tags=np.array(region_tags),root_centers=centers,root_frames=frames)
partial_faces=[parent.faces,inner.shape[0]+parent.faces[:,::-1]]
for ring in boundary_loops(parent.faces):partial_faces.append(annular_rim(ring,ring+len(v)))
partial=Mesh(np.r_[v,inner],np.concatenate(partial_faces),frame);partial.save(D/'partial-root-paired-wall.npz')
prisms=finite_wall_cells(fields,parent.faces);outer_cross,_=self_contact_summary(parent,max_pairs=100000,max_tests=1000000)
atomic_json(D/'paired-carrier-audit.json',dict(prisms=prisms,outer_self_crossings=outer_cross))
construction_pass=False;failure=None
try:
    assert prisms['finite_wall_cells_pass'],('invalid_paired_root_prisms',prisms['invalid_counts_by_interval'])
    assert outer_cross['complete_within_pair_budget'] and outer_cross['defect_pairs']==0,('root_outer_self_crossings',outer_cross['defect_pairs'])
    # Execute the already reviewed finite-network constructor unchanged. Its
    # actual input arrays and runtime are bound in this stage's evidence.
    constructor=ROOT/'server/stage-male-aortic-coronary-ostia.py';text=constructor.read_text();body='inputs={};sheets=[];'+text.split('inputs={};sheets=[];',1)[1];exec(compile(body,str(constructor),'exec'))
    construction_pass=True
except Exception as error:
    failure=str(error);atomic_json(D/'construction-failure.json',dict(error=failure,elapsed_seconds=time.monotonic()-T,source_geometry_changed=False))

shown=D/('finite-wall.npz'if construction_pass else'partial-root-paired-wall.npz');entries=[dict(mesh=str(shown),name='UNINSTALLED three-sinus paired root prototype',color=[.78,.3,.18,.55])]
for i,cusp in enumerate(cusps):
    path=D/f'unchanged-cusp-{i}.npz';cusp.save(path);entries.append(dict(mesh=str(path),name=cusp_ids[i]+' fixed source',color=[.18,.76,.39,1]))
cut_point=root_center.tolist();atomic_json(D/'review-spec.json',dict(frame=frame,resolution=1100,wall_seconds=60,meshes=entries,
    views=[dict(name='Uninstalled three-sinus paired root and fixed cusps',direction=[.4,.7,1]),dict(name='Uncapped longitudinal root/coronary lumen cut',direction=[0,0,1],section=dict(normal=[0,0,1],point=cut_point,clear_outer=True)),dict(name='Root from below; basal cardiac attachment pending',direction=(-root_axis).tolist())]))
report=dict(schema='male-three-sinus-root-prototype/v1',source_catalogue_sha256=source_sha,actual_partial_geometry=str(D/'partial-root-paired-wall.npz'),finite_network_constructed=construction_pass,
    failure=failure,heart_and_coronary_courses_unchanged=True,native_aortic_stations_14_to_24_literal=True,
    basal_cardiac_join_pending=True,leaflet_attachment_ownership_pending=True,three_sinus_regions=True,three_interleaflet_support_regions=True,
    illustrative_thickness=True,installed=False,installable=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T)
atomic_json(D/'prototype-evidence.json',report);print(json.dumps(report))
