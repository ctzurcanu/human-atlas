"""Independent full-segment patency, literal far source and fixed context."""
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_station_topology import recover_station_rings

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-three-sinus-root-prototype-v1b';O=D/'independent-check-v3';O.mkdir(exist_ok=False)
cp=ROOT/'.local-models/geometry-workbench/jobs/b9c8595514cdc1b8671263ec425dbf7bbad8f059c484638ae4b24211c0b2d9dc/prepared-2/proposal.json';a=json.loads(cp.read_text());rows={p['id']:p for p in a['parts']};reader=AtlasReader();wall=Mesh.load(D/'finite-wall.npz');outer=Mesh.load(D/'outer-sheet.npz');blood=Mesh.load(D/'closed-blood-classifier.npz');root=np.load(D/'root-fields.npz');guides=[root['fields'][0][root['rings']].mean(1)];source_rows=[]
for pid in ['ZA:Left coronary artery','ZA:Right coronary artery']:
    mesh=reader.read_part(a,rows[pid],frame=wall.frame);branch=max(recover_station_rings(mesh.vertices,mesh.faces),key=lambda b:max(mesh.vertices[b['rings'][0],1].mean(),mesh.vertices[b['rings'][-1],1].mean()));rings=list(branch['rings'])
    if mesh.vertices[rings[0],1].mean()<mesh.vertices[rings[-1],1].mean():rings=rings[::-1]
    if 'Right coronary' in pid:rings=rings[:13]
    guides.append(np.array([mesh.vertices[r].mean(0)for r in rings]));source_rows.append((pid,mesh,rings))
segments=[(i,j,p,q)for i,curve in enumerate(guides)for j,(p,q)in enumerate(zip(curve[:-1],curve[1:]))]
for i in [1,2]:
    station=int(np.argmin(np.linalg.norm(guides[0]-guides[i][0],axis=1)));target=guides[0][station]+.60*(guides[i][0]-guides[0][station]);segments.extend([(i,-2,guides[0][station],target),(i,-1,target,guides[i][0])])
points=np.array([(p+q)/2 for _,_,p,q in segments]);winding=igl.winding_number(blood.vertices,blood.faces,points);tri=wall.triangles;lo=tri.min(1);hi=tri.max(1);hits=[]
for component,interval,p,q in segments:
    direction=q-p;sel=np.flatnonzero(np.all(hi>=np.minimum(p,q),axis=1)&np.all(lo<=np.maximum(p,q),axis=1));t=tri[sel];e1=t[:,1]-t[:,0];e2=t[:,2]-t[:,0];h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h);good=abs(det)>1e-24;sel=sel[good];t=t[good];e1=e1[good];e2=e2[good];h=h[good];det=det[good];s=p-t[:,0];u=np.einsum('ij,ij->i',s,h)/det;k=np.cross(s,e1);v=np.einsum('j,ij->i',direction,k)/det;mu=np.einsum('ij,ij->i',e2,k)/det;wrong=np.flatnonzero((u>=0)&(v>=0)&(u+v<=1)&(mu>=0)&(mu<=1));hits.extend(dict(component=component,source_interval=interval,wall_face=int(sel[x]),fraction=float(mu[x]))for x in wrong)
bloodcheck=audit(blood,require_closed=True,expected_components=1,expected_genus=[0]);atomic_json(O/'blood-volume-audit.json',bloodcheck)
literal={tuple(sorted(map(tuple,t)))for t in outer.triangles};retention=[]
native=reader.read_part(a,rows['ZA:Ascending aorta'],frame=wall.frame);selected=np.concatenate(recover_station_rings(native.vertices,native.faces)[0]['rings'][14:]);faces=native.faces[np.all(np.isin(native.faces,selected),axis=1)]
for pid,mesh,ff in [('ZA:Ascending aorta stations 14–24',native,faces),*[(pid,m,m.faces[np.all(np.isin(m.faces,np.concatenate(rings)),axis=1)])for pid,m,rings in source_rows]]:
    kept=np.array([tuple(sorted(map(tuple,t)))in literal for t in mesh.vertices[ff]]);retention.append(dict(source=pid,source_faces=len(ff),literal_retained_faces=int(kept.sum()),modified_source_faces=np.flatnonzero(~kept).tolist()))
patency=dict(schema='male-three-sinus-root-independent-patency/v1',whole_root_native_coronary_connector_segments=len(segments),wall_crossing_hits=hits,
    minimum_blood_winding=float(winding.min()),all_midpoints_inside_one_blood_volume=bool(np.all(winding>.99)),closed_blood_volume_passed=bloodcheck['passed'],
    native_face_retention=retention,patency_pass=bool(not hits and np.all(winding>.99)and bloodcheck['passed']),basal_heart_join_verified=False,
    valve_leaflet_tissue_not_subtracted_from_classifier=True,whole_cardiac_outflow_not_claimed=True)
atomic_json(O/'patency.json',patency);atomic_npz(O/'actual-guide-segments.npz',segments=np.array([[p,q]for _,_,p,q in segments]),component=np.array([i for i,_,_,_ in segments]),source_interval=np.array([j for _,j,_,_ in segments]))

physical={r['part']for r in json.loads((ROOT/'.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/preview-proof.json').read_text())['physical_tissues']};own={'ZA:Ascending aorta','ZA:Left coronary artery','ZA:Right coronary artery'};context=[]
for part in a['parts']:
    pid=part['id'];system=part.get('system');primary=pid.startswith('ZA:');is_lung=pid.startswith('TA98PULMONARY:male:')and system=='respiratory';group='heart'if pid in physical else'lungs'if is_lung else'skeleton-muscles'if primary and system in ['skeletal','muscular']else'vessels'if primary and system in ['arterial','venous']and pid not in own else'viscera'if primary and system in ['digestive','respiratory']else None
    if not group or part.get('suppressed'):continue
    bounds=np.array(part['bounds'])
    if np.any(bounds[1]<wall.bounds[0])or np.any(bounds[0]>wall.bounds[1]):continue
    mesh=reader.read_part(a,part,frame=wall.frame)
    if group=='heart' and part.get('tissueDrawGroups'):
        exposed=[np.arange(g['start']//3,(g['start']+g['count'])//3)for g in part['tissueDrawGroups']if g.get('neighbor')not in physical]
        if not exposed:continue
        take=np.concatenate(exposed);mesh=Mesh(mesh.vertices,mesh.faces[take],mesh.frame)
    contacts,summary=pair_contacts(wall,mesh,0.,max_pairs=100000,max_tests=1000000);assert summary['complete_within_pair_budget'];wrong=np.isin(contacts['relation'],['crossing','coplanar-overlap']);pairs=contacts['pairs'][wrong];record=dict(part=pid,group=group,geometry_sha256=mesh.geometry_sha,summary=summary,witness_file=str(O/(hashlib.sha256(pid.encode()).hexdigest()[:16]+'.npz')))
    atomic_npz(record['witness_file'],candidate_faces=pairs[:,0],context_faces=pairs[:,1],frame=np.array(wall.frame));context.append(record)
    atomic_json(O/f'context-progress-{len(context):03d}.json',dict(completed_parts=len(context),elapsed_seconds=time.monotonic()-T))
report=dict(schema='male-three-sinus-root-independent-assembly/v1',patency=patency,context=context,
    crossing_parts=[dict(part=r['part'],group=r['group'],crossing_face_pairs=r['summary']['crossing_face_pairs'])for r in context if r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs']],
    rootwall_topology_pass_does_not_override_packing=True,heart_and_all_sources_fixed=True,installable=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T)
atomic_json(O/'evidence.json',report);print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'],patency=patency,crossing_parts=report['crossing_parts'])))
