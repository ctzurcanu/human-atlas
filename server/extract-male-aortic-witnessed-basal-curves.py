"""Actual native sections and cusp/host intersection arcs for hinge review.

All curve points retain both anatomical source-face witnesses. Display tubes
are diagnostic only, never tissue, an accepted annulus, or an installation.
"""
from pathlib import Path
import sys,json,time,importlib.util
from collections import defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.distance import exact_intersection_segment,barycentric

started=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-witnessed-basal-curves-v1';D.mkdir(exist_ok=False)
cp=ROOT/'.local-models/geometry-workbench/jobs/b9c8595514cdc1b8671263ec425dbf7bbad8f059c484638ae4b24211c0b2d9dc/prepared-2/proposal.json'
atlas=json.loads(cp.read_text());rows={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
cusp_ids=['ZA:Right coronary leaflet','ZA:Left coronary leaflet','ZA:Non-coronary leaflet']
cusps={pid:reader.read_part(atlas,rows[pid],frame=frame)for pid in cusp_ids}
hosts={pid:reader.read_part(atlas,rows[pid],frame=frame)for pid in ['TA98HEART:male:endocardium','TA98HEART:male:myocardium:left-ventricle']}

def display_segments(points,segments,path,radius=.00006):
    vertices=[];faces=[]
    for aa,bb in segments:
        start,end=points[aa],points[bb];axis=end-start;length=np.linalg.norm(axis)
        if length<1e-12:continue
        axis/=length;u=np.cross(axis,np.eye(3)[np.argmin(abs(axis))]);u/=np.linalg.norm(u);v=np.cross(axis,u)
        ring=np.cos(np.arange(6)*np.pi/3)[:,None]*u+np.sin(np.arange(6)*np.pi/3)[:,None]*v;base=len(vertices)
        vertices.extend(start+radius*ring);vertices.extend(end+radius*ring)
        for k in range(6):
            j=(k+1)%6;faces.extend([[base+k,base+j,base+6+j],[base+k,base+6+j,base+6+k]])
    Mesh(np.asarray(vertices),np.asarray(faces,int).reshape(-1,3),frame).save(path)

records=[];entries=[];colors=[[.15,.76,.35,1],[.2,.48,.93,1],[.94,.68,.2,1]]
for j,(pid,cusp)in enumerate(cusps.items()):
    path=D/f'cusp-{j}.npz';cusp.save(path);entries.append(dict(mesh=str(path),name=pid+' unchanged actual cusp',color=colors[j]))
    for k,(host_id,host)in enumerate(hosts.items()):
        contacts,summary=pair_contacts(cusp,host,0.,max_pairs=100000,max_tests=1000000)
        assert summary['complete_within_pair_budget']
        pairs=contacts['pairs'][contacts['relation']=='crossing'];endpoints=[];face_pairs=[];a_bary=[];b_bary=[];errors=[]
        for fa,fb in pairs:
            hit=exact_intersection_segment(cusp.triangles[fa],host.triangles[fb]);points=hit['endpoints']
            if len(points)!=2 or np.linalg.norm(points[1]-points[0])<1e-14:continue
            endpoints.extend(points);face_pairs.append([fa,fb]);errors.extend(hit['coordinate_error_upper_m'])
            a_bary.extend(barycentric(points,np.broadcast_to(cusp.triangles[fa],(2,3,3)))[0]);b_bary.extend(barycentric(points,np.broadcast_to(host.triangles[fb],(2,3,3)))[0])
        endpoints=np.asarray(endpoints);points,inverse=np.unique(endpoints,axis=0,return_inverse=True);segments=inverse.reshape(-1,2)
        graph=defaultdict(set)
        for aa,bb in segments:graph[int(aa)].add(int(bb));graph[int(bb)].add(int(aa))
        pending=set(graph);components=[]
        while pending:
            first=min(pending);todo=[first];ids=set();pending.remove(first)
            while todo:
                at=todo.pop();ids.add(at)
                for nxt in graph[at]:
                    if nxt in pending:pending.remove(nxt);todo.append(nxt)
            xyz=points[sorted(ids)];components.append(dict(controls=len(ids),degrees=sorted(set(len(graph[i])for i in ids)),
                center_m=xyz.mean(0).tolist(),extent_mm=((xyz.max(0)-xyz.min(0))*1000).tolist()))
        tag=f'cusp-{j}-host-{k}';npz=D/(tag+'.npz');atomic_npz(npz,points=points,segments=segments,
            endpoint_points=endpoints,source_face_pairs=np.asarray(face_pairs,int),cusp_barycentric=np.asarray(a_bary),
            host_barycentric=np.asarray(b_bary),coordinate_error_upper_m=np.asarray(errors),frame=np.array(frame))
        display=D/(tag+'-diagnostic-display.npz');display_segments(points,segments,display)
        # LV arcs are shown separately from endocardial arcs to avoid mistaking
        # overlapping finite-layer boundaries for two distinct anatomical hinges.
        records.append(dict(cusp=pid,host=host_id,source_cusp_geometry_sha256=cusp.geometry_sha,source_host_geometry_sha256=host.geometry_sha,
             exact_non_coplanar_segments=len(segments),components=components,witness_file=str(npz),
             anatomical_hinge_classification='Unclassified actual tissue intersection; requires basal/free-edge/coaptation render review',
             display=str(display),maximum_exact_to_float_endpoint_error_m=float(np.max(errors)) if len(errors) else 0.))
for k,(host_id,host)in enumerate(hosts.items()):
    path=D/f'host-{k}.npz';host.save(path)
    spec=dict(frame=frame,resolution=1100,wall_seconds=50,meshes=[*entries,dict(mesh=str(path),name=host_id+' unchanged actual context',color=[.64,.35,.34,.12]),
        *[dict(mesh=r['display'],name='Actual source intersection — hinge unclassified',color=[.95,.1,.62,1])for r in records if r['host']==host_id]],
        views=[dict(name='Unchanged cusps and exact basal host intersections',direction=[.4,.7,1]),dict(name='Exact tissue-intersection side profile',direction=[-1,.1,.4])])
    atomic_json(D/f'installed-host-{k}-review-spec.json',spec)

# The native toroidal valve is a local unverified-rights orientation reference,
# kept in its original frame and never copied or aligned into the current heart.
source=ROOT/'.local-models/male.json';native=json.loads(source.read_text());nparts={p.get('sourceId'):p for p in native['parts']if p.get('system')=='cardiac'}
native_meshes={name:reader.read_part(native,nparts[name],frame='model:male-detail:assembled-metres')for name in ['aortic_semilunar_valve','aortic_arch','ventricles']}
native_entries=[]
for name,mesh in native_meshes.items():
    path=D/(name+'-native-reference.npz');mesh.save(path);native_entries.append(dict(mesh=str(path),name=name+' native local reference only',color=[.12,.75,.45,1]if name=='aortic_semilunar_valve'else[.65,.4,.34,.16]))
atomic_json(D/'native-reference-review-spec.json',dict(frame='model:male-detail:assembled-metres',resolution=1100,wall_seconds=50,meshes=native_entries,
    views=[dict(name='Native annular valve in its original source context',direction=[.4,.7,1]),dict(name='Native root source side profile',direction=[-1,.1,.4])]))

# Literal triangle/plane contours retain source face/corner witnesses. A paired
# contour is not automatically an anatomical hinge or a host connection.
spec=importlib.util.spec_from_file_location('literal_section',ROOT/'server/extract-female-tongue-course-guides.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
valve=native_meshes['aortic_semilunar_valve'];center=valve.vertices.mean(0);_,sv,basis=np.linalg.svd(valve.vertices-center,full_matrices=False);axis=basis[-1]
if axis[1]<0:axis=-axis
u=basis[0];u-=axis*(u@axis);u/=np.linalg.norm(u);v=np.cross(axis,u);plane_basis=np.array([u,v]);lo,hi=(valve.vertices@axis).min(),(valve.vertices@axis).max();sections=[]
for i,level in enumerate(np.linspace(lo+.0003,hi-.0003,17)):
    hit=module.section(valve,axis,float(level),plane_basis);curves=[]
    for k,loop in enumerate(hit['loops']):
        p=loop['positions'];source_faces=loop['source_face'];bary=loop['barycentric'];residual=float(np.linalg.norm(np.einsum('ni,nij->nj',bary,valve.triangles[source_faces])-p,axis=1).max())
        path=D/f'native-section-{i:02d}-{k:02d}.npz';atomic_npz(path,points=p,source_face=source_faces,source_barycentric=bary,frame=np.array(valve.frame));curves.append(dict(witness_file=str(path),controls=len(p),area_mm2=float(loop['area_m2']*1e6),maximum_barycentric_residual_m=residual))
    sections.append(dict(plane_offset_m=hit['plane_offset_m'],contours=curves,failure=hit['failure']))
report=dict(schema='male-aortic-source-witnessed-basal-curves/v1',source_catalogue_sha256=digest_file(cp),source_geometry_changed=False,
             source_reference_scope='Local Zygote mirror; unverified redistribution rights, reference review only; no precise copied replacement mesh',
             annulus_certified=False,actual_installed_intersection_curves=records,native_section_axis=axis.tolist(),native_sections=sections,
             elapsed_seconds=time.monotonic()-started,anatomical_acceptance=False)
atomic_json(D/'evidence.json',report);print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'],curves=[dict(cusp=r['cusp'],host=r['host'],segments=r['exact_non_coplanar_segments'],components=r['components'])for r in records],native_section_contour_counts=[len(s['contours'])for s in sections])))
