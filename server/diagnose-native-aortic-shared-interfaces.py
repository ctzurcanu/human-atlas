"""Diagnostic-only zero-area filtering and exact native ownership sharing."""
from pathlib import Path
import sys, json, time
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, atomic_npz
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from collections import defaultdict

started=time.monotonic(); output=ROOT/'.local-models/geometry-workbench/heart/native-aortic-source-interface-v1';output.mkdir(exist_ok=False)
reader=AtlasReader();records=[]
for sex in ['male','female']:
    atlas=json.loads((ROOT/f'.local-models/{sex}.json').read_text())
    sources={p['sourceId']:p for p in atlas['parts'] if p.get('system')=='cardiac' and p.get('sourceId') in ['aortic_semilunar_valve','aortic_arch','ventricles']}
    clean={}
    for name,part in sources.items():
        original=reader.read_part(atlas,part,frame=f'model:{sex}-detail:assembled-metres')
        used=np.unique(original.faces);vertices,inv=np.unique(original.vertices[used],axis=0,return_inverse=True);remap=np.full(len(original.vertices),-1,int);remap[used]=inv;faces=remap[original.faces]
        triangles=vertices[faces];area=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)
        face_ids=np.flatnonzero(area>0); mesh=Mesh(vertices,faces[face_ids],original.frame)
        clean[name]=mesh
        atomic_npz(output/f'{sex}-{name}-diagnostic.npz',vertices=vertices,faces=mesh.faces,source_face_ids=face_ids,frame=np.array(mesh.frame))
        records.append(dict(sex=sex,part=part['id'],original_geometry_sha256=original.geometry_sha,original_triangles=len(original.faces),zero_area_omitted_for_diagnostic_only=int((area==0).sum()),topology=topology(mesh)))
    valve=clean['aortic_semilunar_valve'];matches=[]
    vmap={tuple(v):i for i,v in enumerate(valve.vertices)}
    for host_name in ['ventricles','aortic_arch']:
        host=clean[host_name];hmap={tuple(v):i for i,v in enumerate(host.vertices)}
        shared=sorted(set(vmap)&set(hmap));edges=[];vf=set()
        for face in valve.faces:
            for aa,bb in zip(face,np.roll(face,-1)):
                if tuple(valve.vertices[aa]) in hmap and tuple(valve.vertices[bb]) in hmap:
                    edges.append(tuple(sorted((int(aa),int(bb)))))
        hfaces={tuple(sorted(tuple(host.vertices[i]) for i in face)) for face in host.faces}
        common=[fi for fi,face in enumerate(valve.faces)if tuple(sorted(tuple(valve.vertices[i]) for i in face)) in hfaces]
        d,fs,qs=igl.point_mesh_squared_distance(valve.vertices,host.vertices,host.faces)
        distance=np.sqrt(d)*1000
        record=dict(host=host_name,exact_shared_vertex_count=len(shared),exact_shared_valve_edge_count=len(set(edges)),exact_shared_valve_triangle_count=len(common),
                    valve_vertex_distance_mm=dict(minimum=float(distance.min()),median=float(np.median(distance)),maximum=float(distance.max()),within_0_01mm=int((distance<=.01).sum()),within_0_25mm=int((distance<=.25).sum())))
        matches.append(record)
        atomic_npz(output/f'{sex}-valve-{host_name}-witnesses.npz',valve_vertices=valve.vertices,valve_faces=valve.faces,shared_valve_edges=np.asarray(sorted(set(edges)),int).reshape(-1,2),valve_host_distance_mm=distance,host_witness_faces=fs,host_witness_points=qs)
    records.append(dict(sex=sex,actual_native_valve_host_interfaces=matches))
report=dict(schema='native-aortic-source-exact-interface/v1',source_geometry_modified=False,annulus_created=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started,records=records)
atomic_json(output/'evidence.json',report);print(json.dumps(report))
