"""Localize existing versus shifted tissue penetrations; no geometry edits."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts,patch_components
from geometry_tools.audit import topology,winding_contains
from geometry_tools.open3d_backend import SurfaceIndex
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1';F='model:local-ta98:assembled-metres'

def main():
    start=time.monotonic();cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';atlas=json.loads(cat.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    c=np.load(B/'course-correction-v3/candidate.npz');new=np.load(B/'bilateral-course-proof-v2/l-course.npz');a=Mesh(c['source_vertices'],new['faces'],F);b=Mesh(new['vertices'],new['faces'],F)
    station=np.full(len(a.vertices),-1,int)
    for j,ring in enumerate(c['rings']):station[ring]=j
    used=np.unique(a.faces);rows=[]
    for pid in ['ZA:Subcutaneous trochanteric bursa.l','ZA:Vastus lateralis muscle.l:Tendon']:
        p=reader.read_part(atlas,parts[pid],frame=F);t=topology(p);samples=[];contacts=[]
        for name,m in [('before',a),('after',b)]:
            q,s=pair_contacts(m,p,0.);assert s['complete_within_pair_budget'];mask=np.isin(q['relation'],['crossing','coplanar-overlap']);pairs=q['pairs'][mask];ci=np.unique(pairs[:,0]);pi=np.unique(pairs[:,1]);groups=patch_components(p,pi)
            row=dict(stage=name,crossing_pairs=len(pairs),course_faces=ci.tolist(),peer_faces=pi.tolist(),course_station_range=[int(station[m.faces[ci]].min()),int(station[m.faces[ci]].max())],peer_patch_components=[dict(face_count=len(g),bounds=p.triangles[g].reshape(-1,3).min(0).tolist()+p.triangles[g].reshape(-1,3).max(0).tolist())for g in groups])
            if t['closed_oriented_manifold']:
                membership=winding_contains(p,m.vertices[used],distance_backend='open3d');near,_=SurfaceIndex(p).closest_points(m.vertices[used],certify=True);inside=np.asarray(membership['inside']);row['inside_control_count']=int(inside.sum());row['maximum_control_penetration_m']=float(near['distance_m'][inside].max())if inside.any()else 0.;row['inside_stations']=sorted(set(station[used[inside]].tolist()))
            samples.append(row);contacts.append(pi)
        oldids=set(contacts[0]);newids=set(contacts[1]);oldvertices=set(p.faces[list(oldids)].ravel().tolist());one_ring=set(np.flatnonzero(np.any(np.isin(p.faces,list(oldvertices)),axis=1)).tolist());outside=sorted(newids-one_ring)
        rows.append(dict(id=pid,mesh_geometry_sha256=p.geometry_sha,peer_topology=t,stages=samples,new_peer_faces=sorted(newids-oldids),new_peer_faces_outside_old_one_ring=outside,interpretation='Real existing tissue penetrations move with the course; any absence of new disconnected patches does not establish anatomical validity. Bursa and tendon are not automatically authorized vessel-entry interfaces.'))
    proof=dict(schema='native-lcfa-context-localization/v1',source_catalogue_sha256=digest_file(cat),side='l',mirror_scope='Right raw actual-contact pair results exactly match left; no independent right topology assumption needed for raw counts.',rows=rows,passed_for_clearance=False,elapsed_seconds=time.monotonic()-start)
    atomic_json(B/'context-proof-v1/localization.json',proof);print(str(B/'context-proof-v1/localization.json'))

if __name__=='__main__':main()
