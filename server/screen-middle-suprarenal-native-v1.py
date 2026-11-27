"""Read-only source quality, end-incidence and fixed-context screening."""
from pathlib import Path
import hashlib, json, sqlite3, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.ta98_authority import corrected_terms

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/middle-suprarenal-native-v1'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='043e699fff254b38460d0f008046a11f2291365f57cd84bd2246a8de4c574c54'

def patch(mesh, faceids):
    faceids=np.unique(faceids)
    if not len(faceids):return dict(face_count=0)
    tris=mesh.faces[faceids];v=mesh.vertices[np.unique(tris)];owners={};neighbors={int(i):set()for i in faceids}
    for fid,tri in zip(faceids,tris):
        for x,y in zip(tri,np.roll(tri,-1)):
            key=tuple(sorted((int(x),int(y))))
            for old in owners.get(key,[]):neighbors[int(fid)].add(old);neighbors[old].add(int(fid))
            owners.setdefault(key,[]).append(int(fid))
    unvisited=set(neighbors);components=[]
    while unvisited:
        todo=[unvisited.pop()];group=[]
        while todo:
            f=todo.pop();group.append(f)
            for other in neighbors[f]&unvisited:unvisited.remove(other);todo.append(other)
        points=mesh.vertices[np.unique(mesh.faces[group])];components.append(dict(faces=len(group),bounds=[points.min(0).tolist(),points.max(0).tolist()]))
    return dict(face_count=len(faceids),bounds=[v.min(0).tolist(),v.max(0).tolist()],components=components)

def prepare():
    assert digest_file(SOURCE)==EXPECTED;B.mkdir(exist_ok=True,parents=True)
    assert not (B/'inventory.json').exists(), 'Preserve existing inventory'
    a=json.loads(SOURCE.read_bytes());authority=ROOT/'.local-models/ta98-authority.sqlite';db=sqlite3.connect(authority);db.row_factory=sqlite3.Row
    table='ta98'
    terms=corrected_terms({r['id']:dict(r)for r in db.execute('select * from '+table)})
    selected={'A12.2.12.074'}
    while True:
        more={k for k,t in terms.items()if t.get('parent_id')in selected}
        if more<=selected:break
        selected|=more
    inventory=dict(schema='ta98-geometry-inventory/v1',root_term='A12.2.12.074',term_count=len(selected),applicable_term_count=len(selected),model='local-ta98',model_sex='male',source_catalogue=str(SOURCE),source_catalogue_sha256=EXPECTED,authority_path=str(authority),authority_sha256=digest_file(authority),authority_corrections_path=str(ROOT/'geometry_tools/ta98_authority.py'),authority_corrections_sha256=digest_file(ROOT/'geometry_tools/ta98_authority.py'),checked_in_report=str(ROOT/'reports/mesh-audit/REVIEW.md'),terms=[dict(authority=terms[k],representation_class='whole-organ',source_meshes=[dict(part_id=p['id'],source_catalogue=str(SOURCE))for p in a['parts']if p.get('ta98Term')==k],representation_provisional=True,distinct_representation_verified=False,missing_distinct_representation=True,planned_representation=dict(scope='Inspect hidden closed outside-surface casts against parent aorta and recipient suprarenal gland. No finite wall or lumen is established.',part_ids=['FJ3472','FJ3586']))for k in sorted(selected)])
    atomic_json(B/'inventory.json',inventory)
    plan=dict(source_catalogue_sha256=EXPECTED,source_parts=['FJ3472','FJ3586'],fixed_skeleton=True,fixed_all_other_geometry=True,geometry_changes_planned=False,allow_source_preview_only_after_render_and_context=True,declared_interfaces=[dict(part='FJ3472',parent='ZA:Abdominal aorta',recipient='TA98PART:male:Suprarenal:left:capsule'),dict(part='FJ3586',parent='ZA:Abdominal aorta',recipient='TA98PART:male:Suprarenal:right:capsule')],interface_scope='Aorta root crossings and recipient capsule entry crossings are evidence to localize, not valid finite junctions. Other tissue crossings block preview acceptance.',tissue_classes=dict(vessel_source='closed outer cast, not finite wall',aorta='retained open outside vessel source',recipient='finite capsule tissue with cortex/medulla inside'),local_only=True,anatomical_acceptance=False)
    atomic_json(B/'assembly-plan.json',plan)
    report=ROOT/'reports/mesh-audit/REVIEW.md';marker='### Middle suprarenal source assembly inventory — 2026-10-26'
    assert marker not in report.read_text()
    with report.open('a')as out:
        out.write('\n\n'+marker+'\n\nBefore possible local preview activation, the complete authority closure is inventoried. Native current skeleton and all neighboring tissues stay fixed. BP3D hidden sources FJ3472/FJ3586 are closed outside casts; no communicating wall or lumen is claimed. Root A12.2.12.074 has '+str(len(selected))+' inclusive term(s). Actual root and adrenal-entry crossing patches, whole source quality and all active neighbors must be inspected before publication.\n\n| TA98 | Parent | Term | Source geometry | Representation |\n|---|---|---|---|---|\n')
        for code in sorted(selected):t=terms[code];out.write('| '+code+' | '+str(t.get('parent_id'))+' | '+t['name_en']+' | FJ3472, FJ3586 (hidden) | Outside cast; no finite wall or verified interface |\n')
    print(json.dumps(dict(inventory_terms=len(selected),out=str(B))))

def screen(side):
    start=time.monotonic();assert digest_file(SOURCE)==EXPECTED;plan=json.loads((B/'assembly-plan.json').read_bytes());atlas=json.loads(SOURCE.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';sid='FJ3472'if side=='left'else'FJ3586';src=reader.read_part(atlas,parts[sid],frame=frame);quality=audit(src,require_closed=True,expected_components=1,expected_genus=[0]);out=B/side;out.mkdir(exist_ok=False);src.save(out/'source.npz');hidden=set(atlas['defaultHidden']);low=src.vertices.min(0)-.002;high=src.vertices.max(0)+.002
    target=next(x for x in plan['declared_interfaces']if x['part']==sid);interfaces={target['parent'],target['recipient']};candidates=[]
    for p in atlas['parts']:
        if p['id']==sid or p.get('suppressed')or p['id']in hidden or not p.get('indexCount')or not p.get('bounds'):continue
        b=np.array(p['bounds']);
        if np.all(b[1]>=low)and np.all(b[0]<=high):candidates.append(p)
    proof=dict(schema='middle-suprarenal-whole-fixed-context/v1',side=side,source_id=sid,source_catalogue_sha256=EXPECTED,source_geometry_sha256=src.geometry_sha,source_attributes_sha256=src.attributes_sha,source_bank_sha256=src.source['bank_sha256'],source_provenance=parts[sid].get('provenance'),source_license=parts[sid].get('license'),source_quality=quality,broadphase_count=len(candidates),peers=[],complete=False,anatomical_acceptance=False,installed=False)
    views=[dict(mesh=str(out/'source.npz'),name=parts[sid]['name'],color=[.85,.12,.10,1.])]
    for index,p in enumerate(candidates):
        mesh=reader.read_part(atlas,p,frame=frame);v=mesh.vertices
        if np.any(v.max(0)<low)or np.any(v.min(0)>high):continue
        hit,contact=pair_contacts(src,mesh,0.);bad=np.isin(hit['relation'],['crossing','coplanar-overlap','touching']);entry=dict(id=p['id'],system=p['system'],contact=contact,declared_interface=p['id']in interfaces)
        if bad.any():entry['source_patch']=patch(src,hit['pairs'][bad,0]);entry['peer_patch']=patch(mesh,hit['pairs'][bad,1]);npz=out/f'peer-{index:03}.npz';mesh.save(npz);entry['mesh_npz']=str(npz)
        proof['peers'].append(entry);proof['elapsed_seconds']=time.monotonic()-start;atomic_json(out/'checkpoint.json',proof,overwrite=(out/'checkpoint.json').exists())
        if p['id']in interfaces:
            npz=out/f'peer-{index:03}.npz'
            if not npz.exists():mesh.save(npz)
            views.append(dict(mesh=str(npz),name=p['name'],color=[.60,.20,.18,1.]if p['id']==target['parent']else[.73,.60,.32,1.]))
    risks=[p['id']for p in proof['peers']if not p['declared_interface']and(p['contact']['crossing_face_pairs']or p['contact']['coplanar_overlap_face_pairs'])];proof.update(complete=True,unapproved_crossing_peers=risks,numerical_preview_gate_pass=quality['passed']and not risks and all(p['contact']['complete_within_pair_budget']for p in proof['peers']),elapsed_seconds=time.monotonic()-start);atomic_json(out/'proof.json',proof)
    center=(src.vertices.min(0)+src.vertices.max(0))/2;atomic_json(out/'render-spec.json',dict(frame=frame,wall_seconds=40,resolution=900,meshes=views,views=[dict(name=side+'-front',direction=[0,0,1]),dict(name=side+'-posterior',direction=[0,0,-1]),dict(name=side+'-oblique',direction=[1 if side=='left'else-1,.2,.4]),dict(name=side+'-uncapped-section',direction=[0,0,1],section=dict(normal=[0,0,1],point=center.tolist(),clear_outer=True))]))
    print(json.dumps(dict(side=side,quality_pass=quality['passed'],peers=len(proof['peers']),risks=risks,elapsed_seconds=proof['elapsed_seconds'],out=str(out))))

def localize(side):
    """Retain all raw contacts; distinguish recipient tissue from fixed obstacles.

    Surface geodesic extrema are diagnostic positions, never inferred mouths.
    No cast geometry, visibility or installed atlas is changed.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    from geometry_tools.model import Mesh
    from geometry_tools.audit import winding_contains
    started=time.monotonic();assert digest_file(SOURCE)==EXPECTED
    out=B/side;proof=json.loads((out/'proof.json').read_bytes());src=Mesh.load(out/'source.npz')
    assert src.geometry_sha==proof['source_geometry_sha256']
    recipient_prefix=f'TA98PART:male:Suprarenal:{side}:'
    declared={'ZA:Abdominal aorta'}|{recipient_prefix+x for x in ['capsule','cortex','medulla']}
    peers={p['id']:p for p in proof['peers']}
    vertices,inverse=np.unique(src.vertices,axis=0,return_inverse=True)
    faces=inverse[src.faces];edges=np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    edges=np.unique(np.sort(edges,axis=1),axis=0)
    lengths=np.linalg.norm(vertices[edges[:,0]]-vertices[edges[:,1]],axis=1)
    graph=coo_matrix((np.tile(lengths,2),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(vertices),len(vertices))).tocsr()
    root_mesh=Mesh.load(peers['ZA:Abdominal aorta']['mesh_npz'])
    hit,root_evidence=pair_contacts(src,root_mesh,0.)
    root_faces=np.unique(hit['pairs'][:,0]);root_vertices=np.unique(faces[root_faces])
    assert len(root_vertices) and root_evidence['complete_within_pair_budget']
    distance=dijkstra(graph,directed=False,indices=root_vertices,min_only=True)
    assert np.isfinite(distance).all()
    farthest=int(np.argmax(distance));tip_cluster=distance>=distance.max()-.002
    tip_point=vertices[farthest];distal_center=vertices[tip_cluster].mean(0)
    root_center=vertices[root_vertices].mean(0)
    root_patch=vertices[root_vertices]
    evidence=dict(schema='middle-suprarenal-root-recipient-localization/v1',side=side,source_id=proof['source_id'],source_catalogue_sha256=EXPECTED,source_geometry_sha256=src.geometry_sha,source_attributes_sha256=src.attributes_sha,source_geometry_unchanged=True,geometry_changes=False,installed=False,anatomical_acceptance=False,source_has_verified_open_mouth=False,diagnostic_basis='Exact-welded source triangle-edge geodesic distance from actual aortic contact faces. A farthest surface extremity is not an anatomical mouth or centerline terminal.',root_contact=dict(contact=root_evidence,source_face_ids=root_faces.tolist(),source_centroid_m=root_center.tolist(),bounds_m=[root_patch.min(0).tolist(),root_patch.max(0).tolist()]),distal_extremity=dict(source_vertex_exact_welded_index=farthest,point_m=tip_point.tolist(),geodesic_distance_from_root_patch_m=float(distance.max()),terminal_two_mm_source_vertex_count=int(tip_cluster.sum()),terminal_two_mm_centroid_m=distal_center.tolist(),distance_from_root_centroid_m=float(np.linalg.norm(tip_point-root_center))),declared_recipient_cells=sorted(declared-{'ZA:Abdominal aorta'}),recipient_occupancy={},localized_contacts=[],strict_unapproved_crossing_peers=proof['unapproved_crossing_peers'])
    for pid,p in peers.items():
        if not p.get('mesh_npz'):continue
        mesh=Mesh.load(p['mesh_npz']);hits,summary=pair_contacts(src,mesh,0.)
        source_faces=np.unique(hits['pairs'][:,0]);indices=np.unique(faces[source_faces])
        row=dict(id=pid,contact=summary,source_face_ids=source_faces.tolist(),classification='declared-aortic-root'if pid=='ZA:Abdominal aorta'else'declared-adrenal-recipient-tissue'if pid in declared else'fixed-obstacle',geodesic_range_from_root_patch_m=[float(distance[indices].min()),float(distance[indices].max())]if len(indices)else None,source_patch_centroid_m=vertices[indices].mean(0).tolist()if len(indices)else None)
        evidence['localized_contacts'].append(row)
        if pid.startswith(recipient_prefix):
            points=np.array([root_center,tip_point,distal_center])
            evidence['recipient_occupancy'][pid]=dict(point_order=['aortic_contact_centroid','distal_surface_extremity','distal_two_mm_cluster_centroid'],result=winding_contains(mesh,points,distance_backend='bvh'))
    evidence['fixed_obstacle_crossing_peers']=[p['id']for p in evidence['localized_contacts']if p['classification']=='fixed-obstacle'and(p['contact']['crossing_face_pairs']or p['contact']['coplanar_overlap_face_pairs'])]
    evidence['anatomically_declared_recipient_context_gate_pass']=not evidence['fixed_obstacle_crossing_peers'] and all(p['contact']['complete_within_pair_budget']for p in evidence['localized_contacts'])
    evidence['elapsed_seconds']=time.monotonic()-started
    atomic_json(out/'localized-interface-proof.json',evidence)
    print(json.dumps(dict(side=side,root_centroid=root_center.tolist(),tip=tip_point.tolist(),root_to_tip_geodesic_mm=float(distance.max()*1000),recipient_occupancy={k:v['result']['inside']for k,v in evidence['recipient_occupancy'].items()},remaining_fixed_obstacles=evidence['fixed_obstacle_crossing_peers'],elapsed_seconds=evidence['elapsed_seconds'])))

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='localize':localize(sys.argv[2])
    else:screen(sys.argv[1])
