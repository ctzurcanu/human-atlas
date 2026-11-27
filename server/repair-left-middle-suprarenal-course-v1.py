"""One local, source-bound branch proposal; never edits or installs an atlas."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.surface_rays import first_surface_hits
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,projected_area

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/middle-suprarenal-native-v1'
OUT=BASE/'left-course-repair-v1'
CATALOGUE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres'

def section_points(mesh,point,normal):
    v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True)
    scalar=(v-point)@normal
    result=clip_scalar_surface(v,inv[mesh.faces],scalar,0.)
    rings=boundary_loops(result['negative_faces'])
    assert len(rings)==1,'Source calibre station is not one unbranched closed profile'
    return result['positions'][rings[0]]

def prepare():
    started=time.monotonic();OUT.mkdir(parents=True,exist_ok=False)
    proof=json.loads((BASE/'left'/'proof.json').read_bytes());loc=json.loads((BASE/'left'/'localized-interface-proof.json').read_bytes())
    atlas=json.loads(CATALOGUE.read_bytes());reader=AtlasReader();parts={p['id']:p for p in atlas['parts']}
    retained={p['id']:p for p in proof['peers']if p.get('mesh_npz')}
    equivalence=[]
    for pid,p in retained.items():
        old=Mesh.load(p['mesh_npz']);current=reader.read_part(atlas,parts[pid],frame=FRAME)
        assert current.geometry_sha==old.geometry_sha and current.attributes_sha==old.attributes_sha,pid+' changed since context screen'
        equivalence.append(dict(id=pid,geometry_sha256=current.geometry_sha,attributes_sha256=current.attributes_sha))
    source=reader.read_part(atlas,parts['FJ3472'],frame=FRAME)
    assert source.geometry_sha==proof['source_geometry_sha256'] and source.attributes_sha==proof['source_attributes_sha256']
    aorta=Mesh.load(retained['ZA:Abdominal aorta']['mesh_npz']);capsule=Mesh.load(retained['TA98PART:male:Suprarenal:left:capsule']['mesh_npz'])
    oldroot=np.array(loc['root_contact']['source_centroid_m']);tip=np.array(loc['distal_extremity']['point_m'])
    # Named lateral aortic origin at the retained source origin's height;
    # choose a directional ray on the actual aortic wall, not a nearest mesh.
    y=oldroot[1];tri=aorta.triangles;points=[]
    for k in range(3):
        a,b=tri[:,k],tri[:,(k+1)%3];d0=a[:,1]-y;d1=b[:,1]-y;take=d0*d1<0
        q=a[take]+((y-a[take,1])/(b[take,1]-a[take,1]))[:,None]*(b[take]-a[take]);points.extend(q)
    cross=np.array(points);assert len(cross)>6
    center=cross.mean(0);center[1]=y;direction=np.array([1.,0.,1.]);direction/=np.linalg.norm(direction)
    root_hit,root_ray=first_surface_hits(aorta,center[None],direction[None],maximum_m=.04,time_limit=3.)
    assert root_hit['hit'][0]
    root=root_hit['points'][0];root_fid=int(root_hit['face_ids'][0]);root_tri=aorta.triangles[root_fid]
    root_normal=np.cross(root_tri[1]-root_tri[0],root_tri[2]-root_tri[0]);root_normal/=np.linalg.norm(root_normal)
    if root_normal@direction<0:root_normal=-root_normal
    # Mid-height medial gland entry: an explicit medial-to-lateral ray hits
    # the first actual finite capsule boundary before the underlying cortex.
    lo,hi=capsule.vertices.min(0),capsule.vertices.max(0)
    medial=np.array([lo[0]-.03,(lo[1]+hi[1])/2,(lo[2]+hi[2])/2])
    entry_hit,entry_ray=first_surface_hits(capsule,medial[None],np.array([[1.,0.,0.]]),maximum_m=.12,time_limit=3.)
    assert entry_hit['hit'][0]
    entry=entry_hit['points'][0];entry_fid=int(entry_hit['face_ids'][0]);entry_tri=capsule.triangles[entry_fid]
    entry_normal=np.cross(entry_tri[1]-entry_tri[0],entry_tri[2]-entry_tri[0]);entry_normal/=np.linalg.norm(entry_normal)
    assert entry_normal[0]<0,'First medial capsule hit must face the incoming artery'
    normal=tip-oldroot;normal/=np.linalg.norm(normal);sections=[]
    for fraction in [.35,.5,.65]:
        ring=section_points(source,oldroot+fraction*(tip-oldroot),normal)
        sections.append(dict(fraction=fraction,points_m=ring.tolist(),area_m2=projected_area(ring),diameter_m=float(np.linalg.norm(ring[:,None]-ring[None,:],axis=2).max())))
    radius=float(np.sqrt(np.median([x['area_m2']for x in sections])/np.pi))
    # One explicit smoothly rounded superior/anterior corridor proposal.
    # It remains in the aorto-adrenal region; no obstacle is moved or waived.
    approach=entry+entry_normal*.010
    anterior=np.array([(root[0]+approach[0])/2,max(root[1],entry[1])+.007,max(root[2],approach[2])+.006])
    waypoints=np.array([root,root+root_normal*.006,anterior,approach,entry])
    plan=dict(schema='left-middle-suprarenal-single-course-plan/v1',source_catalogue=str(CATALOGUE),source_catalogue_sha256=digest_file(CATALOGUE),original_context_catalogue_sha256=proof['source_catalogue_sha256'],metadata_revision_equivalence=equivalence,inventory=str(BASE/'inventory.json'),source_id='FJ3472',source_geometry_sha256=source.geometry_sha,source_attributes_sha256=source.attributes_sha,source_license=proof['source_license'],source_provenance=proof['source_provenance'],frame=FRAME,fixed_skeleton=True,fixed_aorta=True,fixed_recipient=True,fixed_all_other_meshes=True,parent=dict(id='ZA:Abdominal aorta',mesh=retained['ZA:Abdominal aorta']['mesh_npz'],port_point_m=root.tolist(),face_id=root_fid,barycentric=root_hit['barycentric'][0].tolist(),outward_normal=root_normal.tolist(),choice='Left anterolateral directional ray at retained source origin height',ray=root_ray),recipient=dict(id='TA98PART:male:Suprarenal:left:capsule',mesh=retained['TA98PART:male:Suprarenal:left:capsule']['mesh_npz'],port_point_m=entry.tolist(),face_id=entry_fid,barycentric=entry_hit['barycentric'][0].tolist(),outward_normal=entry_normal.tolist(),choice='First capsule surface along medial-to-lateral ray at gland bounding mid-height and AP midpoint',ray=entry_ray),source_profiles=sections,outer_radius_m=radius,inner_radius_m=radius*.8,wall_thickness_scope='Illustrative positive 20%-radius wall, not source-derived or clinically validated histology',waypoints_m=waypoints.tolist(),curve_method='centripetal',station_spacing_m=.001,curve_chord_tolerance_m=.000025,ring_controls=16,max_numerical_seconds=40,max_attempts=1,install=False,anatomical_acceptance=False,interfaces=dict(parent='Named surface port, currently closed/unmodified parent; not a communicating junction',recipient='Named outer capsule port, no hole created; must not claim a patent microvascular gland interface'),references=[dict(url='https://pubmed.ncbi.nlm.nih.gov/11035646/',scope='Adult cadaver evidence for variable middle suprarenal origins; aortic-origin variant selected'),dict(url='https://medicine.uams.edu/neuroscience/education/medical-school-courses/human-structure-module/anatomy-tables/artery-tables/arteries-of-the-abdomen/',scope='Aortic branch supplies medial suprarenal gland'),dict(url='https://www.scielo.cl/scielo.php?lng=en&pid=S0717-95022019000301023&script=sci_arttext',scope='Fetal microdissection shows variable lateral/superolateral courses and close crural relationships; no adult dimensions inferred')],viewed_source_render=str(ROOT/'.local-models/geometry-workbench/review/5cfa15fcb7c713971f1750fded7fbfe8bd1a6929dc9d37cb59cb83ed89fb04db/view-2.png'),elapsed_seconds=time.monotonic()-started)
    atomic_json(OUT/'plan.json',plan)
    print(json.dumps(dict(root=root.tolist(),entry=entry.tolist(),outer_diameter_mm=radius*2000,waypoints=waypoints.tolist(),elapsed_seconds=plan['elapsed_seconds'])))

def rebase():
    """Authorized metadata-only rebase, retaining the rejected SHA guard."""
    import gzip,hashlib
    started=time.monotonic();plan=json.loads((OUT/'plan.json').read_bytes())
    snapshot=ROOT/'.local-models/ta98-review/male-native-pontine-source-controls-local-v1/0-atlas.json.gz'
    payload=gzip.decompress(snapshot.read_bytes());assert hashlib.sha256(payload).hexdigest()==plan['source_catalogue_sha256']
    old=json.loads(payload);current=json.loads(CATALOGUE.read_bytes())
    assert digest_file(CATALOGUE)=='752da7ae0744bce97407f1b582ece9e073f42b83c0a555dcf9d8f73e0f48e564'
    oldparts={p['id']:p for p in old['parts']};newparts={p['id']:p for p in current['parts']};assert oldparts.keys()==newparts.keys()
    changes=[]
    for pid,a in oldparts.items():
        b=newparts[pid];fields=sorted(k for k in a.keys()|b.keys()if a.get(k)!=b.get(k))
        if fields:
            assert fields==['ta98Term'] and 'pontine' in pid,pid+str(fields)
            changes.append(dict(id=pid,fields=fields,old=a.get('ta98Term'),new=b.get('ta98Term')))
    assert len(changes)==4
    metadata_changes=[]
    for key in old.keys()|current.keys():
        if key=='parts' or old.get(key)==current.get(key):continue
        assert key in {'concepts','ta98TermBindings'},key+' unexpectedly changed'
        metadata_changes.append(dict(field=key,old_canonical_sha256=hashlib.sha256(json.dumps(old.get(key),sort_keys=True).encode()).hexdigest(),new_canonical_sha256=hashlib.sha256(json.dumps(current.get(key),sort_keys=True).encode()).hexdigest(),scope='Nonphysical TA98 identity metadata; all part geometry records independently compared'))
    # Cover the old full source extent AND an intentionally generous envelope
    # of the entire planned route, including potential cubic overshoot.
    points=np.array(plan['waypoints_m']);low=points.min(0)-.015;high=points.max(0)+.015
    reader=AtlasReader();hidden=set(old['defaultHidden']);rows=[]
    for pid,a in oldparts.items():
        if not a.get('indexCount'):continue
        bounds=np.asarray(a['bounds']);relevant=pid=='FJ3472' or(pid not in hidden and not a.get('suppressed')and np.all(bounds[1]>=low)and np.all(bounds[0]<=high))
        if not relevant:continue
        before=reader.read_part(old,a,frame=FRAME);after=reader.read_part(current,newparts[pid],frame=FRAME)
        assert np.array_equal(before.vertices,after.vertices)and np.array_equal(before.faces,after.faces),pid
        assert before.geometry_sha==after.geometry_sha and before.attributes_sha==after.attributes_sha,pid
        assert before.attributes.keys()==after.attributes.keys() and all(np.array_equal(before.attributes[k],after.attributes[k])for k in before.attributes),pid
        rows.append(dict(id=pid,geometry_sha256=after.geometry_sha,attributes_sha256=after.attributes_sha,bank_sha256=after.source['bank_sha256']))
    evidence=dict(schema='metadata-rebase-literal-physical-equivalence/v1',old_catalogue_snapshot=str(snapshot),old_catalogue_sha256=hashlib.sha256(payload).hexdigest(),new_catalogue_sha256=digest_file(CATALOGUE),all_part_records=len(oldparts),part_metadata_changes=changes,other_nonphysical_metadata_changes=metadata_changes,all_other_catalogue_fields_equal=True,physical_source_records_all_equal_except_nonphysical_terms=True,fixed_context_envelope_m=[low.tolist(),high.tolist()],literal_decoded_context_count=len(rows),literal_decoded_context=rows,geometry_and_all_attributes_equal=True,initial_precondition_failure_receipt='9257377520954159bd502daec7e06adf',authorization='Root explicitly authorized metadata-only rebase and one previously unexecuted 40-second construction stage',elapsed_seconds=time.monotonic()-started)
    atomic_json(OUT/'rebase-equivalence.json',evidence)
    rebound={**plan,'source_catalogue_sha256':digest_file(CATALOGUE),'metadata_rebase_equivalence':str(OUT/'rebase-equivalence.json'),'metadata_rebase_equivalence_sha256':digest_file(OUT/'rebase-equivalence.json'),'original_plan_sha256':digest_file(OUT/'plan.json')}
    atomic_json(OUT/'bound-plan.json',rebound)
    print(json.dumps(dict(all_records=len(oldparts),decoded_context=len(rows),changed_nonphysical_terms=len(changes),elapsed_seconds=evidence['elapsed_seconds'])))

def candidate():
    from geometry_tools.algorithms.tube_sections import centripetal_curve,sweep_sections,TubeTolerances
    from geometry_tools.algorithms.tube_continuity import loft_faces,annular_rim
    from geometry_tools.tubes import audit_fields
    from geometry_tools.audit import audit,topology,winding_contains
    from geometry_tools.contact import pair_contacts
    from geometry_tools.spatial import BVH
    started=time.monotonic();plan_file=OUT/'bound-plan.json'if(OUT/'bound-plan.json').exists()else OUT/'plan.json';plan=json.loads(plan_file.read_bytes())
    if (OUT/'attempt-started.json').exists():
        assert plan_file.name=='bound-plan.json'and not(OUT/'candidate-wall.npz').exists(),'Only an unexecuted construction stage may resume after explicit metadata rebase'
        assert digest_file(OUT/'rebase-equivalence.json')==plan['metadata_rebase_equivalence_sha256']
        assert not(OUT/'construction-started.json').exists(),'One construction attempt only'
        atomic_json(OUT/'construction-started.json',dict(plan_sha256=digest_file(plan_file),maximum_wall_seconds=40,maximum_attempts=1,retains_initial_guard_failure=True))
    else:atomic_json(OUT/'attempt-started.json',dict(plan_sha256=digest_file(plan_file),maximum_wall_seconds=40,maximum_attempts=1))
    assert digest_file(CATALOGUE)==plan['source_catalogue_sha256'],'Catalogue changed after the explicit port plan'
    atlas=json.loads(CATALOGUE.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    source=reader.read_part(atlas,parts['FJ3472'],frame=FRAME)
    assert source.geometry_sha==plan['source_geometry_sha256'] and source.attributes_sha==plan['source_attributes_sha256']
    path,tangents=centripetal_curve(plan['waypoints_m'],chord_tolerance_m=plan['curve_chord_tolerance_m'],max_step_m=plan['station_spacing_m'])
    n=len(path);k=plan['ring_controls'];angles=np.linspace(0,2*np.pi,k,endpoint=False)
    straight=np.zeros((2,n,k,3));straight[:,:,:,2]=np.linspace(0,np.linalg.norm(np.diff(path,axis=0),axis=1).sum(),n)[None,:,None]
    for j,r in enumerate([plan['outer_radius_m'],plan['inner_radius_m']]):
        straight[j,:,:,0]=r*np.cos(angles);straight[j,:,:,1]=r*np.sin(angles)
    tol=TubeTolerances(absolute_m=1e-8,area_relative=.001,diameter_relative=.001,thickness_relative=.001,planarity_m=1e-9,max_station_turn_rad=np.deg2rad(12))
    fields,frames,transport=sweep_sections(straight,path,tangents,tol,dtype=np.float64)
    origin=np.mean(path,axis=0);encoded=(fields-origin).astype('<f4');actual=encoded.astype(float)+origin
    vertices=actual.reshape(2,-1,3);faces=loft_faces(n,k);count=n*k
    wall_faces=np.vstack([faces,faces[:,::-1]+count,annular_rim(np.arange(k),np.arange(k)+count),annular_rim(np.arange(count-k,count)[::-1],np.arange(count-k,count)[::-1]+count)])
    wall=Mesh(vertices.reshape(-1,3),wall_faces,FRAME,source=dict(parent_source_geometry_sha256=source.geometry_sha,plan_sha256=digest_file(plan_file),construction='Source-calibre circular paired finite wall with annular tissue-only end rims'))
    wall.save(OUT/'candidate-wall.npz');Mesh(vertices[0],faces,FRAME).save(OUT/'candidate-outer.npz');Mesh(vertices[1],faces[:,::-1],FRAME).save(OUT/'candidate-inner.npz')
    atomic_npz(OUT/'candidate-fields.npz',fields=encoded.reshape(2,-1,3),faces=faces,centres=path,tangents=tangents,frames=frames,geometryOrigin=origin,frame=np.array(FRAME))
    proof=dict(schema='left-middle-suprarenal-single-course-candidate/v1',source_catalogue_sha256=plan['source_catalogue_sha256'],plan_sha256=digest_file(plan_file),candidate_geometry_sha256=wall.geometry_sha,source_geometry_sha256=source.geometry_sha,complete=False,installed=False,anatomical_acceptance=False,attempt=1,source_profile_diameter_mm=plan['outer_radius_m']*2000,wall_thickness_m=plan['outer_radius_m']-plan['inner_radius_m'],wall_thickness_clinically_verified=False,transport=transport,field_audit=audit_fields(vertices,faces,FRAME,centres=path,minimum_lumen_radius_m=plan['inner_radius_m']*.97),wall_audit=audit(wall,require_closed=True,expected_components=1,expected_genus=[1]),peers=[])
    low=wall.vertices.min(0)-.002;high=wall.vertices.max(0)+.002;hidden=set(atlas['defaultHidden']);wall_bvh=BVH(wall)
    views=[dict(mesh=str(OUT/'candidate-wall.npz'),name='Left middle suprarenal artery: unjoined finite wall proposal',color=[.86,.10,.08,1.])]
    sideface_count=len(faces)
    for index,p in enumerate(atlas['parts']):
        if p.get('suppressed')or p['id']in hidden or not p.get('indexCount')or not p.get('bounds'):continue
        bounds=np.asarray(p['bounds'])
        if np.any(bounds[1]<low)or np.any(bounds[0]>high):continue
        peer=reader.read_part(atlas,p,frame=FRAME)
        if np.any(peer.vertices.max(0)<low)or np.any(peer.vertices.min(0)>high):continue
        hits,contacts=pair_contacts(wall,peer,0.,first_bvh=wall_bvh)
        face_ids=np.unique(hits['pairs'][:,0]);stations=[]
        for fid in face_ids:
            if fid<sideface_count:station=int(fid//(2*k))
            elif fid<2*sideface_count:station=int((fid-sideface_count)//(2*k))
            else:station=0 if fid<2*sideface_count+2*k else n-1
            stations.append(station)
        root=p['id']==plan['parent']['id'];entry=p['id']==plan['recipient']['id']
        localized=(root and all(s<=2 for s in stations))or(entry and all(s>=n-3 for s in stations))
        row=dict(id=p['id'],geometry_sha256=peer.geometry_sha,attributes_sha256=peer.attributes_sha,contact=contacts,source_wall_face_ids=face_ids.tolist(),contact_station_range=[min(stations),max(stations)]if stations else None,declared_port_only=bool(localized),contact_permitted=bool(localized),containment=None)
        # Surface contacts are never a complete solid-occupancy proof. Query
        # all guide stations inside this peer's actual AABB if it is a valid
        # closed oriented volume. Open fascia stays explicitly unclassified.
        inside_bounds=np.all(path>=peer.vertices.min(0),axis=1)&np.all(path<=peer.vertices.max(0),axis=1)
        if inside_bounds.any() and not root and not entry:
            topo=topology(peer)
            if topo['closed_oriented_manifold']:
                positions=np.flatnonzero(inside_bounds);inside=winding_contains(peer,path[positions],distance_backend='bvh')
                row['containment']=dict(path_station_ids=positions.tolist(),evidence=inside,occupied_station_ids=positions[np.array(inside['inside'])].tolist())
            else:row['containment']=dict(valid_closed_volume=False,scope='Open source sheet: only actual triangle contacts are classified')
        has_contact=contacts['returned_triangle_pairs']>0
        occupied=row['containment']and bool(row['containment'].get('occupied_station_ids'))
        row['fixed_context_blocker']=bool((has_contact and not localized)or occupied)
        proof['peers'].append(row);proof['elapsed_seconds']=time.monotonic()-started
        if has_contact or root or entry:
            path_npz=OUT/f'fixed-peer-{index:04}.npz';peer.save(path_npz)
            color=[.62,.17,.15,.55]if root else[.75,.63,.34,.45]if entry else[.42,.51,.66,.30]
            views.append(dict(mesh=str(path_npz),name=p['name'],color=color))
        atomic_json(OUT/'context-checkpoint.json',proof,overwrite=(OUT/'context-checkpoint.json').exists())
    proof['remaining_fixed_context_blockers']=[p['id']for p in proof['peers']if p['fixed_context_blocker']]
    proof['complete']=True;proof['all_context_pair_budgets_complete']=all(p['contact']['complete_within_pair_budget']for p in proof['peers'])
    proof['numerical_course_gate_pass']=bool(proof['transport']['section_transport_pass']and proof['field_audit']['passed']and proof['wall_audit']['passed']and proof['all_context_pair_budgets_complete']and not proof['remaining_fixed_context_blockers'])
    proof['communicating_parent_root_verified']=False;proof['recipient_hole_verified']=False
    proof['interface_limit']='An unmodified aortic outer carrier and intact gland capsule cannot establish communicating lumens; named port contact is not an anatomical join.'
    proof['elapsed_seconds']=time.monotonic()-started;atomic_json(OUT/'proof.json',proof)
    atomic_json(OUT/'render-spec.json',dict(frame=FRAME,wall_seconds=40,resolution=900,meshes=views,views=[dict(name='left-repaired-front-context',direction=[0,0,1]),dict(name='left-repaired-oblique-context',direction=[1,.2,.4]),dict(name='left-repaired-lateral-context',direction=[1,0,0]),dict(name='left-repaired-diagnostic-open-lumen',direction=[1,.2,.4],section=dict(normal=[0,0,1],point=path[len(path)//2].tolist(),clear_outer=True))]))
    print(json.dumps(dict(numerical_course_gate_pass=proof['numerical_course_gate_pass'],transport_pass=proof['transport']['section_transport_pass'],field_pass=proof['field_audit']['passed'],wall_pass=proof['wall_audit']['passed'],fixed_blockers=proof['remaining_fixed_context_blockers'],checked_peers=len(proof['peers']),stations=n,elapsed_seconds=proof['elapsed_seconds'])))

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='rebase':rebase()
    elif sys.argv[1]=='candidate':candidate()
