"""One new anatomical facial course plan, with fixed native occupancy gates.

The supplied description tables and primary dissection study establish route
relationships.  They supply no copied private course or unregistered vertices.
"""
from pathlib import Path
import importlib.util
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.audit import topology
from geometry_tools.algorithms.tube_sections import centripetal_curve
from geometry_tools.algorithms.tube_station_topology import recover_station_rings, align_station_correspondence

EXPECTED = '5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'
BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion'
D = BASE / 'male-left-facial-cervical-routes-anatomical-v5'
source = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'


def main():
    began = time.monotonic()
    assert digest_file(source) == EXPECTED and not D.exists()
    D.mkdir()
    old = json.loads((BASE / 'male-facial-family-inputs-v1/inputs.json').read_bytes())
    atlas = json.loads(source.read_bytes())
    parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    roles = {r['role']: r['part_id'] for r in old['fixed_sources'] if not r.get('missing_native_source')}
    roles.update({'styloglossus-left':'TA98TONGUE:male:styloglossus:left', 'palatopharyngeus-left':'ZA:Palatopharyngeus muscle.l',
                  'medial-pterygoid-left':'ZA:Medial pterygoid muscle.l', 'mandible':'ZA:Mandible:Bone',
                  'stylohyoid-left':'ZA:Stylohyoid muscle.l', 'ascending-pharyngeal-left':'ZA:Ascending pharyngeal artery.l',
                  'pharyngeal-mucosa':'ZA:Oropharynx'})
    meshes, records = {}, {}
    for role, pid in roles.items():
        mesh = reader.read_part(atlas, parts[pid], frame=old['frame'], catalogue_path=source, catalogue_sha=EXPECTED)
        path = D/(role+'.npz')
        mesh.save(path)
        meshes[role] = mesh
        records[role] = dict(role=role,part_id=pid,source_record=parts[pid],mesh=str(path.relative_to(ROOT)),mesh_sha256=digest_file(path),geometry_sha256=mesh.geometry_sha)
    parent = meshes['facial-left']
    component = recover_station_rings(parent.vertices,parent.faces)[0]
    rings,replay = align_station_correspondence(parent.vertices,parent.faces[component['face_ids']],component['rings'])
    centers = parent.vertices[np.array(rings)].mean(1)
    records['facial-left']['source_station_components']=[dict(rings=[r.tolist()for r in rings],centers_m=centers.tolist(),source_faces=component['face_ids'].tolist(),phase_replay=replay)]
    spec=importlib.util.spec_from_file_location('facial_corridor_helpers',ROOT/'server/plan-facial-cervical-corridors.py')
    helpers=importlib.util.module_from_spec(spec);spec.loader.exec_module(helpers)
    blocked_roles=['mandible','facial-vein-left','submandibular-duct-left','medial-pterygoid-left','styloglossus-left','stylopharyngeus-left','ascending-pharyngeal-left','submental-left']
    classifiers={role:helpers.closed_classifier(meshes[role])for role in blocked_roles+['superior-constrictor-left']}
    target_reports=[]

    def recipient(role,seed,radius,exceptions=()):
        m=meshes[role];assert topology(m)['closed_oriented_manifold']
        lo,hi=m.bounds;seed=np.asarray(seed)
        # Keep the complete source terminal calibre.  Choose an adjacent real
        # tissue volume that avoids fixed bone/veins/ducts, not a name center.
        points=np.stack(np.meshgrid(*[np.linspace(max(lo[i],seed[i]-.006),min(hi[i],seed[i]+.006),17)for i in range(3)],indexing='ij'),-1).reshape(-1,3)
        inside=abs(igl.winding_number(m.vertices,m.faces,points))>.99
        distance=np.sqrt(igl.point_mesh_squared_distance(points,m.vertices,m.faces)[0])
        safe=inside&(distance>radius+.00015)
        rejection_counts={}
        for blocker in blocked_roles+(['superior-constrictor-left']if role=='palatopharyngeus-left'else[]):
            if blocker in exceptions:continue
            occupied=classifiers[blocker];mesh=meshes[blocker]
            ids=np.flatnonzero(safe)
            if not len(ids):break
            d=np.sqrt(igl.point_mesh_squared_distance(points[ids],mesh.vertices,mesh.faces)[0]);bad=d<=radius+.00015
            if occupied is not None:bad|=abs(igl.winding_number(occupied.vertices,occupied.faces,points[ids]))>.99
            rejection_counts[blocker]=int(bad.sum());safe[ids[bad]]=False
        valid=np.flatnonzero(safe)
        report=dict(recipient=role,radius_m=radius,search_seed_m=seed.tolist(),candidate_points=len(points),valid_complete_profile_balls=len(valid),fixed_occupancy_rejections=rejection_counts,
                    source_tissue_geometry_sha256=m.geometry_sha)
        if not len(valid):target_reports.append(report);return None,None
        k=valid[np.argmin(np.linalg.norm(points[valid]-seed,axis=1))]
        witness=dict(recipient_role=role,point_m=points[k].tolist(),minimum_surface_distance_m=float(distance[k]),terminal_outer_radius_m=radius,closed_source_sha256=m.geometry_sha,inside_actual_source_volume=True,perfusion_bed_modeled=False)
        report['selected']=witness;target_reports.append(report)
        return points[k],witness

    palate,palate_witness=recipient('palatopharyngeus-left',[.007,1.6585,.012],.00055)
    # Constrictor perforation is explicitly intended only at the tonsillar
    # recipient.  Neither bone nor medial-pterygoid crossings are waived.
    tonsil,tonsil_witness=recipient('palatine-tonsil-left',[.028,1.663,.033],.00055)
    gland0,gland0_witness=recipient('submandibular-gland-left',[.030,1.622,.038],.00050)
    gland1,gland1_witness=recipient('submandibular-gland-left',[.025,1.614,.051],.00040)
    sc=meshes['superior-constrictor-left'];tri=sc.triangles;face_centers=tri.mean(1)
    # Actual lateral and medial source face witnesses define the change of
    # route topology.  The global bounding-box maximum is not called a border.
    upper=[]
    for x in [.024,.008]:
        ids=np.flatnonzero(abs(face_centers[:,0]-x)<.002)
        assert len(ids)
        ids=ids[np.argsort(face_centers[ids,1])[-20:]]
        face=int(ids[np.argmax(face_centers[ids,1])]);point=face_centers[face]
        upper.append(dict(source_part_id=roles['superior-constrictor-left'],face=face,surface_point_m=point.tolist(),authored_clearance_waypoint_m=(point+[0,.0015,0]).tolist(),qualification='Source upper surface patch in an explicit lateral/medial band; anatomical border footprint remains provisional, so actual packing/sections are required.'))
    # Source-sided corridor witness: choose points between SG and SP using
    # actual source surfaces, independently of the unregistered private GLB.
    sg,sp=meshes['styloglossus-left'],meshes['stylopharyngeus-left'];v=sg.triangles.mean(1);ids=np.flatnonzero((v[:,1]>1.641)&(v[:,1]<1.655)&(v[:,2]<.025));q=v[ids]
    d2,faces,near=igl.point_mesh_squared_distance(q,sp.vertices,sp.faces);mid=(q+near)/2
    valid=(np.sqrt(d2)>.0015)&(np.sqrt(d2)<.012)
    valid&=(abs(igl.winding_number(sg.vertices,sg.faces,mid))<.01)&(abs(igl.winding_number(sp.vertices,sp.faces,mid))<.01)
    corridor_ids=np.flatnonzero(valid)
    interval=None
    if len(corridor_ids):
        k=corridor_ids[np.argmin(np.linalg.norm(mid[corridor_ids]-[.035,1.647,.008],axis=1))]
        interval=dict(styloglossus_face=int(ids[k]),stylopharyngeus_face=int(faces[k]),source_surface_points_m=[q[k].tolist(),near[k].tolist()],authored_interval_point_m=mid[k].tolist(),surface_separation_m=float(np.sqrt(d2[k])))
    contract=dict(schema='anatomical-left-facial-cervical-construction-contract/v5',source_catalogue=str(source.relative_to(ROOT)),source_catalogue_sha256=EXPECTED,
                  authority_inventory='reports/ta98-male/facial-arterial-family-contract-v1.json',authority_inventory_sha256=digest_file(ROOT/'reports/ta98-male/facial-arterial-family-contract-v1.json'),
                  primary_references=[dict(url='https://www.scielo.cl/pdf/ijmorphol/v30n3/art14.pdf',actual_scope='Maistry et al.2012:100 sides from30 fetal and20 adult heads; adult Fig2 onp849 explicitly shows ascending superficial to SC and descending over its superior border. It establishes course topology, not individual target dimensions.'),dict(url='https://pubmed.ncbi.nlm.nih.gov/31291033/',actual_scope='44 hemi-sections; parapharyngeal facial bulb and tonsillar branches. SG/SP interval is frequent; placement remains source-specific.')],
                  private_reference_data_policy='Description relationships used only; no private surface/course or unverified GLB coordinate registration.',
                  family_terms=old['terms'],fixed_sources=list(records.values()),native_upper_margin_witnesses=upper,sg_sp_interval_witness=interval,
                  recipient_profile_feasibility=target_reports,parent_departure_policy='Ascending palatine station2 near actual cervical facial origin; tonsillar station8 in same proximal segment; glandular stations32 and44 adjoining native lateral/inferior gland groove.',
                  legitimate_interfaces=[dict(branch='tonsillar',recipient='ZA:Superior pharyngeal constrictor.l',kind='perforating arterial entry',support='Final localized passage toward native inferior tonsillar tissue; whole muscular entry is not categorically forbidden.'),dict(branch='glandular-0/glandular-1',recipient='ZA:Submandibular gland.l',kind='adjacent arterial wall/tissue entry')],
                  fixed_exclusions=['whole skeleton and mandible','veins and unrelated arteries','submandibular duct','pharyngeal/oral mucosal luminal boundaries'],
                  wall_and_small_calibre_illustrative=True,native_parent_outer_fixed=True,root_and_terminal_section_calibres_not_reduced=True,
                  pharyngeal_airspace_full_volume_unverified=True,native_mucosal_triangle_checks_mandatory=True,installed=False,term_count_added=0)
    atomic_json(ROOT/'reports/ta98-male/facial-anatomical-construction-contract-v5.json',contract)
    atomic_json(D/'construction-contract.json',contract)
    if any(x is None for x in [palate,tonsil,gland0,gland1])or interval is None:
        report=dict(stopped_before_guides=True,reason='No complete recipient profile ball or SG/SP interval in required fixed native source tissues',recipient_feasibility=target_reports,sg_sp_interval_witness=interval,installed=False,elapsed_seconds=time.monotonic()-began)
        atomic_json(D/'stage1-evidence.json',report);print(json.dumps(report));return
    main_interval=np.array(interval['authored_interval_point_m']);top0,top1=[np.array(x['authored_clearance_waypoint_m'])for x in upper]
    planned=[dict(name='ascending-palatine',term='A12.2.05.021',station=2,radius=.00055,points=[centers[2],[.038,1.637,.003],main_interval,top0,top1,[.004,1.675,.009],palate],recipient=palate_witness),
             dict(name='tonsillar',term='A12.2.05.022',station=8,radius=.00055,points=[centers[8],[.039,1.639,.008],main_interval+[.002,0,.004],[.034,1.658,.026],tonsil],recipient=tonsil_witness),
             dict(name='glandular-0',term='A12.2.05.024',station=32,radius=.00050,points=[centers[32],centers[32]+[.0005,-.002,-.0015],[.030,1.617,.033],gland0],recipient=gland0_witness),
             dict(name='glandular-1',term='A12.2.05.024',station=44,radius=.00040,points=[centers[44],centers[44]+[0,-.003,-.001],[.030,1.611,.047],gland1],recipient=gland1_witness)]
    routes=[]
    for row in planned:
        points,_=centripetal_curve(row['points'],chord_tolerance_m=.00001,max_step_m=.0005)
        path=D/(row['name']+'-guide.npz');atomic_npz(path,points=points,outer_radii=np.full(len(points),row['radius']),inner_radii=np.full(len(points),row['radius']-.00009),frame=np.array(old['frame']))
        routes.append(dict(name=row['name'],term=row['term'],parent='facial-left',parent_station=row['station'],guide=str(path.relative_to(ROOT)),guide_sha256=digest_file(path),stations=len(points),length_m=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum()),illustrative_radius_m=row['radius'],illustrative_wall_m=.00009,recipient=row['recipient']))
    fresh_inputs=dict(old,source_catalogue_sha256=EXPECTED,fixed_sources=list(records.values()));atomic_json(D/'inputs.json',fresh_inputs)
    plan=dict(source_catalogue=str(source.relative_to(ROOT)),source_catalogue_sha256=EXPECTED,source_inputs=str((D/'inputs.json').relative_to(ROOT)),source_inputs_sha256=digest_file(D/'inputs.json'),frame=old['frame'],fixed_parent=records['facial-left'],additional_fixed_sources=[],routes=routes,illustrative_parent_wall_m=.00015,
              construction_contract=str((D/'construction-contract.json').relative_to(ROOT)),construction_contract_sha256=digest_file(D/'construction-contract.json'),source_geometry_modified=False,installed=False)
    atomic_json(D/'routes.json',plan)
    report=dict(source_catalogue_sha256=EXPECTED,routes=[dict(name=r['name'],station=r['parent_station'],length_mm=r['length_m']*1000)for r in routes],upper_margin_witnesses=upper,sg_sp_interval=interval,all_terminal_profile_balls_clear=True,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(D/'stage1-evidence.json',report);assert digest_file(source)==EXPECTED;print(json.dumps(report))


if __name__=='__main__':main()
