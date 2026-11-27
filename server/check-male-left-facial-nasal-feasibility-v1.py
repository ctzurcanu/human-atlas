"""Bounded, read-only-to-atlas nasal arterial source and interface feasibility."""
from pathlib import Path
import json
import sqlite3
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import audit

SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'
BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-left-facial-nasal-v1'
INVENTORY = ROOT / 'reports/ta98-male/facial-arterial-family-contract-v1.json'


def main():
    start = time.monotonic()
    assert digest_file(SOURCE) == EXPECTED
    # The first saved run stopped before the distance query; only that exact
    # validated source-export prefix may be resumed, never a completed report.
    assert not (BASE/'feasibility.json').exists()
    BASE.mkdir(exist_ok=True)
    old = json.loads(INVENTORY.read_bytes())
    atlas = json.loads(SOURCE.read_bytes())
    parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    rows = {x['role']: x for x in old['fixed_sources'] if not x.get('missing_native_source')}
    roles = ['facial-left', 'angular-left', 'superior-labial-left', 'nasal-septal-cartilage', 'nasal-mucosa', 'nasalis-left', 'nasal-skin-left']
    source_rows = []
    meshes = {}
    for role in roles:
        row = rows[role]
        mesh = reader.read_part(atlas, parts[row['part_id']], frame=old['frame'], catalogue_path=SOURCE, catalogue_sha=EXPECTED)
        assert mesh.geometry_sha == row['geometry_sha256'], ('Source changed', role)
        output=BASE/(role+'.npz')
        if output.exists():
            assert Mesh.load(output).geometry_sha==mesh.geometry_sha
        else:
            mesh.save(output)
        meshes[role] = mesh
        source_rows.append(dict(role=role, part_id=row['part_id'], geometry_sha256=mesh.geometry_sha,
                                old_geometry_identical=True, mesh=str((BASE / (role+'.npz')).relative_to(ROOT))))
    extras = {'nasal-skin-lateral-left': "ZA:Nasal region.l:Skin-5'", 'alar-cartilage-left': 'ZA:Major alar cartilage.l',
              'lateral-septal-process-left': 'ZA:Lateral process of nasal septal cartilage.l'}
    for role, pid in extras.items():
        mesh = reader.read_part(atlas, parts[pid], frame=old['frame'], catalogue_path=SOURCE, catalogue_sha=EXPECTED)
        output=BASE/(role+'.npz')
        if output.exists():
            assert Mesh.load(output).geometry_sha==mesh.geometry_sha
        else:
            mesh.save(output)
        meshes[role] = mesh
        source_rows.append(dict(role=role, part_id=pid, geometry_sha256=mesh.geometry_sha,
                                mesh=str((BASE / (role+'.npz')).relative_to(ROOT))))
    authority = sqlite3.connect('file:' + str(ROOT / '.local-models/ta98-authority.sqlite') + '?mode=ro', uri=True)
    term_rows = [dict(zip(['id','name_en','name_la','parent_id'], x)) for x in authority.execute("select id,name_en,name_la,parent_id from ta98 where id in ('A12.2.05.027','A12.2.05.028') order by id")]
    assert [x['parent_id'] for x in term_rows] == ['A12.2.05.026','A12.2.05.020']
    contract = dict(schema='male-left-facial-nasal-contract/v1', source_catalogue_sha256=EXPECTED,
        source_catalogue=str(SOURCE.relative_to(ROOT)), frame=old['frame'],
        whole_family_inventory=str(INVENTORY.relative_to(ROOT)), whole_family_inventory_sha256=digest_file(INVENTORY),
        all_ten_authority_terms=old['terms'], planned_terms=term_rows, fixed_sources=source_rows,
        planned_parts=[dict(term='A12.2.05.027', type='finite-vascular-tissue', parent='ZA:Superior labial artery.l', recipient='left anterior-inferior septal mucoperichondrium', airway_excluded=True),
                       dict(term='A12.2.05.028', type='finite-vascular-tissue', parent='distal ZA:Facial artery.l / proximal continuing ZA:Angular artery.l component 0', recipient='left lateral nose and ala tissue', airway_excluded=True)],
        interface_rules=['Actual common blood side opening and shared reciprocal tissue faces required; mere cylinder overlap rejected.',
                         'Cartilage, skeleton, skin and native recipient tissue fixed; target artery belongs in tissue, never airway.',
                         'Angular seven-component carrier is not a single branch and will not be relabeled wholesale.',
                         'All far native source faces remain literal; every changed owner must retain face ancestry.'],
        fixed_global=['skeleton','heart','both lingual families','both cystic families','all native nasal tissues'],
        primary_references=[dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11548982/', type='60-hemiface human cadaveric dissection', use='Lateral nasal origin beside nose/alar base; observed mean diameter 0.74mm, length 10.6mm; no exact patient coordinates'),
                            dict(url='https://pubmed.ncbi.nlm.nih.gov/27568352/', type='Injected fresh-frozen cadaveric dissection', use='Septal superior-labial contribution to anterior septal vascular supply, distinct from posterior septal artery')],
        source_rights='Native Z-Anatomy geometry CC BY-SA 4.0, retained local. Supplied private geometry is reference only and not copied.',
        installed=False, acceptance=False)
    for output in [BASE/'contract.json',ROOT/'reports/ta98-male/left-facial-nasal-contract-v1.json']:
        if output.exists():
            assert json.loads(output.read_bytes())==contract
        else:
            atomic_json(output,contract)
    angular = meshes['angular-left']
    components = rows['angular-left']['source_station_components']
    bvh = {k:BVH(meshes[k]) for k in ['nasal-skin-left','nasal-skin-lateral-left','nasalis-left','nasal-septal-cartilage','alar-cartilage-left','nasal-mucosa']}
    main_component = Mesh(angular.vertices, angular.faces[np.asarray(components[0]['source_faces'])], angular.frame)
    checks=[]
    for i, comp in enumerate(components):
        rings=np.asarray(comp['rings'])
        centers=angular.vertices[rings].mean(1)
        mesh=Mesh(angular.vertices,angular.faces[np.asarray(comp['source_faces'])],angular.frame)
        output=BASE/f'angular-component-{i}.npz'
        if output.exists():
            assert Mesh.load(output).geometry_sha==mesh.geometry_sha
        else:
            mesh.save(output)
        near={}
        for target, tree in bvh.items():
            a,_=tree.closest_points(centers)
            minimum=int(np.argmin(a['distance_m']))
            near[target]=dict(minimum_center_distance_m=float(a['distance_m'][minimum]), station=minimum,
                              artery_center_m=centers[minimum].tolist(), closest_tissue_point_m=a['points'][minimum].tolist(),
                              target_face=int(a['face_ids'][minimum]))
        parent=None
        if i:
            a,p=pair_contacts(mesh,main_component,0.,max_pairs=10000,max_tests=300000)
            parent=dict(proof=p,first_three_pairs=a['pairs'][:3].tolist())
        checks.append(dict(component=i, source_faces=len(comp['source_faces']), stations=len(rings),
                           end_centers_m=centers[[0,-1]].tolist(), length_m=float(np.linalg.norm(np.diff(centers,axis=0),axis=1).sum()),
                           diameters_m=(2*np.linalg.norm(angular.vertices[rings]-centers[:,None],axis=2).mean(1)[[0,-1]]).tolist(),
                           target_proximity=near, parent_contact=parent))
    topology={k:audit(meshes[k],self_intersections=False) for k in ['nasal-mucosa','nasal-skin-left','nasal-skin-lateral-left','nasal-septal-cartilage','alar-cartilage-left','superior-labial-left']}
    result=dict(schema='male-left-facial-nasal-feasibility/v1',source_catalogue_sha256=EXPECTED,
                angular_components=checks,source_topology=topology,
                septal_target_bounds=meshes['nasal-mucosa'].bounds.tolist(),
                left_septal_mucosa_surface_present=bool(meshes['nasal-mucosa'].vertices[:,0].max()>0.001),
                source_native_parent_blood_space=False, source_native_parent_status='Open single source skins; no independently established native blood lumen or finite wall.',
                installed=False, installable=False, elapsed_seconds=time.monotonic()-start)
    atomic_json(BASE/'feasibility.json',result)
    assert digest_file(SOURCE)==EXPECTED
    print(json.dumps(dict(output=str(BASE), angular_component_count=len(checks), left_septal_mucosa_surface_present=result['left_septal_mucosa_surface_present'],
                          closest_lateral_tissue=[dict(component=x['component'],length_mm=x['length_m']*1000,nasalis_mm=x['target_proximity']['nasalis-left']['minimum_center_distance_m']*1000,skin_mm=x['target_proximity']['nasal-skin-lateral-left']['minimum_center_distance_m']*1000) for x in checks],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
