"""Freeze actual left basal venous mouths and ownership before one construction."""
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader

BASE = Path('.local-models/geometry-workbench/circulatory-completion')
OUT = BASE / 'male-intrapulmonary-circulation-current-v1/collector-left-common-basal-v1'
SOURCE = Path('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')
SOURCE_SHA = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
FRAME = 'model:local-ta98:assembled-metres'


def main():
    began = time.monotonic()
    OUT.mkdir(exist_ok=True)
    assert not (OUT / 'contract-before-construction.json').exists()
    assert digest_file(SOURCE) == SOURCE_SHA
    atlas = json.loads(SOURCE.read_text())
    parts = {p['id']: p for p in atlas['parts']}
    old = BASE / 'pulmonary-circulation-frontier-v1'
    curves = {c['index']: c for c in json.loads((old / 'native-curves-v3/evidence.json').read_text())['curves']}
    sections = np.load(old / 'native-curves-v3/literal-section-centres.npz', allow_pickle=False)
    reader = AtlasReader()
    host = reader.read_part(atlas, parts['ZA:Left inferior pulmonary vein'], frame=FRAME)
    host.save(OUT / 'fixed-current-lipv.npz')
    plan_path = BASE / 'male-lipv-native-port-restoration-v1/plan.json'
    port = json.loads(plan_path.read_text())['ports'][0]
    tree = cKDTree(host.vertices)
    arrays = {}
    mouth = {}
    for side in ['outer', 'inner']:
        error, ids = tree.query(np.asarray(port['target_' + side]))
        assert float(error.max()) < 3e-9
        q = host.vertices[ids]
        arrays['parent_' + side] = q
        arrays['parent_' + side + '_compact_ids'] = ids
        arrays['parent_' + side + '_raw_ids'] = np.asarray(host.source['source_vertex_ids'])[ids]
        arrays['parent_' + side + '_effective_normal_vectors'] = host.attributes['normals'][ids] * host.attributes['normalSigns'][ids, None]
        arrays['parent_' + side + '_uvs'] = host.attributes['uvs'][ids]
        mouth[side] = dict(compact_ids=ids.tolist(), raw_ids=arrays['parent_' + side + '_raw_ids'].tolist(), storage_target_error_m=float(error.max()))
    capids = np.asarray(port['annulus_source_faces'], int)
    arrays['parent_annulus_source_face_ids'] = capids
    arrays['parent_annulus_oriented_corners'] = host.vertices[host.faces[capids]]
    arrays['parent_annulus_compact_faces'] = host.faces[capids]
    center = arrays['parent_outer'].mean(0)
    normal = np.linalg.svd(arrays['parent_outer'] - center, full_matrices=False)[2][-1]
    direction = sections[curves[131]['key'] + '_centers'][5] - center
    if normal @ direction < 0:
        normal *= -1
    arrays['parent_center'] = center
    arrays['parent_distal_normal'] = normal
    source_meshes = {}
    rows = []
    graph = {60: (131, 2), 61: (60, 84), 62: (60, 72), 63: (62, 24), 132: (131, 36), 133: (132, 24), 134: (131, 60)}
    chosen = [60, 61, 62, 63, 131, 132, 133, 134]
    for index in chosen:
        curve = curves[index]
        ident = curve['part_id']
        if ident not in source_meshes:
            source_meshes[ident] = reader.read_part(atlas, parts[ident], frame=FRAME)
            source_meshes[ident].save(OUT / ('source-' + ('inferior' if 'Inferior' in ident else 'superior') + '-basal.npz'))
        mesh = source_meshes[ident]
        rings = sections[curve['key'] + '_literal_sections'].copy()
        root = 0 if curve['larger_radius_endpoint'] == 0 else len(rings) - 1
        # Whole source profiles are held. Only the explicitly crossing root collar
        # of curve60 is replaced; first fully distal profile is derived, not fitted.
        signed_min = ((rings - center) @ normal).min(1)
        trim = 0
        if index == 60:
            assert root == 0
            trim = int(np.flatnonzero(signed_min > 1e-7)[0])
            assert trim == 8
        points = rings.reshape(-1, 3)
        original_faces = sections[curve['key'] + '_source_face_ids']
        error, mapped = cKDTree(points).query(mesh.vertices[mesh.faces[original_faces]].reshape(-1, 3))
        assert float(error.max()) == 0
        f = mapped.reshape(-1, 3)
        keep = np.all(f // 12 >= trim, axis=1)
        f = f[keep] - 12 * trim
        retained = rings[trim:].copy()
        source_unchanged = retained.copy()
        if index == 131:
            # The two circular rims have a real polygon-phase mismatch. Replace
            # only native strip0 with an exact-rim bridge; profile1+ is literal.
            dist, correspondence = cKDTree(retained[0]).query(arrays['parent_outer'])
            assert len(set(correspondence)) == 12 and float(dist.max()) < .0002
            retained[0, correspondence] = arrays['parent_outer']
        inner = retained.mean(1)[:, None] + .8 * (retained - retained.mean(1)[:, None])
        if index == 131:
            inner[0, correspondence] = arrays['parent_inner']
        key = 'curve_' + str(index)
        arrays[key + '_outer'] = retained.reshape(-1, 3)
        arrays[key + '_inner'] = inner.reshape(-1, 3)
        arrays[key + '_faces'] = f
        arrays[key + '_source_face_ids'] = original_faces[keep]
        arrays[key + '_retained_native_sections'] = source_unchanged
        anchor = None
        if index in graph:
            parent, station = graph[index]
            anchor = sections[curves[parent]['key'] + '_centers'][station]
            arrays[key + '_internal_anchor'] = anchor
        rows.append(dict(curve=index, part_id=ident, native_component=curve['component'], key=key,
                         source_geometry_sha256=mesh.geometry_sha, source_bank_sha256=mesh.source['bank_sha256'],
                         removed_root_stations=list(range(trim)), replaced_bridge_strip_ids=[0] if index == 131 else [], retained_first_station=trim,
                         original_profile_count=len(rings), retained_faces=len(f), native_root_endpoint=root,
                         retained_root_endpoint=0 if root == 0 else 1,
                         parent_curve=graph[index][0] if index in graph else None,
                         parent_station=graph[index][1] if index in graph else None,
                         internal_anchor_m=None if anchor is None else anchor.tolist(),
                         source_control_motion_max_m=float(np.linalg.norm(retained-source_unchanged, axis=2).max()),
                         root_profile_min_distal_plane_m=float(signed_min[trim]),
                         outer_source_calibre_retained=True))
    archive = OUT / 'construction-inputs.npz'
    np.savez_compressed(archive, **arrays)
    contract = dict(schema='left-common-basal-collector-construction-contract/v1', source_catalogue_sha256=SOURCE_SHA,
                    frame=FRAME, target_term='A12.3.02.405', target_identity='common basal venous confluence upstream of superior/inferior basal tributaries, not V6',
                    authoritative_inventory='male-intrapulmonary-circulation-current-v1/venous-inventory-v2.json',
                    authoritative_inventory_sha256=digest_file(OUT.parent / 'venous-inventory-v2.json'),
                    input_archive_sha256=digest_file(archive), parent_part_id='ZA:Left inferior pulmonary vein',
                    parent_source_record=parts['ZA:Left inferior pulmonary vein'], parent_geometry_sha256=host.geometry_sha,
                    parent_bank_sha256=host.source['bank_sha256'], parent_annulus=mouth,
                    parent_annulus_source_face_ids=capids.tolist(), component_graph=rows,
                    whole_source_cohort='All8 source components in both left basal source objects; all other83-component-family anatomy fixed.',
                    physical_wall='Outer union minus one connected inner blood union. 0.8 radial fraction follows existing illustrative root convention; no histological thickness or TA98 coat claim.',
                    physical_ownership='First repair complete communicating finite basal wall; root confluence .405 and downstream .406/.410 ownership is provisional until actual shared-cell partition. No child term is an alias of the whole network.',
                    parent_fixed='All5960 parent faces, every source position/normal/UV/sign and all heart-facing interfaces unchanged. Its basal annulus is the explicit reciprocal attachment; V6 annulus unchanged.',
                    source_preservation='Curve60rootprofiles0..7 are consumed by the junction collar. Curve131nativefirststrip0 is replaced by an exact bridge between current canonical parent rim and literal sourceprofile1; source root polygon phase differs0.173275mm. Every other native exterior control is literal. No radius shrink.',
                    fixed_skeleton=True, fixed_heart=True, fixed_bronchi=True, fixed_lungs=True,
                    fixed_other_vessels=True, no_registration_or_pose_change=True,
                    mandatory_gates=['One connected closed finite wall and one connected capped blood classifier; positive wall volume; actual self0 and no lumen disk.', 'Nine actual external blood ports: existing parent annulus plus8native distal mouths, source calibre held.', 'Full source-field replay outside declared collar, root reciprocal coords/normalwords/UV and fixed V6 mouth.', 'Independent fixed lung/bronchus/heart/arterial context review; source congenital variations not normalized by proximity.'],
                    primary_sources=['https://pmc.ncbi.nlm.nih.gov/articles/PMC7139082/', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7715996/', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9735388/'],
                    no_intra_intersegmental_identity_from_nearest_plane=True, installed=False,
                    anatomical_acceptance=False, elapsed_seconds=time.monotonic()-began)
    path = OUT / 'contract-before-construction.json'
    path.write_text(json.dumps(contract, indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT), contract_sha256=digest_file(path), inputs_sha256=digest_file(archive),
                          parent_annulus=mouth, components=len(rows), elapsed_seconds=contract['elapsed_seconds'])))


if __name__ == '__main__':
    main()
