"""First native-frame cervical facial branch courses and tissue recipients."""
from pathlib import Path
import json
import sys
import time
import argparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import audit, topology
from geometry_tools.algorithms.tube_sections import centripetal_curve


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', default='v1')
    args = parser.parse_args()
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/circulatory-completion'
    previous = base / 'male-facial-family-inputs-v1/inputs.json'
    inputs = json.loads(previous.read_bytes())
    source = ROOT / inputs['source_catalogue']
    assert digest_file(source) == inputs['source_catalogue_sha256']
    directory = base / ('male-left-facial-cervical-routes-' + args.version)
    assert not directory.exists()
    directory.mkdir()
    records = {r['role']: r for r in inputs['fixed_sources']}
    meshes = {k: Mesh.load(ROOT / r['mesh']) for k, r in records.items() if not r.get('missing_native_source')}
    catalogue = json.loads(source.read_bytes())
    parts = {p['id']: p for p in catalogue['parts']}
    reader = AtlasReader()
    for role, pid in [('styloglossus-left', 'TA98TONGUE:male:styloglossus:left'), ('palatopharyngeus-left', 'ZA:Palatopharyngeus muscle.l')]:
        mesh = reader.read_part(catalogue, parts[pid], frame=inputs['frame'], catalogue_path=source, catalogue_sha=inputs['source_catalogue_sha256'])
        target = directory / (role + '.npz')
        mesh.save(target)
        meshes[role] = mesh
        records[role] = dict(role=role, part_id=pid, mesh=str(target.relative_to(ROOT)), mesh_sha256=digest_file(target), geometry_sha256=mesh.geometry_sha, source_record=parts[pid])
    parent = meshes['facial-left']
    parent_audit = audit(parent, expected_components=1)
    atomic_json(directory / 'native-parent-audit.json', parent_audit)
    assert parent_audit['passed'], parent_audit['failures']
    native = records['facial-left']
    stations = np.asarray(native['source_station_components'][0]['rings'])
    centers = parent.vertices[stations].mean(1)

    def recipient(role, seed, radius):
        mesh = meshes[role]
        gate = topology(mesh)
        assert gate['closed_oriented_manifold'], (role, gate)
        seed = np.asarray(seed)
        low, high = mesh.bounds
        points = np.stack(np.meshgrid(*[np.linspace(max(low[i], seed[i] - .003), min(high[i], seed[i] + .003), 11) for i in range(3)], indexing='ij'), -1).reshape(-1, 3)
        distance = np.sqrt(igl.point_mesh_squared_distance(points, mesh.vertices, mesh.faces)[0])
        winding = np.abs(igl.winding_number(mesh.vertices, mesh.faces, points))
        valid = (winding > .99) & (distance > radius + .00015)
        assert valid.any(), (role, 'No complete terminal profile fits within actual source tissue')
        ids = np.flatnonzero(valid)
        k = ids[np.argmin(np.linalg.norm(points[ids] - seed, axis=1))]
        return points[k], dict(recipient_role=role, point_m=points[k].tolist(), minimum_surface_distance_m=float(distance[k]), terminal_outer_radius_m=radius,
                               closed_source_sha256=mesh.geometry_sha, inside_actual_source_volume=True, perfusion_bed_modeled=False)

    if args.version == 'v1':
        palate, palate_evidence = recipient('soft-palate', [.013, 1.659, .027], .00055)
    else:
        # v1 proved the original generic palate envelope is a ~0.5 mm slab:
        # the same 1.1 mm profile cannot be placed inside it.  Keep calibre;
        # end this initial artery in its actual upper palatal-muscle recipient,
        # retaining the broader soft-palate perfusion interface as unfinished.
        palate, palate_evidence = recipient('palatopharyngeus-left', [.007, 1.6585, .012], .00055)
        palate_evidence['scope'] = 'Upper native palatal-muscle recipient only; generic thin Soft palate envelope cannot contain this profile. Additional soft-palate distribution remains unmodeled.'
    tonsil, tonsil_evidence = recipient('palatine-tonsil-left', [.030, 1.664, .032], .00055)
    gland0, gland0_evidence = recipient('submandibular-gland-left', [.028, 1.627, .043], .00050)
    gland1, gland1_evidence = recipient('submandibular-gland-left', [.026, 1.619, .047], .00040)
    planned = [
        dict(name='ascending-palatine', term='A12.2.05.021', parent_station=4, radius=.00055, waypoints=[[.035, 1.640, .013], [.030, 1.650, .016], [.019, 1.657, .014], palate], recipient=palate_evidence),
        dict(name='tonsillar', term='A12.2.05.022', parent_station=16, radius=.00055, waypoints=[[.041, 1.639, .024], [.039, 1.651, .028], [.033, 1.661, .031], tonsil], recipient=tonsil_evidence),
        dict(name='glandular-0', term='A12.2.05.024', parent_station=27, radius=.00050, waypoints=[[.032, 1.628, .034], gland0], recipient=gland0_evidence),
        dict(name='glandular-1', term='A12.2.05.024', parent_station=37, radius=.00040, waypoints=[[.030, 1.618, .043], gland1], recipient=gland1_evidence),
    ]
    routes = []
    for row in planned:
        points, _ = centripetal_curve([centers[row['parent_station']], *row['waypoints']], chord_tolerance_m=.00001, max_step_m=.0005)
        path = directory / (row['name'] + '-guide.npz')
        atomic_npz(path, points=points, outer_radii=np.full(len(points), row['radius']), inner_radii=np.full(len(points), row['radius'] - .00009), frame=np.array(inputs['frame']))
        routes.append(dict(name=row['name'], term=row['term'], parent='facial-left', parent_station=row['parent_station'], guide=str(path.relative_to(ROOT)), guide_sha256=digest_file(path),
                           stations=len(points), length_m=float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum()), illustrative_radius_m=row['radius'], illustrative_wall_m=.00009,
                           recipient=row['recipient']))
    result = dict(schema='source-bound-left-facial-cervical-route-contract/v1', source_catalogue=inputs['source_catalogue'], source_catalogue_sha256=inputs['source_catalogue_sha256'],
                  inventory_contract='reports/ta98-male/facial-arterial-family-contract-v1.json', inventory_contract_sha256=digest_file(ROOT / 'reports/ta98-male/facial-arterial-family-contract-v1.json'),
                  source_inputs=str(previous.relative_to(ROOT)), source_inputs_sha256=digest_file(previous), frame=inputs['frame'], fixed_parent=native, additional_fixed_sources=[records[k] for k in ['styloglossus-left', 'palatopharyngeus-left']],
                  planned_parts=['native-facial-parent', *[r['name'] for r in routes]], routes=routes,
                  all_four_branch_mouths='Actual parent inner/outer side openings, no retained internal classifier caps; common source profile properly transported along each smooth authored course',
                  native_parent_outer_coordinates_fixed=True, illustrative_parent_wall_m=.00015,
                  tissue_interfaces=['ascending-palatine terminal within actual upper native palatal muscle; wider generic palate distribution unfinished', 'tonsillar passage through superior constrictor to native inferior tonsillar tissue', 'two glandular mouths within same-side submandibular gland'],
                  exclusions=['No oral/pharyngeal lumen entry', 'No new nerves/bone/eyes/vein crossings', 'No geometry edits to recipients or global protected anatomy'],
                  residual_interfaces=inputs['residual_upstream_interfaces'], primary_references=inputs['primary_references'],
                  calibre_limitation='Cervical small-branch radii and wall thickness illustrative; study supplies anatomical course relationships but does not establish these individual dimensions. Calibre is held through bends.',
                  missing_distal_family_terms=['A12.2.05.027', 'A12.2.05.028'], source_geometry_modified=False, installed=False, anatomical_acceptance=False, elapsed_seconds=time.monotonic()-began)
    atomic_json(directory / 'routes.json', result)
    atomic_json(ROOT / ('reports/ta98-male/facial-cervical-construction-contract-' + args.version + '.json'), result)
    assert digest_file(source) == inputs['source_catalogue_sha256']
    print(json.dumps(dict(routes=[dict(name=r['name'], length_mm=r['length_m']*1000, terminal=r['recipient']['point_m'], recipient_clearance_mm=r['recipient']['minimum_surface_distance_m']*1000) for r in routes],
                          parent_self_crossings=parent_audit['self_intersections']['defect_pairs'], elapsed_seconds=result['elapsed_seconds'])))


if __name__ == '__main__':
    main()
