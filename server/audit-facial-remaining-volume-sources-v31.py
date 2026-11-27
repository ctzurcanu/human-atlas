"""Read-only literal-source classifier audit; never fill anatomical openings."""
from pathlib import Path
import importlib.util
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import manifold3d as md
from scipy.spatial import ConvexHull
from geometry_tools.audit import topology, audit
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh


def compact_topology(t):
    return {k: (len(v) if k == 'invalid_vertex_links' else v) for k, v in t.items()}


def main():
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/circulatory-completion'
    work = base / 'male-left-facial-cervical-finite-network-cervical-c2-v17'
    out = work / 'remaining-physical-source-audit-v31'
    assert not out.exists()
    out.mkdir()
    binding_path = work / 'source-binding-v19/evidence.json'
    binding = json.loads(binding_path.read_bytes())
    source = ROOT / binding['snapshot']
    expected = 'd73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
    revision = '31fab36e6897adaec9c8ddca201261607902203d75c64f3b214e57d329f460c2'
    assert digest_file(source) == expected
    assert digest_file(ROOT / binding['current_catalogue']) == expected
    atlas = json.loads(source.read_bytes())
    reader = AtlasReader()
    wall = Mesh.load(work / 'finite-wall.npz')
    origin = wall.bounds.mean(0)
    helper_spec = importlib.util.spec_from_file_location('facial_context_helper_v31', ROOT / 'server/check-facial-cervical-feasibility.py')
    helper = importlib.util.module_from_spec(helper_spec)
    helper_spec.loader.exec_module(helper)
    normal_spec = importlib.util.spec_from_file_location('source_exact_normalization_v31', ROOT / 'server/normalize-facial-source-classifiers-v26.py')
    normal = importlib.util.module_from_spec(normal_spec)
    normal_spec.loader.exec_module(normal)
    plan = json.loads((base / 'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes())
    parent, _ = helper.occupancy_classifier(Mesh.load(ROOT / plan['fixed_parent']['mesh']), True)
    assert parent is not None
    parent_solid = helper.solid(parent, origin)
    occupied = {}
    field_bindings = []
    for route in plan['routes']:
        field_path = work / (route['name'] + '-fields.npz')
        z = np.load(field_path)
        outer = Mesh(z['fields'][0], z['faces'], wall.frame)
        classifier, mode = helper.occupancy_classifier(outer, True)
        assert classifier is not None, (route['name'], mode)
        occupied[route['name']] = helper.solid(classifier, origin) - parent_solid
        field_bindings.append(dict(branch=route['name'], path=str(field_path.relative_to(ROOT)), sha256=digest_file(field_path), classifier_mode=mode))
    context = json.loads((work / 'registered-context-v25/evidence.json').read_bytes())
    prior_native = json.loads((work / 'native-classifiers-v20/evidence.json').read_bytes())
    objects = ['Mandible', 'Superficial investing cervical fascia.l', 'Mylohyoid muscle.ol']
    results = []
    for name in objects:
        assert time.monotonic() - began < 45, '45-second saved-stage budget'
        parts = [p for p in atlas['parts'] if p.get('sourceId') == name and not p.get('suppressed')]
        assert parts, name
        vertices, faces, members = [], [], []
        offset = 0
        for p in parts:
            mesh = reader.read_part(atlas, p, frame=wall.frame, catalogue_path=source, catalogue_sha=expected)
            members.append(dict(part_id=p['id'], source_id=p.get('sourceId'), system=p['system'], material=p.get('material'), groups=p.get('groups'), geometry_sha256=mesh.geometry_sha, attributes_sha256=mesh.attributes_sha, topology=compact_topology(topology(mesh)), explicit_pose={k:p.get(k) for k in ('geometryOrigin','sourceOffset','geometryScale','geometryInstance')}))
            vertices.append(mesh.vertices)
            faces.append(mesh.faces + offset)
            offset += len(mesh.vertices)
        raw = np.concatenate(vertices)
        v, inverse = np.unique(raw, axis=0, return_inverse=True)
        f = inverse[np.concatenate(faces)]
        components, normalization = normal.exact_components(v, f)
        checks = []
        complete_closed = True
        for used, cf, source_faces in components:
            cm = Mesh(v[used], cf, wall.frame)
            t = topology(cm)
            c = dict(geometry_sha256=cm.geometry_sha, topology=compact_topology(t), source_triangle_coordinates_literal=bool(np.array_equal(cm.triangles, v[f[source_faces]])))
            if t['closed_oriented_manifold']:
                a = audit(cm, require_closed=True, max_pairs=1000000, max_tests=12000000)
                c['audit'] = a
                complete_closed &= bool(a['passed'])
            else:
                complete_closed = False
            checks.append(c)
        row = dict(source_object_id=name, source_members=members, joined_geometry_sha256=Mesh(v, f, wall.frame).geometry_sha, normalization={k: (len(value) if isinstance(value,list) else value) for k,value in normalization.items()}, components=checks, complete_literal_closed_classifier=complete_closed, native_caps_or_repairs=0, native_arrays_unchanged=True)
        if name == 'Mylohyoid muscle.ol':
            assert all(p['system'] == 'attachments' and p['material'] == 'ZA:attachments:Origin mastication' and '2: Muscular insertions' in p['groups'] for p in parts)
            row.update(classification='source-authored-nonvolume-muscular-origin-attachment-site', occupied_volume_obligation_not_applicable=True, full_muscle_not_represented_by_this_site=True, actual_saved_surface_checks=[r for r in context['rows'] if r['part_id'] in {p['id'] for p in parts}])
        else:
            # This is an explicit conservative domain only, not a native solid,
            # anatomical cap, repaired source, or wall-thickness assertion.
            h = ConvexHull(v)
            hf = h.simplices.copy()
            tri = v[hf]
            normals = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
            inward = np.einsum('ij,ij->i', normals, h.equations[:,:3]) < 0
            hf[inward] = hf[inward][:, [0,2,1]]
            hull_mesh = Mesh(v, hf, wall.frame)
            hull_solid = helper.solid(hull_mesh, origin)
            assert hull_solid.status() == md.Error.NoError
            domain_checks = []
            for branch, solid in occupied.items():
                intersection = solid ^ hull_solid
                assert intersection.status() == md.Error.NoError
                domain_checks.append(dict(branch=branch, conservative_domain_intersection_empty=bool(intersection.is_empty()), conservative_domain_overlap_mm3=float(intersection.volume()), an_overlap_is_not_an_actual_collision=True))
            row.update(classification='unresolved-open-or-nonmanifold-physical-source' if not complete_closed else 'literal-closed-source', conservative_domain='convex hull of every literal whole-source point, not anatomy', conservative_domain_checks=domain_checks, source_hull_geometry_sha256=hull_mesh.geometry_sha, all_branches_conservatively_excluded=all(r['conservative_domain_intersection_empty'] for r in domain_checks))
        results.append(row)
        atomic_json(out / 'progress.json', dict(rows=results, elapsed_seconds=time.monotonic()-began), overwrite=True)
    evidence_paths = [work / p for p in ('default-runtime-context-v29/held-default-outcome-v30.json','classifier-normalization-v26/evidence.json','conservative-envelopes-v21/evidence.json','registered-context-v25/evidence.json')]
    result = dict(schema='remaining-facial-native-physical-source-audit/v31', source_catalogue=str(source.relative_to(ROOT)), source_catalogue_sha256=expected, registry_revision=revision, candidate_geometry_sha256=wall.geometry_sha, frame=wall.frame, source_evidence=[dict(path=str(p.relative_to(ROOT)),sha256=digest_file(p)) for p in evidence_paths], conduit_inputs=field_bindings, rows=results, source_native_shapes_unchanged=True, native_openings_not_filled=True, no_build_pack_install=True, no_anatomical_acceptance=True, elapsed_seconds=time.monotonic()-began, complete=True)
    assert digest_file(ROOT / binding['current_catalogue']) == expected
    assert Mesh.load(work / 'finite-wall.npz').geometry_sha == wall.geometry_sha
    atomic_json(out / 'evidence.json', result)
    print(json.dumps(dict(output=str((out/'evidence.json').relative_to(ROOT)), elapsed_seconds=result['elapsed_seconds'], rows=[dict(source_object_id=r['source_object_id'],classification=r['classification'],complete_literal_closed_classifier=r['complete_literal_closed_classifier'],conservative_domain_checks=r.get('conservative_domain_checks')) for r in results])))


if __name__ == '__main__':
    main()
