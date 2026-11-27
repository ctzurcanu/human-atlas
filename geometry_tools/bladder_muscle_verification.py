"""Independent replay of emitted bladder muscle cells and their source peers."""
from pathlib import Path
import json
import numpy as np

from .common import atomic_json, digest_file, identity
from .constraints import check_constraint
from .model import AtlasReader, Mesh
from .contact import pair_contacts


def verify(delivery_directory, output):
    folder, output = Path(delivery_directory).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    spec = json.loads((folder/'delivery.json').read_text())
    if not spec['passed'] or digest_file(spec['source_catalogue']) != spec['source_catalogue_sha256']:
        raise ValueError('Passed output and unchanged source catalogue are required')
    atlas = json.loads(Path(spec['source_catalogue']).read_text())
    records = {p['id']: p for p in atlas['parts']}
    reader, actual, parts = AtlasReader(), {}, {}
    for chart, group in spec['groups'].items():
        if digest_file(group['bank']) != group['bank_sha256']:
            raise ValueError('An emitted geometry bank changed')
        for pid, part in group['parts'].items():
            actual[pid] = Mesh.load(folder/chart/(pid.replace(':', '_')+'.npz'))
            parts[pid] = part
    def get(pid):
        if pid not in actual:
            actual[pid] = reader.read_part(atlas, records[pid], frame=spec['frame'])
        return actual[pid]
    interfaces, checks, ledger, seen = [], [], [], set()
    def face_indices(pid, peer):
        return [i for draw in parts.get(pid, records.get(pid, {})).get('tissueDrawGroups', []) if draw.get('neighbor') == peer
                for i in range(draw['start']//3, (draw['start']+draw['count'])//3)]
    for pid in spec['cells']:
        for use in spec['cells'][pid]:
            peer = use.get('neighbor')
            if not peer or tuple(sorted([pid, peer])) in seen:
                continue
            seen.add(tuple(sorted([pid, peer])))
            a, b = get(pid), get(peer)
            af = face_indices(pid, peer)
            bf = face_indices(peer, pid)
            reference = not bf
            if reference:
                index = {tuple(sorted(map(tuple, t))): i for i, t in enumerate(b.triangles)}
                bf = [index.get(tuple(sorted(map(tuple, t))), -1) for t in a.triangles[af]]
                if not af or -1 in bf:
                    raise ValueError('A source contact is not a literal native facet: '+peer)
            constraint = dict(kind='shared-interface', first=pid, second=peer,
                              parameters=dict(first_faces=af, second_faces=bf),
                              evidence=dict(source=str(folder/'delivery.json')))
            check = check_constraint(constraint, a, b)
            checks.append(dict(first=pid, second=peer, source_surface_contact=reference, **check))
            interfaces.append(constraint)
            if spec['chart_for'].get(pid) != spec['chart_for'].get(peer) and not reference:
                normals = a.attributes['normals'].astype(int)*a.attributes.get('normalSigns', np.ones(len(a.vertices), int))[:, None]
                field = {tuple(sorted(map(tuple, t))): {tuple(q): tuple(map(int, w)) for q, w in zip(t, n)}
                         for t, n in zip(a.triangles[af], normals[a.faces[af]])}
                ledger.append(dict(first=pid, second=peer, triangles=len(af),
                                   coordinate_normal_signature=identity([(k, sorted(v.items())) for k, v in sorted(field.items())]),
                                   actual_decoded_positions_and_opposed_words_match=check['passed']))
    old_fixed, altered = [], []
    modified_surfaces = {'BLADDER:male:neck:adventitia', 'BLADDER:male:neck:neck-external-longitudinal', 'ZA:Prostate'}
    for pid in spec['cells']:
        if pid not in records:
            continue
        old = reader.read_part(atlas, records[pid], frame=spec['frame'])
        if pid in modified_surfaces:
            altered.append(pid)
            continue
        match = old.surface_sha == get(pid).surface_sha
        old_fixed.append(dict(part_id=pid, literal_oriented_surface=match))
    new_ids = [pid for pid in spec['cells'] if pid not in records]
    hidden = set(atlas.get('defaultHidden', []))
    packing, failures = [], []
    for pid in new_ids:
        a = get(pid)
        low, high = a.bounds
        for record in atlas['parts']:
            peer = record['id']
            if record.get('suppressed') or peer in hidden or peer == pid:
                continue
            bound = np.asarray(record['bounds'], float)
            if np.any(bound[0] > high) or np.any(bound[1] < low):
                continue
            b = get(peer)
            if np.any(b.bounds[0] > high) or np.any(b.bounds[1] < low):
                continue
            arrays, result = pair_contacts(a, b, 0., max_pairs=1000000)
            declared = face_indices(pid, peer)
            bad = result['crossing_face_pairs']
            overlaps = arrays['pairs'][arrays['relation'] == 'coplanar-overlap']
            undeclared_overlap = int(sum(int(face) not in set(declared) for face in overlaps[:, 0]))
            passed = not result['truncated'] and bad == 0 and undeclared_overlap == 0
            packing.append(dict(first=pid, second=peer, crossing_face_pairs=bad,
                                coplanar_pairs=len(overlaps), undeclared_coplanar_pairs=undeclared_overlap,
                                complete=not result['truncated'], passed=passed))
            if not passed:
                failures.append(dict(first=pid, second=peer, crossing_face_pairs=bad, undeclared_overlap=undeclared_overlap))
        print('packing', pid, 'cumulative failures', len(failures), flush=True)
    passed = all(c['passed'] for c in checks) and all(c['literal_oriented_surface'] for c in old_fixed) and not failures
    output.mkdir(parents=True)
    report = dict(delivery_sha256=digest_file(folder/'delivery.json'), passed=passed, interface_checks=interfaces,
                  decoded_seams=checks, external_interfaces=ledger, exact_old_surfaces=old_fixed,
                  changed_source_surfaces=altered, new_parts=new_ids, packing=packing, failures=failures,
                  complete_active_body_bounds_shortlist=True, hidden_review_alternatives_excluded=True,
                  source_lining_and_renal_roots_fixed=all(c['literal_oriented_surface'] for c in old_fixed),
                  anatomical_acceptance=False, installed=False)
    atomic_json(output/'verification.json', report)
    return report
