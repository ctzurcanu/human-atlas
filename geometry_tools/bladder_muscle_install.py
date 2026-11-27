"""Prepare the reviewed, source-bound male bladder muscle revision locally."""
from pathlib import Path
import copy, gzip, json
import numpy as np

from .common import ROOT, atomic_bytes, atomic_json, digest_file, identity
from .model import AtlasReader, Mesh
from .constraints import list_constraints


def prepare(db, delivery_directory, verification_directory, output):
    folder, verification_directory, output = map(lambda p: Path(p).resolve(), (delivery_directory, verification_directory, output))
    if output.exists():
        raise FileExistsError(output)
    delivery = json.loads((folder/'delivery.json').read_text())
    verified = json.loads((verification_directory/'verification.json').read_text())
    if not delivery['passed'] or not verified['passed'] or verified['delivery_sha256'] != digest_file(folder/'delivery.json'):
        raise ValueError('The actual emitted tissue and interface checks must pass')
    source = Path(delivery['source_catalogue'])
    if digest_file(source) != delivery['source_catalogue_sha256']:
        raise ValueError('The installed source assembly changed')
    before = json.loads(source.read_text()); atlas = copy.deepcopy(before)
    original = {p['id']: p for p in before['parts']}
    parts = {p['id']: p for p in atlas['parts']}
    marker = atlas['finiteUrinaryIntegration']
    source_reader = AtlasReader()
    new_ids = verified['new_parts']
    posterior = 'BLADDER:male:neck:vesicoprostatic'
    template = original['BLADDER:male:neck:neck-external-longitudinal']
    for chart, group in delivery['groups'].items():
        bank = Path(group['bank'])
        if digest_file(bank) != group['bank_sha256']:
            raise ValueError('The verified output bank changed')
        compressed = bank.with_suffix('.bin.gz')
        if not compressed.exists():
            atomic_bytes(compressed, gzip.compress(bank.read_bytes(), mtime=0))
        elif gzip.decompress(compressed.read_bytes()) != bank.read_bytes():
            raise ValueError('A compressed delivery bank differs')
        chunk = len(atlas['chunks'])
        atlas['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')), bytes=bank.stat().st_size,
                                    sha256=digest_file(bank), gzip='/local-models/'+str(compressed.relative_to(ROOT/'.local-models')), gzipBytes=compressed.stat().st_size))
        for pid, packed in group['parts'].items():
            if pid in parts:
                record = parts[pid]
            else:
                record = copy.deepcopy(template)
                role = pid.rsplit(':', 1)[1]
                record.update(id=pid, conceptId=pid, sourceId='constructed-bladder-muscle:'+role,
                              name='Vesicoprostaticus' if pid == posterior else 'Pubovesicalis ('+('left' if role.endswith('left') else 'right')+')',
                              ta98Term='A08.3.01.021M' if pid == posterior else 'A08.3.01.020', ta98Kind='exact',
                              urinaryRegion='neck' if pid == posterior else 'attachment', tissueWallLayer=role,
                              material='URINARY:male:bladder:'+role,
                              provenance=dict(label='Source-bound local bladder muscle construction',
                                              detail='Finite tissues with reciprocal source wall contacts; native skeleton and clinical urinary surfaces retained. Illustrative dimensions; anatomical review remains ongoing.'),
                              description='Posterior bladder-neck muscle connecting to the prostatic base.' if pid == posterior else 'Paired anterior bladder support muscle extending from the outer longitudinal neck muscle to the posterior pubic surface.')
                atlas['parts'].append(record); parts[pid] = record
                atlas['materials'][record['material']] = {**copy.deepcopy(atlas['materials'][template['material']]), 'geometryOrigin': packed['geometryOrigin']}
                atlas['concepts'].append(dict(id=pid, name=record['name'], elements=[pid], ta98Term=record['ta98Term'], ta98Kind='exact',
                                              partOf=['URINARY:male:bladder:A08.3.01.010']))
            record.update(copy.deepcopy(packed)); record['chunk'] = chunk
            if pid != 'ZA:Prostate':
                record['coverageLimitation'] = 'Finite source-bound wall and muscle attachments. Coat depths, regional contours, mucosal folds and complete anatomical acceptance remain under review.'
        existing_chart = next((c for c in marker['charts'] if c['name'] == chart), None)
        if existing_chart:
            existing_chart.update(partIds=list(group['parts']), bankSha256=group['bank_sha256'])
        else:
            marker['charts'].append(dict(name=chart, partIds=list(group['parts']), bankSha256=group['bank_sha256']))
    marker['partIds'] = [pid for chart in marker['charts'] for pid in chart['partIds']]
    muscle_terms = {'A08.3.01.020': [p for p in new_ids if ':pubovesical-' in p], 'A08.3.01.021M': [posterior]}
    for term, ids in muscle_terms.items():
        cid = 'URINARY:male:bladder:'+term
        atlas['concepts'].append(dict(id=cid, name='Pubovesicalis' if term.endswith('020') else 'Vesicoprostaticus',
                                     elements=ids, ta98Term=term, ta98Kind='exact', partOf=['URINARY:male:bladder:A08.3.01.010']))
        marker['termGroups'].append(dict(term=term, conceptId=cid, partIds=ids))
    concepts = {c['id']: c for c in atlas['concepts']}
    for group in marker['termGroups']:
        addition = new_ids if group['term'] in ['A08.3.01.001', 'A08.3.01.010'] else [posterior] if group['term'] == 'A08.3.01.005' else []
        if addition:
            group['partIds'] = [*group['partIds'], *addition]
            concepts[group['conceptId']]['elements'] = list(group['partIds'])
    parts['ZA:Urinary bladder']['replacementParts'].extend(new_ids)
    concepts['ZA:Urinary bladder']['elements'].extend(new_ids)
    # Clinical curves stay literal; only their source-triangle ordinals change.
    decoded = {}
    for chart, group in delivery['groups'].items():
        for pid in group['parts']:
            decoded[pid] = Mesh.load(folder/chart/(pid.replace(':', '_')+'.npz'))
    lookup = {pid: {tuple(sorted(map(tuple, t))): i for i, t in enumerate(mesh.triangles)} for pid, mesh in decoded.items()}
    for group in marker['termGroups']:
        for boundary in group.get('openingBoundaries', []):
            for i, pid in enumerate(boundary['samplePartIds']):
                old = source_reader.read_part(before, original[pid], frame=delivery['frame'])
                triangle = old.triangles[boundary['faceIds'][i]]
                face = lookup[pid][tuple(sorted(map(tuple, triangle)))]
                actual = decoded[pid].triangles[face]
                weights = np.asarray(boundary['barycentric'][i])
                boundary['barycentric'][i] = [float(weights[np.flatnonzero(np.all(triangle == point, axis=1))[0]]) for point in actual]
                boundary['faceIds'][i] = face
    replaced = set(delivery['cells'])
    marker['externalInterfaces'] = [r for r in marker['externalInterfaces'] if r['first'] not in replaced and r['second'] not in replaced]+verified['external_interfaces']
    source_contacts = []
    for pid in muscle_terms['A08.3.01.020']:
        check = next(i for i in verified['interface_checks'] if i['first'] == pid and i['second'].startswith('ZA:Hip bone.'))
        faces = check['parameters']['first_faces']; peer_faces = check['parameters']['second_faces']
        mesh = decoded[pid]
        words = mesh.attributes['normals'].astype(int)*mesh.attributes['normalSigns'][:, None]
        field = [(t.tolist(), words[mesh.faces[f]].tolist()) for f, t in zip(faces, mesh.triangles[faces])]
        source_contacts.append(dict(partId=pid, chart=delivery['chart_for'][pid], sourceContact=dict(first=pid, second=check['second'],
                                   firstFaces=faces, secondFaces=peer_faces, triangles=1, representation='fixed-native-surface-contact',
                                   coordinate_normal_signature=identity(field), actual_decoded_positions_and_opposed_words_match=True)))
    marker['muscleExtensions'] = dict(sourcePlanSha256=delivery['arrangement_sha256'], posteriorPartId=posterior, pubicAttachments=source_contacts)
    marker['sourceDeliverySha256'] = digest_file(folder/'delivery.json')
    # Face-order changes are explicit source-bound replacements of old seam
    # parameters. The live database and all original constraints stay literal.
    target_meshes = {db.mesh(pid, 'local-ta98')['id']: pid for pid in replaced if pid in original}
    rebindings = []
    for constraint in list_constraints(db, 'local-ta98', list(target_meshes)):
        if constraint['kind'] != 'shared-interface' or constraint['parameters'].get('verification_pending'):
            continue
        first, second = [db.mesh(constraint[k])['part_id'] for k in ['first_mesh_id', 'second_mesh_id']]
        parameters = dict(constraint['parameters'])
        for role, pid, peer in [('first', first, second), ('second', second, first)]:
            parameters[role+'_faces'] = [i for g in parts[pid]['tissueDrawGroups'] if g.get('neighbor') == peer
                                       for i in range(g['start']//3, (g['start']+g['count'])//3)]
        rebindings.append(dict(id=constraint['id'], first=first, second=second, expected_parameters_sha256=identity(constraint['parameters']),
                               parameters=parameters, evidence=dict(source=str(verification_directory/'verification.json'),
                               source_catalogue_sha256=digest_file(source), operation='Source-conserving tissue partition and explicit reciprocal triangle-index rebinding')))
    atlas['triangles'] = sum(p.get('triangles', p['indexCount']//3) for p in atlas['parts'])
    changed = [p['id'] for p in atlas['parts'] if p != original.get(p['id'])]
    output.mkdir(parents=True)
    atomic_json(output/'atlas.json', atlas)
    job = dict(kind='catalogue-candidate', parameters=dict(model='local-ta98', proposed_catalogue=str(output/'atlas.json'),
               expected_catalogue_sha256=digest_file(source), changed_parts=changed,
               inventory=str(ROOT/'.local-models/geometry-workbench/bladder-completion/source-v5/male-inventory.json'),
               assembly=dict(whole_organ='Male bladder with finite pubovesical and vesicoprostatic muscles',
                             adjoining_interfaces='Source-bound outer muscle, adventitia, fixed pubic facets and prostatic base; lining, ostia, renal roots and skeleton fixed'),
               shared_interface_rebindings=rebindings, interface_checks=verified['interface_checks'], max_pairs=1000000,
               evidence=dict(delivery=str(folder/'delivery.json'), verification=str(verification_directory/'verification.json'),
                             arrangement=delivery['arrangement'], source_manifest=str(ROOT/'.local-models/geometry-workbench/bladder-completion/male-muscular-attachment-sources-v1/sources.json'))))
    atomic_json(output/'job.json', job)
    return dict(catalogue=str(output/'atlas.json'), changed_parts=changed, source_bound_seam_rebindings=len(rebindings), new_parts=new_ids)
