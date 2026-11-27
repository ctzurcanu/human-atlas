"""Presentation-only catalogues; authoritative installed sources remain literal."""
from pathlib import Path
import copy
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / '.local-models/organ-models/spinal-cord'
SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas'

# These eleven source objects form one illustrative transverse section. Their
# legacy concepts have no TA98 identity, so the authority inventory alone cannot
# remove them from the lightweight whole-body presentation.
MALE_SECTION_DETAILS = {
    'ZA:' + name for name in (
        'Anterior horn of spinal cord', 'Anterior spinocerebellar tract',
        'Cuneate fasciculus', 'Gracile fasciculus',
        'Lateral vestibulospinal tract', 'Posterior horn of spinal cord',
        'Posterior spinocerebellar tract', 'Rubrospinal tract',
        'Spinal reticular process', 'Tectospinal tract',
        'White matter of spinal cord',
    )
}

def mesh_count(parts):
    return len({p['id'] for p in parts if not p.get('suppressed') and p.get('indexCount', 0) > 0})

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, separators=(',', ':'), allow_nan=False) + '\n')
    temporary.replace(path)

def compact_catalogue(a):
    used = sorted({p['chunk'] for p in a['parts']})
    prior = a.get('sourceChunkIndices')
    assert prior is None or len(prior) == len(a['chunks'])
    a['sourceChunkCount'] = a.get('sourceChunkCount', len(a['chunks']))
    a['sourceChunkIndices'] = [prior[i] if prior is not None else i for i in used]
    mapping = {old: new for new, old in enumerate(used)}
    a['chunks'] = [a['chunks'][i] for i in used]
    for p in a['parts']:
        p['chunk'] = mapping[p['chunk']]
    a['triangles'] = sum(p['indexCount'] // 3 for p in a['parts'] if not p.get('suppressed'))
    return a

def main():
    inventory = json.loads((ROOT / 'reports/ta98-spinal-cord-inventory.json').read_text())
    male_detail_ids = {
        mesh['part_id']
        for term in inventory['terms']
        if term['authority']['id'].startswith('A14.1.02.')
        and term['authority']['id'] != 'A14.1.02.001'
        for mesh in term['source_meshes']['male']
    }
    source_path = SOURCE / 'female/atlas.json'
    source = json.loads(source_path.read_text())
    marker = source['finiteSpinalIntegration']
    owned = marker['partIds'] + marker['coveringPartIds']
    ids = set(owned + marker['referencePartIds'])
    by_id = {p['id']: p for p in source['parts']}
    first = by_id[owned[0]]
    bank_path = ROOT / '.local-models' / source['chunks'][first['chunk']]['url'].removeprefix('/local-models/')
    assert digest(bank_path) == marker['bankSha256']
    raw = np.memmap(bank_path, dtype=np.uint8, mode='r')
    layers = ('grey', 'white', 'pia')
    layer_ids = {layer: f'SPINALOVERVIEW:female:{layer}' for layer in layers}
    membership = {pid: by_id[pid]['tissueWallLayer'] for pid in owned}
    charts = {layer: {} for layer in layers}
    for pid in owned:
        p = by_id[pid]
        assert p['chunk'] == first['chunk'] and p['positions'] == first['positions'] and p['vertexCount'] == first['vertexCount']
        indices = np.ndarray((p['indexCount'],), dtype='<u4', buffer=raw, offset=p['indices'])
        for group in p['tissueDrawGroups']:
            neighbor = group.get('neighbor')
            if group.get('computationalOnly') or membership.get(neighbor) == membership[pid]:
                continue
            assert not neighbor or neighbor in membership
            key = membership.get(neighbor)
            charts[membership[pid]].setdefault(key, []).append(indices[group['start']:group['start'] + group['count']])
    arrays, groups, layer_ranges = [], {}, {}
    offset = 0
    for layer in layers:
        groups[layer] = []
        local = 0
        for neighbor, pieces in charts[layer].items():
            values = np.concatenate(pieces)
            arrays.append(values)
            groups[layer].append(dict(start=local, count=len(values), **({'neighbor': layer_ids[neighbor]} if neighbor else {})))
            local += len(values)
        layer_ranges[layer] = (offset, local)
        offset += local
    original_indices = np.concatenate(arrays)
    vertices, inverse = np.unique(original_indices, return_inverse=True)
    positions = np.ndarray((first['vertexCount'], 3), dtype='<f4', buffer=raw, offset=first['positions'])[vertices].copy()
    normals = np.ndarray((first['vertexCount'], 3), dtype='<i2', buffer=raw, offset=first['normals'])[vertices].copy()
    signs = np.ndarray((first['vertexCount'],), dtype='i1', buffer=raw, offset=first['normalSigns'])[vertices].copy()
    OUT.mkdir(parents=True, exist_ok=True)
    compact_bank = OUT / 'overview.bin'
    with compact_bank.open('wb') as stream:
        offsets = {}
        for name, values in [('positions', positions), ('normals', normals), ('normalSigns', signs), ('indices', inverse.astype('<u4'))]:
            stream.write(b'\0' * ((-stream.tell()) % 4))
            offsets[name] = stream.tell()
            stream.write(values.tobytes())
    overview_parts, overview_concepts = [], []
    terms = {'grey': 'A14.1.02.020', 'white': 'A14.1.02.201', 'pia': 'A14.1.01.309'}
    names = {'grey': 'Grey substance', 'white': 'White substance', 'pia': 'Spinal pia mater'}
    for layer in layers:
        p = copy.deepcopy(next(by_id[i] for i in owned if membership[i] == layer))
        begin, count = layer_ranges[layer]
        p.update(id=layer_ids[layer], conceptId=layer_ids[layer], name=names[layer], vertexCount=len(vertices), indexCount=count,
                 positions=offsets['positions'], normals=offsets['normals'], normalSigns=offsets['normalSigns'],
                 indices=offsets['indices'] + begin * 4, tissueDrawGroups=groups[layer], tissueWallAssembly='spinal-cord-overview:female',
                 groups=['Spinal cord', names[layer]], regions=[], chunk=len(source['chunks']))
        p.pop('spinalCellRegion', None)
        p['coverageLimitation'] = 'Whole-body overview consolidated from the detailed cord. Nuclear, tract and segment controls are in the separate detailed model.'
        used = inverse[begin:begin + count]
        world = positions[used].astype(np.float64) + np.asarray(p['geometryOrigin'])
        p['bounds'] = [world.min(axis=0).tolist(), world.max(axis=0).tolist()]
        overview_parts.append(p)
        overview_concepts.append(dict(id=p['conceptId'], name=p['name'], elements=[p['id']], ta98Term=terms[layer], ta98Kind='exact', partOf=marker['wholeConceptId']))
        # Literal corner coordinates, packed normals and signs survive compaction.
        assert np.array_equal(vertices[used], original_indices[begin:begin + count])
    # Each cross-coat chart remains reciprocal; same-coat charts have been removed.
    for p in overview_parts:
        for g in p['tissueDrawGroups']:
            if 'neighbor' in g:
                neighbor = next(q for q in overview_parts if q['id'] == g['neighbor'])
                assert sum(h['count'] for h in neighbor['tissueDrawGroups'] if h.get('neighbor') == p['id']) == g['count']
    detail = {key: copy.deepcopy(source[key]) for key in ('version', 'sex', 'source', 'chunks', 'materials', 'concepts', 'finiteSpinalIntegration') if key in source}
    detail.update(version='Spinal cord · detailed', scope='organ', parts=[copy.deepcopy(by_id[i]) for i in owned + marker['referencePartIds']],
                  defaultHidden=list(marker['referencePartIds']))
    detail['concepts'] = [{**c, 'elements': [i for i in c['elements'] if i in ids]} for c in detail['concepts'] if any(i in ids for i in c['elements'])]
    detail.pop('bowelWallPrototype', None)
    detail.pop('delivery', None)
    compact_catalogue(detail)
    write_json(OUT / 'atlas.json', detail)
    detail_meshes = mesh_count(detail['parts'])
    results = {}
    for sex in ('female', 'male'):
        original = json.loads((SOURCE / sex / 'atlas.json').read_text())
        body = copy.deepcopy(original)
        whole = marker['wholeConceptId'] if sex == 'female' else 'LOCAL:male:spinal_cord'
        if sex == 'female':
            removed = {p['id'] for p in body['parts'] if p['id'].startswith('TA98CORD:') or p['id'] == marker['sourcePartId']}
            body['parts'] = [p for p in body['parts'] if p['id'] not in removed] + copy.deepcopy(overview_parts)
            body['concepts'] = [{**c, 'elements': [i for i in c['elements'] if i not in removed]} for c in body['concepts'] if c['id'] != whole and any(i not in removed for i in c['elements'])]
            body['concepts'] += copy.deepcopy(overview_concepts) + [dict(id=whole, name='Spinal cord', ta98Term='A14.1.02.001', ta98Kind='exact', elements=[layer_ids['grey'], layer_ids['white']])]
            body['defaultHidden'] = [i for i in body.get('defaultHidden', []) if i not in removed]
            body['chunks'].append(dict(url='/local-models/organ-models/spinal-cord/overview.bin', bytes=compact_bank.stat().st_size))
            body.pop('finiteSpinalIntegration', None)
            body['organOverviews'] = [dict(wholeConceptId=whole, partIds=list(layer_ids.values()), systemPath=['Central nervous system', 'Torso & pelvis', 'Spinal cord'])]
            overview_ids = list(layer_ids.values())
        else:
            overview_ids = next(c['elements'] for c in body['concepts'] if c['id'] == whole)
            for p in body['parts']:
                if p['id'] in MALE_SECTION_DETAILS:
                    assert p['system'] == 'nervous' and p['conceptId'] == p['id']
                    assert 1.442 < p['bounds'][0][1] < p['bounds'][1][1] < 1.449
            # Small donor cross-section objects are details, not extra copies
            # of the cord inside a whole-body demonstration.
            detailed_concepts = {c['id'] for c in body['concepts'] if c.get('ta98Term', '').startswith('A14.1.02.') and c['id'] != whole}
            # Legacy donor concepts often have no embedded TA98 term. Their
            # authoritative inventory records still identify these cord details.
            removed = {p['id'] for p in body['parts'] if (p['conceptId'] in detailed_concepts or p['id'] in male_detail_ids or p['id'] in MALE_SECTION_DETAILS) and p['id'] not in overview_ids}
            body['parts'] = [p for p in body['parts'] if p['id'] not in removed]
            body['concepts'] = [{**c, 'elements': [i for i in c['elements'] if i not in removed]} for c in body['concepts'] if any(i not in removed for i in c['elements'])]
            body['defaultHidden'] = [i for i in body.get('defaultHidden', []) if i not in removed and i not in overview_ids]
            next(c for c in body['concepts'] if c['id'] == whole).update(ta98Term='A14.1.02.001', ta98Kind='exact')
            body['organOverviews'] = [dict(wholeConceptId=whole, partIds=overview_ids, systemPath=['Central nervous system', 'Torso & pelvis'])]
        body['separateOrganModels'] = [dict(id='spinal-cord-detail', name='Spinal cord · detailed', model='local-spinal-cord', meshCount=detail_meshes, parentConceptId=whole, parentName='Spinal cord', parentTerm='A14.1.02.001', overviewPartIds=overview_ids)]
        assert not set(p['id'] for p in detail['parts']) & set(p['id'] for p in body['parts'])
        compact_catalogue(body)
        write_json(OUT / f'{sex}-body.json', body)
        results[sex] = dict(sourceSha256=digest(SOURCE / sex / 'atlas.json'), sourceMeshes=mesh_count(original['parts']), bodyMeshes=mesh_count(body['parts']), overviewMeshes=len(overview_ids), removedDetailMeshes=mesh_count([p for p in original['parts'] if p['id'] in removed]))
    write_json(OUT / 'evidence.json', dict(schema='separate-organ-presentation/v1', sourceInventorySha256=digest(ROOT / 'reports/ta98-spinal-cord-inventory.json'), sourceBankSha256=marker['bankSha256'], overviewBankSha256=digest(compact_bank), overviewBytes=compact_bank.stat().st_size, detailMeshes=detail_meshes, detailsShareOriginalBank=True, sourceCoordinatesAndNormalWordsPreserved=True, reciprocalCoatCharts=True, models=results))
    print(json.dumps(dict(output=str(OUT), overviewBytes=compact_bank.stat().st_size, detailMeshes=detail_meshes, models=results)))

if __name__ == '__main__':
    main()
