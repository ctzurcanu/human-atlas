"""Anatomical dimension and construction obligations for the 71 cardiac sites.

Collections select typed children. Empty spaces and orifices are never tissue
solids. Variant adult foramen ovale is represented by its closed septal site.
"""

CURVES = {
    'A12.1.00.006', 'A12.1.00.008', 'A12.1.00.009', 'A12.1.00.010',
    'A12.1.00.011', 'A12.1.01.006', 'A12.1.01.013',
}
SURFACES = {
    'A12.1.00.030', 'A12.1.01.005', 'A12.1.01.007', 'A12.1.01.012',
    'A12.1.02.015', 'A12.1.04.017',
}
OPENINGS = {
    'A12.1.01.004', 'A12.1.01.009', 'A12.1.01.010', 'A12.1.01.011',
    'A12.1.02.002', 'A12.1.02.009', 'A12.1.03.004', 'A12.1.04.002',
    'A12.1.04.010',
}
POINTS = {'A12.1.02.016', 'A12.1.04.018'}
SPACES = {'A12.1.04.009', 'A12.1.07.001', 'A12.1.07.002'}
COLLECTIONS = {
    'A12.1.00.018': ['A12.1.01.002', 'A12.1.03.002'],
    'A12.1.00.020': ['A12.1.02.021', 'A12.1.04.011'],
    'A12.1.00.022': ['A12.1.02.017', 'A12.1.02.018', 'A12.1.02.019',
                      'A12.1.04.007', 'A12.1.04.008'],
    'A12.1.06.002': ['A12.1.06.003', 'A12.1.06.004', 'A12.1.06.005',
                      'A12.1.06.006', 'A12.1.06.007', 'A12.1.06.008'],
    'A12.1.08.001': ['A12.1.08.002', 'A12.1.08.003', 'A12.1.08.004',
                      'A12.1.08.005', 'A12.1.08.010'],
    'A12.1.08.005': ['A12.1.08.006', 'A12.1.08.007', 'A12.1.08.008',
                      'A12.1.08.009'],
}
EXISTING_MESHES = {'A12.1.04.007', 'A12.1.04.008'}
SOURCE_CIDS = {
    'A12.1.00.025': 20644, 'A12.1.00.026': 23777,
    'A12.1.01.002': 20454, 'A12.1.01.008': 18054,
    'A12.1.02.008': 16459, 'A12.1.02.020': 19988,
    'A12.1.02.021': 18687, 'A12.1.03.002': 23809,
    'A12.1.04.011': 22800,
    'A12.1.06.003': 18229, 'A12.1.06.004': 17776,
    'A12.1.06.005': 22137, 'A12.1.06.006': 20044,
    'A12.1.06.007': 18445, 'A12.1.06.008': 23412,
    'A12.1.08.002': 23932, 'A12.1.08.006': 23251,
}

def required_type(term):
    if term in COLLECTIONS: return 'collection'
    if term in POINTS: return 'point-set'
    if term in CURVES: return 'boundary-curve'
    if term in SURFACES: return 'surface-region'
    if term in OPENINGS: return 'opening-rim'
    if term in SPACES: return 'space-region'
    return 'finite-mesh'

if __name__ == '__main__':
    from pathlib import Path
    import argparse, collections, hashlib, json
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', required=True)
    args = ap.parse_args()
    out = root / '.local-models/geometry-workbench/heart' / args.stage
    out.mkdir(exist_ok=False)
    bundle = json.loads((root / '.local-models/ta98-review/bowel-wall-local-atlas/male/cardiac-label-sites.json').read_text())
    inventory = json.loads((root / '.local-models/geometry-workbench/heart/male-cardiopulmonary-final-install-v525/inventory.json').read_text())
    authority = {r['authority']['id']: r['authority'] for r in inventory['terms']}
    rows = []
    for group in bundle['integration']['termGroups']:
        term = group['term']
        kind = required_type(term)
        rows.append(dict(term=term, name=group['name'], parent=authority[term]['parent_id'],
                         authority_materiality=authority[term]['type_of_entity'], variant=authority[term]['variant'],
                         required_type=kind, current_type=group['representation'],
                         source_carriers=group['partIds'], independent_reference_cid=SOURCE_CIDS.get(term),
                         existing_distinct_mesh=term in EXISTING_MESHES,
                         children=COLLECTIONS.get(term),
                         construction='source-indexed native facet footprint' if kind=='surface-region' else
                         'empty channel/rim; never a solid disk' if kind=='opening-rim' else
                         'bounded empty-space envelope; never filled tissue' if kind=='space-region' else
                         'closed finite geometry with actual scale and anatomical placement review' if kind=='finite-mesh' else
                         'ordered native surface path' if kind=='boundary-curve' else
                         'actual leaflet attachment junctions' if kind=='point-set' else 'typed child union'))
    assert len(rows)==71
    counts = dict(collections.Counter(r['required_type'] for r in rows))
    spec = dict(schema='cardiac-typed-reference-plan/v1', sex='male-first', terms=rows,
                counts=counts, anatomy_accepted=False, original_tissue_and_skeleton_fixed=True,
                secondary_transfer='literal shared geometry bank with current proper female similarity',
                references=['https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/12.1.08%20TA98%20EN.htm',
                            'https://ifaa.unifr.ch/Public/TNAEntryPage/auto/unit/LAEN/TAH3643%20Unit%20EN.htm',
                            'https://pubmed.ncbi.nlm.nih.gov/22109383/'])
    (out/'plan.json').write_text(json.dumps(spec,indent=2)+'\n')
    report = root/'TA98-CARDIAC-ASSEMBLY.md'
    heading='## Typed reference correction — 21 October 2026, '+args.stage
    text=report.read_text()
    assert heading not in text
    table='\n'.join('| '+r['term']+' | '+r['name']+' | '+r['parent']+' | '+r['required_type']+' | '+', '.join(r['source_carriers'])+' |' for r in rows)
    report.write_text(text+'\n'+heading+'\n\nThe 71 provisional controls below are audited before construction. The complete 113-term authority inventory above remains binding. No material structure is accepted as a point. Collections select their actual children; cavities and orifices remain empty. A default adult foramen-ovale site is closed. Source-space registration alone is not anatomical acceptance.\n\n| TA98 | Structure | Authority parent | Required representation | Current carrier(s) |\n|---|---|---|---|---|\n'+table+'\n\nCurrent skeleton, heart, valves, vessel roots and surrounding organs are fixed for this reference pass. New reference geometry requires rendered and section review; the rejected v639 serosal bank must not be installed.\n')
    print(json.dumps(dict(stage=str(out),counts=counts,terms=len(rows))),flush=True)
