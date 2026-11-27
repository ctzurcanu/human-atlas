"""Bind the complete spinal authority inventory before constructing tissue.

This stage reads the installed atlases, not old source-registration reports.
No geometry, visibility, camera, hierarchy or live catalogue is changed.
"""
from pathlib import Path
import argparse, json, sqlite3, time
from geometry_tools.common import ROOT, atomic_json, digest_file

REFERENCE = 'https://nba.uth.tmc.edu/neuroscience/m/s2/chapter03.html'
LAB = 'https://nba.uth.tmc.edu/NEUROANATOMY/L2/Lab02p11_index.html'
ap = argparse.ArgumentParser()
ap.add_argument('--output', required=True)
args = ap.parse_args()
started = time.monotonic()
out = Path(args.output).resolve()
out.mkdir(parents=True, exist_ok=False)
db = sqlite3.connect(ROOT / '.local-models/ta98-authority.sqlite')
db.row_factory = sqlite3.Row
authority = [dict(r) for r in db.execute('select * from ta98 order by id')]
metadata = json.loads((ROOT / 'app/data/ta98-metadata.json').read_text())
paths = {s: ROOT / f'.local-models/ta98-review/bowel-wall-local-atlas/{s}/atlas.json' for s in ['female', 'male']}
atlases = {s: json.loads(p.read_text()) for s, p in paths.items()}
coverings = {'A14.1.01.' + s for s in ['001','002','003','101','111','112','201','202','203','206','218','219','301','309','310','311','401','402','403']}
roots = {'A14.2.00.' + s for s in ['027','028','029','030','031','032','033','034','035','036','037']}
terms = [r for r in authority if r['id'].startswith('A14.1.02.') or r['id'] in coverings | roots]

def classification(row):
    code, name = row['id'], row['name_en'].lower()
    if row['immaterial'] or any(s in name for s in ['canal', 'ventricle', 'cistern', ' space']):
        return 'space-or-lumen; boundary-owned void, never a solid plug'
    if any(s in name for s in ['fissure','sulcus']):
        return 'boundary curve and actual adjoining surface relief'
    if 'nucleus' in name or 'nuclei' in name or 'lamina' in name or code.endswith('.302'):
        return 'cytoarchitectural tissue region; level-restricted 3D territory, not a point'
    if any(s in name for s in ['tract','fibres','fascicul','funicul']):
        return 'longitudinal fibre territory; regional extent and transitions required'
    if any(s in name for s in ['dura','arachnoid','pia','meninx','meninges']):
        return 'covering tissue; finite layer with explicitly owned inner/outer interfaces'
    if code in {'A14.1.02.014','A14.1.02.015','A14.1.02.016','A14.1.02.017','A14.1.02.018'}:
        return 'segment-defined whole-cord region; nerve-root levels, not vertebral levels'
    if any(s in name for s in ['root','ramus','rami','nerve','plexus','cauda']):
        return 'nerve tissue / root interface; native branch ownership preserved'
    return 'finite 3D tissue or assembly of its explicitly owned tissue cells'

def source_records(sex, code):
    atlas = atlases[sex]
    concept_terms = {c['id']: (c.get('ta98Term'), c.get('ta98Kind')) for c in atlas['concepts']}
    rows = []
    for p in atlas['parts']:
        key = p['conceptId']
        match = metadata.get('byConceptMatch', {}).get(key, {})
        reviewed = metadata.get('reviewedCatalogueIdentities', {}).get(key, {})
        term = p.get('ta98Term') or reviewed.get('term') or concept_terms.get(key, (None, None))[0] or match.get('term')
        if term != code:
            continue
        rows.append(dict(part_id=p['id'], concept_id=key, name=p['name'],
            hidden=p.get('suppressed',False) or p['id'] in atlas.get('defaultHidden',[]) or p.get('defaultHidden',False),
            triangles=p['indexCount']//3, source_id=p.get('sourceId'), bounds=p['bounds'],
            bank_url=atlas['chunks'][p['chunk']]['url'], record=p,
            limitation='Source identity is not proof of regional extent, finite tissue or anatomical acceptance.'))
    return rows

rows = []
for term in terms:
    sources = {s: source_records(s, term['id']) for s in paths}
    rows.append(dict(authority=term, representation_class=classification(term), source_meshes=sources,
        missing_representation={s: not any(not p['hidden'] for p in sources[s]) for s in sources},
        acceptance='unverified; current source or planned construction'))

plan = dict(schema='spinal-cord-assembly-plan/v1',
    authority=dict(path='.local-models/ta98-authority.sqlite',sha256=digest_file(ROOT/'.local-models/ta98-authority.sqlite')),
    source_catalogues={s:dict(path=str(p),sha256=digest_file(p)) for s,p in paths.items()},
    terms=rows,cord_term_count=sum(r['authority']['id'].startswith('A14.1.02.') for r in rows),
    references=[dict(url=REFERENCE,scope='Macroscopic shape, regional grey/white variation, roots and meninges'),
        dict(url=LAB,scope='Thoracic horns, Clarke column and above-T6 cuneate extent; illustrations are consulted, not copied')],
    assembly=dict(primary='Native female spinal cord, native dura and cauda equina, native ipsilateral root sources.',
        transfer='Only after female surface and finite tissue review: transport the common tissue charts into the native male cord; keep male skeleton fixed.',
        fixed='Every vertebra, skull, muscle, vessel and existing native peripheral nerve outside its declared root attachment.',
        source_frame='Metres in each installed assembled model, decoded by AtlasReader exactly once.',
        topology_repair='At the conus tip only, collapse the five collinear terminal controls to the existing central tip; maximum 0.321 mm. Remove only resulting zero-area faces, preserve all other source controls and attributes.',
        whole_surface_gate='One closed genus-zero cord, no proper self-crossings, fixed-context preservation and actual surface / transverse-section inspection before tissue partition.',
        planned_parts=['Native outer-cord envelope','Cervical/thoracic/lumbar/sacral/coccygeal regions located from native root attachments',
            'Variable butterfly-shaped grey matter and surrounding white matter with reciprocal interfaces',
            'Anterior/posterior horns and intermediate/commissural grey tissue',
            'Anterior/lateral/posterior white funiculi, with level-specific tract extents',
            'Pia, arachnoid and dura as separate coverings','Anterior/posterior roots and cauda equina'],
        interfaces=['Grey/white: one geometric interface with opposed ownership; no coincident separate envelopes.',
            'Central canal: void bounded by tissue; no solid cylinder representing the lumen.',
            'Pia follows the cord; subarachnoid interval remains a space; dura/root sleeves retain native exits.',
            'Rootlets join anterolateral/posterolateral cord entry zones; root names do not assign vertebral levels.',
            'Cord continues into medulla cranially and tapers at conus; filum and cauda are distinct.',
            'Unverified microscopic nuclei/tract territories remain missing until their correct region and extent are constructed.']),
    elapsed_seconds=time.monotonic()-started,geometry_changed=False,installed=False,anatomical_acceptance=False)
atomic_json(out/'inventory.json',plan)
report = ROOT/'reports/ta98-spinal-cord-inventory.json'
atomic_json(report,plan)
lines=['# Spinal cord: authority inventory and assembly', '',
    'Female first; transfer the reviewed common assembly to the male native frame. Skeleton fixed.', '',
    f'{plan["cord_term_count"]} cord terms; {len(rows)-plan["cord_term_count"]} applicable covering/root terms. A source match is not anatomical acceptance.', '',
    'Source plan: reports/ta98-spinal-cord-inventory.json. Local stage: '+str(out.relative_to(ROOT))+'.', '',
    'An opening is a boundary-owned void. Sulci are curves and surface relief. Nuclei and tracts require finite regional tissue; a point or generic outer envelope is insufficient.', '',
    '| TA98 | Parent | Name | Representation | Female active sources | Male active sources |',
    '| --- | --- | --- | --- | --- | --- |']
for row in rows:
    t=row['authority']
    src=lambda s:', '.join(p['part_id'] for p in row['source_meshes'][s] if not p['hidden']) or '**Missing**'
    lines.append('| '+' | '.join([t['id'],t['parent_id'] or '',t['name_en'],row['representation_class'],src('female'),src('male')])+' |')
lines += ['', 'Whole-surface repair and all adjoining interfaces are defined in the JSON plan before construction.', '',
    'References: [UTHealth spinal-cord anatomy]('+REFERENCE+'), [UTHealth thoracic sections]('+LAB+').']
(ROOT/'TA98-SPINAL-CORD.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(cord_terms=plan['cord_term_count'],covering_and_root_terms=len(rows)-plan['cord_term_count'],
    female_active_term_sources=sum(not r['missing_representation']['female'] for r in rows),
    male_active_term_sources=sum(not r['missing_representation']['male'] for r in rows),
    elapsed_seconds=time.monotonic()-started,geometry_changed=False)))
