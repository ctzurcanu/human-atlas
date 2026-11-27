"""Source-bound bladder inventory and assembly review, before tissue construction."""
from pathlib import Path
import json
import sqlite3

import numpy as np

from .audit import audit
from .common import ROOT, atomic_json, digest_file
from .constraints import closure
from .db import Database
from .model import AtlasReader, Mesh
from .registry import inventory, load_registered
from .review import prepare as prepare_review


ROOT_TERM = 'A08.3.01.001'
CATALOGUES = {sex: ROOT / f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json'
              for sex in ['male', 'female']}
MODELS = {'male': 'local-ta98', 'female': 'local-female-ta98'}
REFERENCES = [
    'https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/08.3.01%20TA98%20EN.htm',
    'https://www.histologyguide.org/slideview/MHS-214-bladder/16-slide-1.html',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC4566164/',
    'https://journals.biologists.com/dev/article/134/20/3763/64555/The-development-of-the-bladder-trigone-the-center',
]


def classifications():
    with sqlite3.connect(ROOT / '.local-models/ta98-authority.sqlite') as c:
        c.row_factory = sqlite3.Row
        terms = {r['id']: dict(r) for r in c.execute('SELECT * FROM ta98')}
    codes = {ROOT_TERM}
    while True:
        children = {k for k, t in terms.items() if t['parent_id'] in codes}
        if children <= codes:
            break
        codes |= children
    layers = {8, 9, 10, 12, 13, 17, 18, 19, 22, 23}
    result = {}
    for code in codes:
        if code == ROOT_TERM:
            kind = 'whole-organ'
        elif code.startswith('A09.') or code in ['A08.3.01.026', 'A08.3.01.027M', 'A08.3.01.028F']:
            kind = 'opening'
        elif code in ['A08.3.01.002', 'A08.3.01.025']:
            kind = 'boundary-curve'
        elif code[-3:].isdigit() and int(code[-3:]) in layers:
            kind = 'wall-layer'
        else:
            kind = 'tissue-region'
        result[code] = kind
    return result


def assembly(meshes, frame, sources):
    vertices, faces, start = [], [], 0
    for mesh in meshes:
        if mesh.frame != frame:
            raise ValueError('Do not combine different source frames')
        vertices.append(mesh.vertices); faces.append(mesh.faces + start)
        start += len(mesh.vertices)
    return Mesh(np.concatenate(vertices), np.concatenate(faces), frame,
                source=dict(kind='diagnostic-assembly', source_parts=sources,
                            material_attributes_retained_in_individual_exports=True))


def crop(mesh, bounds):
    """Uncapped review crop only; original exported mesh remains untouched."""
    low, high = np.asarray(bounds)
    t = mesh.triangles
    keep = np.all(t.max(1) >= low, axis=1) & np.all(t.min(1) <= high, axis=1)
    if not keep.any():
        return None
    return Mesh(mesh.vertices, mesh.faces[keep], mesh.frame, mesh.attributes,
                dict(mesh.source, diagnostic_uncapped_crop=True))


def prepare(output):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    db = Database(); db.assert_current_catalogues()
    classes = classifications(); atomic_json(output / 'classifications.json', classes)
    inventories = {}
    for sex, model in MODELS.items():
        inv = inventory(db, ROOT_TERM, model, classifications=classes)
        atomic_json(output / f'{sex}-inventory.json', inv); inventories[sex] = inv

    # This report is written before source cleanup or new surfaces are made.
    report = ['# Urinary bladder: full TA98 construction inventory', '',
              'Male is the construction master; transfer follows review of the complete male assembly.',
              'Current whole-organ source correspondences do not establish named finite layers.', '',
              '| TA98 | Parent | Term | Representation | Male source mappings | Female source mappings |',
              '| --- | --- | --- | --- | --- | --- |']
    for first, second in zip(inventories['male']['terms'], inventories['female']['terms']):
        t = first['authority']
        assert t['id'] == second['authority']['id']
        source = lambda item: '; '.join(s['part_id'] + (' (retired)' if s['suppressed'] else '')
                                      for s in item['source_meshes']) or 'Missing'
        report.append(f"| {t['id']} | {t['parent_id']} | {t['name_en']} | {first['representation_class']} | {source(first)} | {source(second)} |")
    report += ['', 'The apex is a point landmark and the interureteric crest is a ridge; neither is an independent filler solid. '
               'The trigone is a selected inner-wall region bounded by two ureteric ostia and the internal urethral ostium. '
               'Filling/voiding entries describe states of the same opening, not additional holes. '
               'Male/female sex-qualified branches remain distinct; the cached parent of female filling/voiding references requires comparison with the authority before mapping.', '',
               '## Assembly and adjoining interfaces', '',
               'One shared tissue bank must hold the whole wall, detrusor/trigonal and neck muscle subdivisions, connective coats and mucosa. '
               'Peritoneal serosa/subserosa covers the reviewed superior surface only; uncovered external tissue is adventitia. '
               'Mucosal folds taper out at the trigone. Named features must use distinct, anatomically located selections and ordinary labels.', '',
               'Required interfaces: bilateral oblique intramural ureters with patent ostia; shared urethral neck; '
               'median umbilical ligament attachment; pubovesical and sex-specific muscle attachments; '
               'prostate or anterior vagina, rectum, pelvic floor, peritoneum, native vessels/nerves and fixed pelvic skeleton. '
               'Vessels and nerves may occupy the wall but cannot enter the urinary lumen. '
               'No skeleton motion, independent donor overlay, blind lumen disk or overlapping layer faces is allowed.', '',
               'All new anatomical representations and interfaces remain unfinished until actual surfaces, sections, finite cells, '
               'source/normal replay, urinary patency and whole pelvic packing are verified. Source review meshes remain local.', '',
               'Sources: ' + ', '.join(f'[{i + 1}]({url})' for i, url in enumerate(REFERENCES)), '',
               f'Bound source inventories: `{output.relative_to(ROOT)}/male-inventory.json` and `female-inventory.json`.', '']
    report_path = ROOT / 'TA98-BLADDER-PARTS.md'
    if report_path.exists():
        raise FileExistsError(report_path)
    report_path.write_text('\n'.join(report))

    reader = AtlasReader(); sources = {}; review_requests = {}
    for sex, model in MODELS.items():
        folder = output / sex; folder.mkdir()
        a = json.loads(CATALOGUES[sex].read_text())
        ids = [f'TA98PART:ZA:Urinary bladder:A08.3.01.{n:03d}' for n in [3, 4, 5]] if sex == 'male' else ['LOCAL:female:bladder:187']
        pieces = [load_registered(db, pid, model, reader=reader)[1] for pid in ids]
        whole = assembly(pieces, pieces[0].frame, ids)
        whole.save(folder / 'current-bladder.npz')
        for i, m in enumerate(pieces):
            m.save(folder / f'current-region-{i}.npz')
        bounds = whole.bounds; center = bounds.mean(0)
        context_bounds = bounds + np.array([[-.015, -.025, -.015], [.015, .035, .015]])
        entries = [dict(mesh=str(folder / 'current-bladder.npz'), name=f'{sex}: current bladder', color=[.72, .4, .34, 1])]
        for p in a['parts']:
            n = p.get('name', '').lower()
            if p.get('suppressed'):
                continue
            if p['id'].endswith('ureter-wall:adventitia') or n in ['urethra', 'prostatic urethra', 'prostate']:
                m = reader.read_part(a, p, frame=whole.frame, catalogue_path=CATALOGUES[sex], catalogue_sha=digest_file(CATALOGUES[sex]))
                clipped = crop(m, context_bounds)
                if clipped is not None:
                    dest = folder / f'context-{len(entries)}.npz'; clipped.save(dest)
                    entries.append(dict(mesh=str(dest), name=p['name'] + ' (uncapped context crop)',
                                        color=[.6, .64, .8, 1] if 'ureter' in n else [.69, .55, .28, 1]))
        checks = audit(whole, require_closed=True, max_pairs=1000000)
        source = dict(catalogue=str(CATALOGUES[sex]), catalogue_sha256=digest_file(CATALOGUES[sex]),
                      whole_geometry_sha256=whole.geometry_sha, sources=ids, checks=checks,
                      closure=closure(db, ids, model), context_meshes=entries)
        atomic_json(folder / 'source-audit.json', source); sources[sex] = source
        views = [dict(name='front', direction=[0, 0, 1]), dict(name='right', direction=[1, 0, 0]),
                 dict(name='coronal-source-section', direction=[.2, .1, 1],
                      section=dict(normal=[0, 0, 1], point=center.tolist(), clear_outer=True))]
        request = prepare_review(db, dict(frame=whole.frame, meshes=entries, views=views, resolution=800))
        review_requests[sex] = request; atomic_json(folder / 'review-request.json', request)

    raw_path = ROOT / '.local-models/male.json'; raw = json.loads(raw_path.read_text())
    native = next(p for p in raw['parts'] if p['id'] == 'LOCAL:male:bladder:1942')
    m = reader.read_part(raw, native, frame='reference:male-native:unregistered', catalogue_path=raw_path, catalogue_sha=digest_file(raw_path))
    m.save(output / 'male-native-reference.npz')
    t = m.triangles; keep = np.any(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]) != 0, axis=1)
    cleaned = Mesh(m.vertices, m.faces[keep], m.frame, m.attributes,
                   dict(m.source, omitted_zero_area_source_faces=np.flatnonzero(~keep).tolist(),
                        use='Local reference only; not a registered delivery source'))
    cleaned.save(output / 'male-native-clean-reference.npz')
    atomic_json(output / 'male-native-reference-audit.json', dict(
        raw_geometry_sha256=m.geometry_sha, zero_area_source_faces=np.flatnonzero(~keep).tolist(),
        checks=audit(cleaned, require_closed=True, max_pairs=1000000),
        frame=m.frame, anatomical_registration_verified=False, not_for_exact_copy=True))
    request = prepare_review(db, dict(frame=m.frame,
        meshes=[dict(mesh=str(output / 'male-native-clean-reference.npz'), name='Retained male hollow reference', color=[.72, .4, .34, 1])],
        views=[dict(name='front', direction=[0, 0, 1]), dict(name='coronal-reference-section', direction=[0, 0, 1],
               section=dict(normal=[0, 0, 1], point=cleaned.bounds.mean(0).tolist(), clear_outer=True))], resolution=800))
    review_requests['male-native-reference'] = request
    result = dict(schema='bladder-workbench-preparation/v1', inventories={sex: len(inv['terms']) for sex, inv in inventories.items()},
                  current_source_checks={sex: s['checks'] for sex, s in sources.items()},
                  reviews=review_requests, anatomical_acceptance=False, atlas_modified=False)
    atomic_json(output / 'preparation.json', result)
    return result


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', required=True)
    args = p.parse_args(); result = prepare(args.output)
    print(json.dumps(dict(inventories=result['inventories'], reviews={k: v['request_id'] for k, v in result['reviews'].items()}), indent=2))
