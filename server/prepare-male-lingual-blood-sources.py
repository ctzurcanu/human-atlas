"""Bind the missing lingual arterial family to the installed native anatomy.

This stage exports fixed sources and anatomical route landmarks. It does not
install vessels or treat a landmark as a branch mouth or perfusion interface.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.tube_station_topology import recover_station_rings, align_station_correspondence, StationTopologyError

SOURCE = ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '9cf8b30456620f7b7bd6f52311ba818fc94f935beee24e6e07cb953d14106f64'
OUT = ROOT/'.local-models/geometry-workbench/circulatory-completion/male-lingual-source-plan-v1'
SHA = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    began = time.monotonic()
    assert SHA(SOURCE) == EXPECTED
    assert not OUT.exists(), 'Retain prior evidence; use a new named stage'
    inventory = ROOT/'reports/ta98-blood-circulation-inventory-installed-v4.json'
    census = json.loads(inventory.read_bytes())
    terms = {r['id']: r for r in census['terms']}
    term_ids = [f'A12.2.05.{n:03d}' for n in range(15, 20)]
    assert all(terms[t]['male_applicable'] for t in term_ids)
    a = json.loads(SOURCE.read_bytes())
    parts = {p['id']: p for p in a['parts']}
    reader = AtlasReader()
    frame = 'model:local-ta98:assembled-metres'
    roles = {
        'external-carotid-left': 'ZA:External carotid artery.l',
        'external-carotid-right': 'ZA:External carotid artery.r',
        'greater-horn-left': 'TA98HYOID:A02.1.16.004:l',
        'greater-horn-right': 'TA98HYOID:A02.1.16.004:r',
        'hyoglossus-left': 'ZA:Hyoglossus muscle.l',
        'hyoglossus-right': 'ZA:Hyoglossus muscle.r',
        'genioglossus-left': 'ZA:Genioglossus muscle.l',
        'genioglossus-right': 'ZA:Genioglossus muscle.r',
        'sublingual-gland-left': 'ZA:Sublingual gland.l',
        'sublingual-gland-right': 'ZA:Sublingual gland.r',
        'tongue-mucosa': 'TA98TONGUE:male:mucosa',
        'inferior-longitudinal-left': 'TA98TONGUE:male:intrinsic:inferior-longitudinal-left',
        'inferior-longitudinal-right': 'TA98TONGUE:male:intrinsic:inferior-longitudinal-right',
    }
    OUT.mkdir(parents=True)
    rows = []
    for role, pid in roles.items():
        p = parts[pid]
        assert not p.get('suppressed')
        m = reader.read_part(a, p, frame=frame, catalogue_path=SOURCE, catalogue_sha=EXPECTED)
        target = OUT/(role+'.npz')
        m.save(target)
        row = dict(role=role, part_id=pid, fixed=True, source_record=p,
                   bank=a['chunks'][p['chunk']], bank_sha256=SHA(reader.asset(a['chunks'][p['chunk']]['url'])),
                   mesh=str(target.relative_to(ROOT)), mesh_sha256=SHA(target), geometry_sha256=m.geometry_sha,
                   bounds_m=m.bounds.tolist())
        if role.startswith('external-carotid-'):
            try:
                components = recover_station_rings(m.vertices, m.faces)
                recovered = []
                for item in components:
                    rings, replay = align_station_correspondence(m.vertices, m.faces[item['face_ids']], item['rings'])
                    recovered.append(dict(source_faces=item['face_ids'].tolist(), rings=[r.tolist() for r in rings],
                                          centres_m=[m.vertices[r].mean(0).tolist() for r in rings],
                                          phase_replay=replay))
                row['source_station_components'] = recovered
            except StationTopologyError as error:
                row['source_station_recovery_pending'] = str(error)
        rows.append(row)
    assert SHA(SOURCE) == EXPECTED
    result = dict(schema='lingual-arterial-construction-inputs/v1', frame=frame,
        source_catalogue=str(SOURCE.relative_to(ROOT)), source_catalogue_sha256=EXPECTED,
        inventory=str(inventory.relative_to(ROOT)), inventory_sha256=SHA(inventory),
        checked_in_inventory='TA98-CIRCULATION.md', checked_in_inventory_sha256=SHA(ROOT/'TA98-CIRCULATION.md'),
        whole_assembly='Bilateral external-carotid to lingual, suprahyoid, dorsal-lingual, sublingual and deep-lingual arteries',
        terms=[{k: terms[t][k] for k in ['id','name_en','name_la','parent_id','planned_representation_class','resolved_direct_part_ids']} for t in term_ids],
        fixed_sources=rows,
        interfaces=[dict(role='proximal-blood-entry', host='external-carotid-left/right', verified=False,
                         requirement='Actual source-profile branch mouth with finite inner/outer walls and no internal membrane'),
                    dict(role='tissue-recipients', hosts=['tongue root and intrinsic muscles','suprahyoid tissues','sublingual gland and floor of mouth'],
                         verified=False, requirement='Branches may enter recipient tissue; vessel lumens must not open into the oral cavity')],
        anatomical_references=[
            dict(url='https://pubmed.ncbi.nlm.nih.gov/27189234/', role='Human cadaveric course classification; medial hyoglossus course and variants'),
            dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC8899351/', role='Human dissection relationships to external carotid, hyoid, hyoglossus and terminal branches'),
            dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC4996103/', role='Human intralingual artery distribution and depth constraints')],
        private_reference_geometry_copied=False, source_geometry_modified=False,
        route_or_mouth_anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-began)
    (OUT/'inputs.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT), fixed_sources=len(rows), terms=term_ids,
                          recovered_carotids=[r['role'] for r in rows if r.get('source_station_components')],
                          elapsed_seconds=result['elapsed_seconds'])))


if __name__ == '__main__':
    main()
