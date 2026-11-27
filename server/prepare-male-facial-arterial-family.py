"""Source-bound facial-family inventory and first construction feasibility inputs.

This exports existing sources; it never changes the live atlas.  All recipient
surfaces remain fixed and are distinct from vessel ownership and blood mouths.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.tube_station_topology import recover_station_rings, align_station_correspondence, StationTopologyError
from geometry_tools.algorithms.tube_continuity import boundary_loops

SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '5f03a758aeea29e282be5a5dc077a663eeacc1c5e4141acf57b57e630ca3642a'
INVENTORY = ROOT / 'reports/ta98-blood-circulation-inventory-installed-v9.json'
OUT = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-facial-family-inputs-v1'
REPORT = ROOT / 'reports/ta98-male/facial-arterial-family-contract-v1.json'


def main():
    start = time.monotonic()
    assert digest_file(SOURCE) == EXPECTED
    assert not OUT.exists() and not REPORT.exists(), 'Retain previous stages'
    census = json.loads(INVENTORY.read_bytes())
    all_terms = {t['id']: t for t in census['terms']}
    selected = {'A12.2.05.020'}
    while True:
        expanded = selected | {k for k, t in all_terms.items() if t['parent_id'] in selected}
        if expanded == selected:
            break
        selected = expanded
    terms = [{k: all_terms[tid][k] for k in ('id', 'name_en', 'name_la', 'parent_id', 'planned_representation_class', 'active_exact_source_part_ids', 'resolved_direct_part_ids', 'representation_status')} for tid in sorted(selected)]
    assert len(terms) == 10
    catalogue = json.loads(SOURCE.read_bytes())
    parts = {p['id']: p for p in catalogue['parts']}
    bases = {
        'facial': 'Facial artery', 'angular': 'Angular artery',
        'submental': 'Submental artery', 'superior-labial': 'Superior labial artery',
        'inferior-labial': 'Inferior labial artery',
        'submandibular-gland': 'Submandibular gland',
        'palatine-tonsil': 'Palatine tonsil',
        'stylopharyngeus': 'Stylopharyngeus muscle',
        'superior-constrictor': 'Superior pharyngeal constrictor',
        'middle-constrictor': 'Middle pharyngeal constrictor',
        'orbicularis-oris': 'Orbicularis oris muscle',
        'submandibular-duct': 'Submandibular duct',
        'facial-vein': 'Facial vein',
    }
    roles = {f'{role}-{side}': f'ZA:{base}.{suffix}' for role, base in bases.items() for side, suffix in [('left', 'l'), ('right', 'r')]}
    roles.update({'soft-palate': 'ZA:Soft palate', 'nasal-septal-cartilage': 'ZA:Nasal septal cartilage', 'nasal-mucosa': 'ZA:Mucosa of nasal cavity'})
    for side, suffix in [('left', 'l'), ('right', 'r')]:
        roles[f'nasalis-{side}'] = f'ZA:Nasalis muscle.{suffix}:Superficial'
        roles[f'nasal-skin-{side}'] = f'ZA:Nasal region.{suffix}:Skin-5'
    reader = AtlasReader()
    frame = 'model:local-ta98:assembled-metres'
    OUT.mkdir(parents=True)
    rows = []
    for role, pid in roles.items():
        if pid not in parts:
            rows.append(dict(role=role, part_id=pid, missing_native_source=True))
            continue
        part = parts[pid]
        assert not part.get('suppressed'), pid
        mesh = reader.read_part(catalogue, part, frame=frame, catalogue_path=SOURCE, catalogue_sha=EXPECTED)
        output = OUT / (role + '.npz')
        mesh.save(output)
        row = dict(role=role, part_id=pid, fixed=True, source_record=part, mesh=str(output.relative_to(ROOT)),
                   mesh_sha256=digest_file(output), geometry_sha256=mesh.geometry_sha,
                   bank=catalogue['chunks'][part['chunk']], bank_sha256=reader.file_sha(reader.asset(catalogue['chunks'][part['chunk']]['url'])),
                   bounds_m=mesh.bounds.tolist(), triangles=len(mesh.faces))
        try:
            loops = boundary_loops(mesh.faces)
            row['source_boundary_loops'] = [dict(indices=r.tolist(), center_m=mesh.vertices[r].mean(0).tolist(), controls=len(r)) for r in loops]
        except (ValueError, AssertionError) as exc:
            row['source_boundary_recovery_pending'] = str(exc)
        if role.split('-left')[0].split('-right')[0] in ['facial', 'angular', 'submental', 'superior-labial', 'inferior-labial']:
            try:
                recovered = []
                for component in recover_station_rings(mesh.vertices, mesh.faces):
                    rings, replay = align_station_correspondence(mesh.vertices, mesh.faces[component['face_ids']], component['rings'])
                    recovered.append(dict(source_faces=component['face_ids'].tolist(), rings=[r.tolist() for r in rings],
                                          centers_m=[mesh.vertices[r].mean(0).tolist() for r in rings], phase_replay=replay))
                row['source_station_components'] = recovered
            except (StationTopologyError, ValueError, AssertionError) as exc:
                row['station_recovery_pending'] = str(exc)
        rows.append(row)
    # The organ is the recipient of blood only: exporting its unchanged tonsil
    # mesh neither edits the lymphatic system nor counts a lymphatic term.
    interfaces = [
        dict(parent='facial-left/right', child='ascending-palatine-left/right', term='A12.2.05.021', recipient='soft-palate', course='Near cervical facial origin, ascend lateral pharynx between styloglossus and stylopharyngeus; enter soft-palate tissue', verified=False),
        dict(parent='facial-left/right', child='tonsillar-left/right', term='A12.2.05.022', recipient='palatine-tonsil-left/right', course='Cervical facial branch passes through superior constrictor into inferior tonsillar pole; keep oral cavity excluded', verified=False),
        dict(parent='facial-left/right', child='glandular-left/right', term='A12.2.05.024', recipient='submandibular-gland-left/right', course='Short independently authored branches enter same-side native gland; gland duct and oral lumen remain excluded', verified=False),
        dict(parent='superior-labial-left/right', child='nasal-septal-left/right', term='A12.2.05.027', recipient='anterior-inferior-nasal-septal-tissue', course='Ascend from upper-lip arterial source into columella/anterior-inferior septal mucoperichondrium; not into cartilage or nasal airway', verified=False),
        dict(parent='distal-facial/angular-left/right', child='lateral-nasal-left/right', term='A12.2.05.028', recipient='nasalis-left/right and nasal-skin-left/right', course='Same-side lateral nose/alar tissue, distal to superior-labial branch; keep eye, cartilage and airway fixed/excluded', verified=False),
    ]
    report = dict(schema='source-bound-facial-arterial-family-contract/v1', source_catalogue=str(SOURCE.relative_to(ROOT)), source_catalogue_sha256=EXPECTED,
                  inventory=str(INVENTORY.relative_to(ROOT)), inventory_sha256=digest_file(INVENTORY), frame=frame,
                  assembly='Both facial arterial families, native parent/continuations/known branches and five missing named branch types',
                  terms=terms, missing_terms=[t['id'] for t in terms if not t['active_exact_source_part_ids']],
                  fixed_sources=rows, explicit_interfaces=interfaces,
                  ownership='Source faces retained once. Whole communicating paired wall first; source-profile side openings before disjoint finite owner cells and one shared emitted bank. No overlapping cylinders or lumen-spanning caps.',
                  fixed_global_structures=['whole skeleton', 'skin', 'eyes', 'nerves', 'fitted heart', 'bilateral lingual finite network', 'cystic artery'],
                  residual_upstream_interfaces=['Existing ECA-to-facial, facial-to-angular, submental and labial source mouth continuity not yet verified or repaired. A passing new branch alone will not establish those interfaces.'],
                  primary_references=[
                      dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11548982/', role='Human cadaveric study of facial branch origin/diameter and relation to facial landmarks; qualitative branch courses, not copied coordinates'),
                      dict(url='https://pubmed.ncbi.nlm.nih.gov/31291033/', role='44 human hemi-sections: facial artery bulb and parapharyngeal branches, pterygoid hamulus relationship'),
                      dict(url='https://pubmed.ncbi.nlm.nih.gov/27568352/', role='Human nasal vascular supply: facial superior-labial septal contribution distinct from posterior septal and ethmoidal arteries')],
                  private_geometry_copied=False, source_geometry_modified=False, installed=False, geometry_acceptance=False,
                  elapsed_seconds=time.monotonic()-start)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(REPORT, report)
    atomic_json(OUT / 'inputs.json', report)
    assert digest_file(SOURCE) == EXPECTED
    print(json.dumps(dict(report=str(REPORT.relative_to(ROOT)), terms=len(terms), missing=report['missing_terms'], exported_sources=sum(not r.get('missing_native_source') for r in rows),
                          native_station_sources=[dict(role=r['role'], components=len(r['source_station_components']), stations=[len(x['rings']) for x in r['source_station_components']]) for r in rows if r.get('source_station_components')],
                          elapsed_seconds=report['elapsed_seconds'])))


if __name__ == '__main__':
    main()
