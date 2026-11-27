"""Bind the complete pulmonary circulation family before any geometry edits."""
import json
import time
from pathlib import Path

from geometry_tools.common import digest_file

SOURCE = Path('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')
CENSUS = Path('reports/ta98-blood-circulation-inventory-installed-v12.json')
OUT = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
EXPECTED = 'd73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'


def main():
    began = time.monotonic()
    assert digest_file(SOURCE) == EXPECTED
    atlas = json.loads(SOURCE.read_text())
    census = json.loads(CENSUS.read_text())
    evidence = json.loads((OUT / 'evidence.json').read_text())
    terms = [t for t in census['terms'] if t['id'].startswith(('A12.2.01.', 'A12.3.02.'))]
    assert len(terms) == 88
    active_airways = atlas['pulmonaryBranchIntegration']['partIds']
    territories = atlas['pulmonarySegmentIntegration']['partIds']
    arterial = {
        'A12.2.01.105': 'ZA:Anterior segmental artery of right lung',
        'A12.2.01.106': 'ZA:Anterior segmental artery of right lung',
        'A12.2.01.108': 'ZA:Posterior segmental artery of right lung',
        'A12.2.01.109': 'ZA:Posterior segmental artery of right lung',
        'A12.2.01.206': 'ZA:Anterior segmental artery of left lung',
        'A12.2.01.207': 'ZA:Anterior segmental artery of left lung',
        'A12.2.01.209': 'ZA:Posterior segmental artery of left lung',
        'A12.2.01.210': 'ZA:Posterior segmental artery of left lung',
    }
    vein_groups = [
        ('103 104', 'ZA:Apical vein of right lung'),
        ('106 107', 'ZA:Anterior vein of right lung'),
        ('109 110', 'ZA:Posterior vein of right lung'),
        ('111', 'ZA:Right superior pulmonary vein'),
        ('203 204', 'ZA:Superior vein of right lung'),
        ('205', 'ZA:Right inferior pulmonary vein'),
        ('207 208 209', 'ZA:Superior basal vein of right lung'),
        ('303 304', 'ZA:Apicoposterior vein of left lung'),
        ('306 307', 'ZA:Anterior vein of left lung'),
        ('403 404', 'ZA:Superior vein of left lung'),
        ('405', 'ZA:Left inferior pulmonary vein'),
        ('407 408 409', 'ZA:Superior basal vein of left lung'),
    ]
    venous = {'A12.3.02.' + suffix: source for suffixes, source in vein_groups for suffix in suffixes.split()}
    inventory = []
    for term in terms:
        row = {key: term[key] for key in ('id', 'name_en', 'name_la', 'parent_id', 'planned_representation_class',
                                          'representation_status', 'active_exact_source_part_ids',
                                          'resolved_descendant_only_part_ids')}
        row['prospective_source_candidate_ids'] = ([arterial[term['id']]] if term['id'] in arterial else
                                                    [venous[term['id']]] if term['id'] in venous else [])
        row['prospective_source_identity_not_accepted'] = bool(row['prospective_source_candidate_ids'])
        inventory.append(row)
    contract = dict(
        schema='whole-pulmonary-circulation-family-contract/v2', source_catalogue=str(SOURCE),
        source_catalogue_sha256=EXPECTED, source_audit_sha256=digest_file(OUT / 'evidence.json'),
        authority_cache_sha256=census['authority_cache_sha256'], total_terms=88,
        currently_directly_available=53, missing_terms=35,
        arterial_division_candidates=arterial, venous_division_and_collector_candidates=venous,
        fixed_current_airway_ids=active_airways, fixed_current_territory_ids=territories,
        fixed_skeleton=True, fixed_heart=True, fixed_other_vessels=True,
        original_ZA_airways_are_hidden_reference_alternatives=True,
        existing_vessel_to_installed_HRA_airway_registration_unestablished=True,
        tissues_and_interfaces=[
            'Arterial roots to fixed main pulmonary arteries; finite wall and continuous blood compartment require separate evidence.',
            'Segmental arterial source courses to fixed segmental airway/territory identities; never occupy airway lumen.',
            'Venous intrasegmental/intersubsegmental and intersegmental identities require bronchial and drainage topology, not a nearest-plane classifier.',
            'Venous collector portions to unchanged pulmonary-vein trunks and left-atrial interfaces; no generic carrier alias counted as new tissue.',
            'Pulmonary trunk sinus is a space, supravalvular ridge a boundary/three-dimensional feature, not independent tube branches.',
            'Ligamentum arteriosum is connective tissue; a persistent ductus is an explicit variant and must not be inserted as an ordinary adult lumen.',
        ],
        planned_stages=[
            'Recover literal source tube stations and a qualified endpoint proximity graph; retain all source identities and arrays.',
            'Review four arterial source families with their original ZA guides and fixed installed HRA guides; propose eight division assignments only when takeoff/course evidence agrees.',
            'Trace all 23 missing venous divisions/collectors against bronchial identity and native drainage graph; expose ambiguity and missing collector geometry.',
            'Only after anatomical source review, partition shared faces or construct finite communicating tissue with declared matched ports.',
            'Verify full current-context packing, tissue occupancy, source normals, open lumina and rendered sections before a guarded install.',
        ],
        primary_reference_urls=[
            'https://pmc.ncbi.nlm.nih.gov/articles/PMC5035324/',
            'https://pmc.ncbi.nlm.nih.gov/articles/PMC8304484/',
            'https://pmc.ncbi.nlm.nih.gov/articles/PMC9735388/',
            'https://pmc.ncbi.nlm.nih.gov/articles/PMC7139082/',
        ],
        inventory=inventory, geometry_changed=False, catalogue_changed=False, anatomy_accepted=False,
        elapsed_seconds=time.monotonic()-began,
    )
    target = OUT / 'family-contract-v2.json'
    assert not target.exists(), 'Do not overwrite retained evidence'
    target.write_text(json.dumps(contract, indent=2)+'\n')
    report = Path('reports/mesh-audit/REVIEW.md')
    marker = '## Pulmonary circulation whole-family inventory 2026-10-23 — source d73c260c'
    assert marker not in report.read_text(), 'Do not duplicate the checked-in inventory'
    lines = ['\n\n' + marker, '',
             'Before construction: all 88 applicable pulmonary-trunk/arterial and pulmonary-venous terms are inventoried below. '
             'Current direct reference availability is 53/88; it does not establish accepted anatomy. '
             'Thirty-one missing vessel divisions/collectors have candidate native ZA source geometry to investigate. '
             'Four other missing terms concern the pulmonary-trunk sinus/ridge, ligamentum arteriosum and persistent-ductus variant. '
             'The installed HRA bronchial/territory assembly and the active ZA vessels are different registered sources; '
             'vascular interfaces between them are not established. No native geometry is changed by this inventory.', '',
             'Authority/candidate hashes, complete fixed assembly and interface contract: '
             '`.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1/family-contract-v2.json`.', '',
             '| TA98 | Term | Parent | Planned representation | Existing direct source meshes | Missing representation / candidate only |',
             '| --- | --- | --- | --- | --- | --- |']
    for t in inventory:
        existing = ', '.join('`'+x+'`' for x in t['active_exact_source_part_ids']) or '—'
        prospective = ', '.join('`'+x+'`' for x in t['prospective_source_candidate_ids'])
        state = 'No direct representation' if t['representation_status'] == 'missing-direct-representation' else t['representation_status']
        if prospective: state += '; source candidate only: ' + prospective
        lines.append('| '+ ' | '.join([t['id'], t['name_en'], t['parent_id'], t['planned_representation_class'], existing, state])+' |')
    with report.open('a') as handle:
        handle.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(contract=str(target), contract_sha256=digest_file(target), inventory_terms=len(inventory),
                          report=str(report), elapsed_seconds=time.monotonic()-began)))


if __name__ == '__main__':
    main()
