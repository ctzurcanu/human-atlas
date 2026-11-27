"""Correct a prefix-based scope error using the exact installed-v8 term rows."""
from pathlib import Path
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file

start = time.monotonic()
folder = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-retinal-vein-placement-feasibility-v1'
inventory_path = ROOT / 'reports/ta98-blood-circulation-inventory-installed-v8.json'
inventory = json.loads(inventory_path.read_text())
assert digest_file(inventory_path) == '3e4361bf5117f627630af990480313f04543544511896e8123fe8e5eec421471'
assert sum(row['male_applicable'] for row in inventory['terms']) == 1149
terms = {row['id']: row for row in inventory['terms']}
scope_path, evidence_path = folder/'scope-and-authority.json', folder/'evidence.json'
scope, evidence = [json.loads(path.read_text()) for path in [scope_path, evidence_path]]
assert scope['schema'] == 'retinal-cross-chapter-blood-scope/v1'
old_hashes = {path.name: digest_file(path) for path in [scope_path, evidence_path]}
atomic_json(folder/'scope-and-authority-before-scope-correction.json', scope)
atomic_json(folder/'evidence-before-scope-correction.json', evidence)
for row in scope['terms']:
    actual = terms[row['id']]
    assert actual['parent_id'] == row['parent_id']
    assert actual['male_applicable'] and actual['scope'] == 'blood-circulation'
    assert actual['representation_status'] == 'missing-direct-representation'
    row.pop('inside_canonical_A12_chapter_denominator')
    row.pop('additional_applicable_blood_anatomy')
    row.update(authority_code_prefix=row['id'].split('.')[0],
               inside_existing_1149_male_blood_closure=True,
               additional_to_existing_denominator=False, inventory_row=actual)
for name in ['canonical_A12_denominator', 'denominator_is_universal_blood_anatomy',
             'additional_A15_blood_terms', 'additional_A15_status']:
    scope.pop(name)
scope.update(
    schema='retinal-cross-chapter-blood-scope/v2', male_blood_parent_closure_denominator=1149,
    denominator_follows_authority_parent_closure_across_code_chapters=True,
    A15_retinal_venules_already_inside_existing_denominator=7,
    additional_terms_to_add_to_denominator=0,
    A15_retinal_venule_status='Missing distinct geometry; included in the existing denominator, neither available nor installed',
    inventory_path=str(inventory_path), inventory_sha256=digest_file(inventory_path),
    inventory_source_catalogue_sha256=inventory['source_catalogue_sha256'],
    inventory_scope_policy=inventory['scope_policy'])
atomic_json(scope_path, scope, overwrite=True)
evidence.update(
    scope_sha256=digest_file(scope_path), inventory_path=str(inventory_path),
    inventory_sha256=digest_file(inventory_path),
    scope_note='All seven A15 retinal venules are already included in the existing 1,149 male blood parent-closure terms and remain unrepresented.')
atomic_json(evidence_path, evidence, overwrite=True)
receipt = dict(
    schema='retinal-scope-correction/v1',
    error='The initial audit mistook an A15 code prefix for exclusion from the A12-rooted blood parent closure.',
    correction='All seven A15.2.04.044–050 are already included in the 1,149 male-applicable terms; add zero to the denominator.',
    original_artifact_hashes=old_hashes,
    corrected_artifact_hashes={path.name: digest_file(path) for path in [scope_path, evidence_path]},
    exact_inventory_path=str(inventory_path), exact_inventory_sha256=digest_file(inventory_path),
    exact_authority_sha256=inventory['authority_cache_sha256'],
    geometry_registrations_and_source_feasibility_unchanged=True,
    live_inventory_builder_changed=False, elapsed_seconds=time.monotonic()-start)
atomic_json(folder/'scope-correction-receipt.json', receipt)
print(json.dumps(receipt))
