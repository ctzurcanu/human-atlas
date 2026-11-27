"""Prepare a hash-bound local geometry job after independent packed replay.

This writes job inputs only. It never runs geometry, installs, or deploys.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1'
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
H = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def run(out):
    started = time.monotonic()
    proof = json.loads((out / 'independent-compiled-replay.json').read_text())
    assert proof.get('complete') and proof.get('passed')
    assert H(CAT) == proof['source_catalogue_sha256']
    assert H(out / 'proposed-atlas.json') == proof['proposed_catalogue_sha256']
    inventory_source = WORK / 'complete-authority-effective-inventory.json'
    inventory = json.loads(inventory_source.read_text())
    assert inventory['schema'] == 'ta98-geometry-inventory/v1'
    assert inventory['term_count'] == len(inventory['terms']) == 24
    assert inventory['source_catalogue_sha256'] == proof['source_catalogue_sha256']
    inventory_path = out / 'registry-authority-inventory.json'
    inventory_path.write_bytes(inventory_source.read_bytes())
    assembly = json.loads((out / 'assembly-plan-before-union.json').read_text())
    assert assembly.get('whole_organ') and assembly.get('adjoining_interfaces')
    assert '2.376 mm gap' in ' '.join(assembly['adjoining_interfaces'])
    assert assembly['finite_wall'] is False and assembly['physiological_patency_claim'] is False
    changed = proof['changed_existing_parts'] + proof['new_parts']
    assert len(changed) == len(set(changed)) == 8
    evidence = [out / name for name in ['packing-proof.json', 'native-tree-proof.json',
        'compiled-merged-audit.json', 'diagnostic-blood-classifier-audit.json',
        'independent-compiled-replay.json', 'assembly-plan-before-union.json']]
    assert all(path.is_file() for path in evidence)
    job = dict(kind='catalogue-candidate', parameters=dict(model='local-ta98',
        proposed_catalogue=str(out / 'proposed-atlas.json'),
        proposed_catalogue_sha256=H(out / 'proposed-atlas.json'),
        expected_catalogue_sha256=proof['source_catalogue_sha256'], changed_parts=changed,
        require_closed=False, closed_parts=[], inventory=str(inventory_path),
        inventory_sha256=H(inventory_path), assembly=assembly,
        evidence=list(map(str, evidence)), max_pairs=600000, max_tests=2000000))
    path = out / 'catalogue-job.json'
    dump(path, job)
    dump(out / 'job-preparation-proof.json', dict(complete=True, passed=True,
        inventory_sha256=H(inventory_path), catalogue_job_sha256=H(path),
        complete_authority_terms=24, changed_parts=changed, geometry_or_live_changes=False,
        elapsed_seconds=time.monotonic() - started))
    print(str(path))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    run(args.directory.resolve())
