"""Prepare a two-record local catalogue job from actual seam-only evidence.

No geometry job, database change, installation or deployment runs here.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
INVENTORY = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1/native-carotid-terminal-cohort-v7b/complete-authority-inventory.json'
IDS = ['TA98PART:ZA:Internal carotid artery.l:A12.2.06.015',
       'TA98PART:ZA:Internal carotid artery.r:A12.2.06.015']
H = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def run(out, evidence):
    started = time.monotonic()
    proposed = out / 'proposed-atlas.json'
    before, after = json.loads(CAT.read_text()), json.loads(proposed.read_text())
    old, new = {p['id']: p for p in before['parts']}, {p['id']: p for p in after['parts']}
    assert old.keys() == new.keys() and len(after['parts']) == len(new)
    changed = [pid for pid in old if old[pid] != new[pid]]
    assert set(changed) == set(IDS), 'Job scope must be exactly the two existing cerebral ICA records'
    for key in ['concepts', 'ta98TermBindings', 'materials', 'defaultHidden']:
        assert before[key] == after[key], 'Seam-only work changed unrelated catalogue metadata: ' + key
    assert after['chunks'][:len(before['chunks'])] == before['chunks']
    inventory = json.loads(INVENTORY.read_text())
    source_sha = H(CAT)
    assert inventory['schema'] == 'ta98-geometry-inventory/v1'
    assert inventory['term_count'] == len(inventory['terms']) == 129
    assert inventory['source_catalogue_sha256'] == source_sha
    assert evidence and all(path.is_file() for path in evidence)
    input_hashes = []
    for path in evidence:
        if path.suffix == '.json':
            value = json.loads(path.read_text())
            declared_source = value.get('source_catalogue_sha256', value.get('master_sha256'))
            assert declared_source in [None, source_sha], 'Evidence belongs to another source revision'
        input_hashes.append(dict(path=str(path), sha256=H(path)))
    inventory_path = out / 'registry-authority-inventory.json'
    inventory_path.write_bytes(INVENTORY.read_bytes())
    assembly = dict(schema='bilateral-ICA-MCA-canonical-seam-only/v9',
        whole_organ='Bilateral existing cerebral ICA portions with only their twelve distal controls and normal words adopting the corresponding native MCA proximal profiles.',
        source_catalogue_sha256=source_sha, changed_parts=IDS,
        adjoining_interfaces=[
            'Each ICA distal twelve-control profile adopts the ipsilateral native MCA canonical positions and normal words; source MCA geometry is fixed.',
            'All other ICA positions and normal words, original ICA UV fields and the proximal ICA mouths remain literal.',
            'The source ICA/MCA corner UV charts remain separate. This cohort changes only the ICA side of the existing geometric boundary.',
            'No ACA root, anterior communicating course, skeleton, brain, nerve or other arterial surface is changed or claimed repaired.',
        ],
        constraints=['Only twelve explicitly corresponding distal ICA controls per side may move.',
                     'All source MCA records and geometry remain literal.',
                     'The fixed skeleton and all other source part records remain literal.'],
        evidence_sha256=input_hashes, finite_wall=False, physiological_patency_claim=False,
        anatomical_acceptance=False, installed=False)
    assembly_path = out / 'catalogue-job-assembly.json'
    dump(assembly_path, assembly)
    parameters = dict(model='local-ta98', proposed_catalogue=str(proposed),
        proposed_catalogue_sha256=H(proposed), expected_catalogue_sha256=source_sha,
        changed_parts=IDS, require_closed=False, closed_parts=[],
        inventory=str(inventory_path), inventory_sha256=H(inventory_path),
        assembly=assembly, evidence=[str(p) for p in evidence] + [str(assembly_path)],
        max_pairs=600000, max_tests=2000000)
    spec = out / 'catalogue-job-spec.json'
    assert not spec.exists(), 'Do not overwrite an existing prepared job spec'
    dump(spec, dict(kind='catalogue-candidate', parameters=parameters))
    dump(out / 'catalogue-job-preparation-proof.json', dict(complete=True, passed=True,
        source_catalogue_sha256=source_sha, proposed_catalogue_sha256=H(proposed),
        job_spec_sha256=H(spec), changed_parts=IDS, unchanged_existing_parts=len(old) - 2,
        complete_authority_terms=129, evidence_sha256=input_hashes,
        geometry_job_started=False, database_changed=False, installed=False,
        elapsed_seconds=time.monotonic() - started))
    print(str(spec))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, action='append', required=True)
    args = parser.parse_args()
    run(args.directory.resolve(), [p.resolve() for p in args.evidence])
