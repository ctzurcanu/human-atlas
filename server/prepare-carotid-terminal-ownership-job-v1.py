"""Prepare a reviewed source-ownership job; no geometry/DB/install execution."""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
INVENTORY = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1/native-carotid-terminal-cohort-v7b/complete-authority-inventory.json'
H = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def run(out, review_path):
    started = time.monotonic()
    proof_path = out / 'independent-compiled-replay.json'
    proof = json.loads(proof_path.read_text())
    review = json.loads(review_path.read_text())
    proposed = out / 'proposed-atlas.json'
    source = H(CAT)
    assert proof.get('complete') and proof.get('passed')
    assert proof['source_catalogue_sha256'] == source == review['source_catalogue_sha256']
    assert proof['proposed_catalogue_sha256'] == H(proposed) == review['proposed_catalogue_sha256']
    assert review['independent_verification_sha256'] == H(proof_path)
    assert review['accepted_for_local_preview'] is True
    assert review.get('views') and all(v.get('inspected') is True and H(v['path']) == v['sha256'] for v in review['views'])
    assert len(proof['changed_parts']) == 6 and proof['new_concepts'] == 0
    assert proof['all_assembled_native_source_corners_literal']
    inventory = json.loads(INVENTORY.read_text())
    assert inventory['term_count'] == len(inventory['terms']) == 129
    assert inventory['source_catalogue_sha256'] == source
    ip = out / 'registry-authority-inventory.json'; ip.write_bytes(INVENTORY.read_bytes())
    assembly = dict(schema='native-terminal-source-ownership/v12',
        whole_organ='Unchanged bilateral native MCA-labelled surfaces partitioned into terminal ICA common stems and native M1 arterial exteriors.',
        adjoining_interfaces=[
            'The two owners on each side share the literal native twelve-control ring17, its normal words and UVs. All original oriented source face corners and both original exterior mouths remain literal.',
            'The new common-stem piece belongs to the existing cerebral ICA .015 concept; the new M1 piece belongs to its existing native M1 concept with an explicit A12.2.07.047 binding.',
            'The original approximately 193 micrometre ICA/MCA seam phase mismatch remains unchanged. This is source ownership correction, not a repaired communicating lumen.',
            'ACA roots, ACom course, native foreign-vessel contacts, nerves, brain and the skeleton remain unchanged. There are no new surface crossings because the assembled source exteriors are literally identical.',
        ], source_catalogue_sha256=source, provisional_boundary_station=17,
        universal_clinical_boundary_claim=False, finite_wall=False, lumen_repair_claim=False,
        global_new_TA98_reference_term_count_claim=0, review_sha256=H(review_path),
        anatomical_acceptance=False, installed=False)
    ap = out / 'catalogue-job-assembly.json'; dump(ap, assembly)
    evidence = [proof_path, out / 'packing-proof.json', ap, review_path]
    evidence.extend(Path(v['path']) for v in review['views'])
    spec = dict(kind='catalogue-candidate', parameters=dict(model='local-ta98',
        proposed_catalogue=str(proposed), proposed_catalogue_sha256=H(proposed),
        expected_catalogue_sha256=source, changed_parts=proof['changed_parts'],
        require_closed=False, closed_parts=[], inventory=str(ip), inventory_sha256=H(ip),
        assembly=assembly, evidence=[str(p) for p in evidence], max_pairs=600000, max_tests=2000000))
    job = out / 'catalogue-job-spec.json'; assert not job.exists(); dump(job, spec)
    dump(out / 'catalogue-job-preparation-proof.json', dict(complete=True, passed=True,
        source_catalogue_sha256=source, proposed_catalogue_sha256=H(proposed),
        job_spec_sha256=H(job), changed_parts=proof['changed_parts'], complete_authority_terms=129,
        geometry_job_started=False, database_changed=False, installed=False,
        elapsed_seconds=time.monotonic() - started))
    print(str(job))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--review', type=Path, required=True)
    args = parser.parse_args()
    run(args.directory.resolve(), args.review.resolve())
