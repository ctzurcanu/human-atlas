"""Independently verify a completed local vascular install and its actual body bank.

Reads immutable source/proposal records and decodes the ordinary viewer catalogue;
does not call installation or revision helpers. Run under geometry_tools.runtime.
"""
import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file
from geometry_tools.model import AtlasReader

FEMALE_SHA = '8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'
CARDIAC_BUNDLES = {
    'male': '375bed821531d0bb6c85f5b330362161c9e50dcef2af40f670292844a03e3040',
    'female': 'f88dd44e18df5fe343d3a9724a68ef4d76d17803b7f0b388f9eb26d5772bded5',
}


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--source-sha256', required=True)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--proposal-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    proposal = args.stage / 'proposed-atlas.json'
    current = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(args.source) == args.source_sha256, 'Source snapshot changed'
    assert digest_file(proposal) == args.proposal_sha256, 'Reviewed proposal changed'
    assert digest_file(current) == args.proposal_sha256, 'Install incomplete or live source changed'
    before = json.loads(args.source.read_bytes())
    actual = json.loads(current.read_bytes())
    assert actual == json.loads(proposal.read_bytes())
    params = json.loads((args.stage / 'catalogue-candidate-job.json').read_bytes())['parameters']
    assert params['expected_catalogue_sha256'] == args.source_sha256
    allowed = set(params['changed_parts'])
    old = {row['id']: row for row in before['parts']}
    new = {row['id']: row for row in actual['parts']}
    assert len(old) == len(before['parts']) and len(new) == len(actual['parts'])
    changed = {key for key in old.keys() | new.keys() if old.get(key) != new.get(key)}
    assert changed == allowed, ('Unexpected physical record changes', changed ^ allowed)
    assert all(new[key] == row for key, row in old.items() if key not in allowed)
    skeleton = {key for key, row in old.items() if row.get('system') == 'skeletal'}
    assert skeleton and not (skeleton & allowed), 'Skeleton was modified'
    heart_proof = json.loads((ROOT / '.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/preview-proof.json').read_bytes())
    heart = {row['part'] for row in heart_proof['physical_tissues']}
    assert len(heart) == 33
    heart.update(key for name, values in before['finiteCardiacIntegration'].items()
                 if name.endswith('PartIds') for key in values)
    assert not (heart & allowed) and all(new[key] == old[key] for key in heart)
    assert actual['chunks'][:len(before['chunks'])] == before['chunks']
    female = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
    assert digest_file(female) == FEMALE_SHA
    for sex, expected in CARDIAC_BUNDLES.items():
        bundle = ROOT / f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/cardiac-label-sites.json'
        assert digest_file(bundle) == expected
        assert json.loads(bundle.read_bytes())['integration']['labeledTermCount'] == 113
    body_path = ROOT / '.local-models/organ-models/spinal-cord/male-body.json'
    body = json.loads(body_path.read_bytes())
    mapping = body['sourceChunkIndices']
    assert len(mapping) == len(body['chunks']) and len(set(mapping)) == len(mapping)
    assert body['sourceChunkCount'] == len(actual['chunks'])
    rows = {row['id']: row for row in body['parts']}
    reader = AtlasReader()
    decoded = []
    for key in sorted(allowed):
        row = rows[key]
        assert dict(row, chunk=mapping[row['chunk']]) == new[key], ('Ordinary body differs', key)
        ordinary = reader.read_part(body, row, frame='model:local-ta98:assembled-metres', catalogue_path=body_path)
        canonical = reader.read_part(actual, new[key], frame=ordinary.frame, catalogue_path=current)
        assert ordinary.geometry_sha == canonical.geometry_sha
        assert ordinary.attributes_sha == canonical.attributes_sha
        assert ordinary.source['bank_sha256'] == canonical.source['bank_sha256']
        decoded.append({'part_id': key, 'geometry_sha256': ordinary.geometry_sha,
                        'bank_sha256': ordinary.source['bank_sha256']})
    terms = sorted({new[key]['ta98Term'] for key in allowed - old.keys()})
    for term in terms:
        assert any(c.get('ta98Term') == term and c.get('ta98Kind') == 'exact'
                   and set(c['elements']) & allowed for c in body['concepts']), ('Missing direct control', term)
    report = dict(schema='independent-installed-vascular-family/v1', passed=True,
                  source_catalogue_sha256=args.source_sha256,
                  installed_catalogue_sha256=args.proposal_sha256,
                  exact_changed_records=sorted(changed), new_owners=len(allowed - old.keys()),
                  canonical_new_owner_terms=terms, unchanged_source_records=len(old) - len(allowed & old.keys()),
                  unchanged_skeleton_records=len(skeleton), unchanged_cardiac_records=len(heart),
                  unchanged_female_catalogue_sha256=FEMALE_SHA,
                  unchanged_cardiac_bundle_sha256=CARDIAC_BUNDLES,
                  ordinary_body_actual_bank_decoding=decoded,
                  anatomical_completion=False, elapsed_seconds=time.monotonic() - started)
    atomic_json(args.output, report)
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ('exact_changed_records', 'ordinary_body_actual_bank_decoding')}, sort_keys=True))


if __name__ == '__main__':
    main()
