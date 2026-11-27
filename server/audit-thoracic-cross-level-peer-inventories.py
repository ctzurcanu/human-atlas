#!/usr/bin/env python3
"""Aggregate the saved, bounded right intercostal source-mouth screens.

This reports geometric candidates only; it does not declare interfaces or fit meshes.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')
BASE = ROOT / 'thoracic-right-cross-level-move-dependencies-v1.json'
OUT = ROOT / 'thoracic-right-cross-level-peer-aggregate-v1.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    base = json.loads(BASE.read_text())
    peers = [ROOT / f'thoracic-right-cross-level-peer-{i}-mouths-v1.json' for i in range(13)]
    screens = [json.loads(p.read_text()) for p in peers]
    assert base['complete'] and all(x['complete'] and not x['missingBounds'] for x in screens)
    assert all(x['catalogueSha256'] == base['catalogueSha256'] and x['forestSha256'] == base['forestSha256'] for x in screens)
    ids = {level['sourceId'] for level in base['levels']} | {x['sourceId'] for x in screens}
    expected_peers = {c['id'] for level in base['levels'] for port in level['ports'] for c in port['candidates']}
    assert {x['sourceId'] for x in screens} == expected_peers
    external = {}
    summaries = []
    for x in screens:
        unresolved = []
        for row in x['rows']:
            if not row['matches']:
                unresolved.append(row['loop'])
            for m in row['matches']:
                if m['id'] not in ids:
                    external.setdefault(m['id'], []).append({'source': x['sourceId'], 'sourceLoop': row['loop'], 'targetLoop': m['loop'], 'maximumCurveGapMm': m['maximumCurveGapMm']})
        summaries.append({'id': x['sourceId'], 'mouths': len(x['rows']), 'unmatchedLoops': unresolved})
    result = {'complete': True, 'anatomicalAcceptance': False, 'catalogueSha256': base['catalogueSha256'], 'forestSha256': base['forestSha256'], 'baseSha256': digest(BASE), 'peerScreenSha256': [digest(p) for p in peers], 'scopedOwnerCount': len(ids), 'peerMouthCount': sum(s['mouths'] for s in summaries), 'unmatchedPeerMouthCount': sum(len(s['unmatchedLoops']) for s in summaries), 'externalOwnerCount': len(external), 'peerSummaries': summaries, 'externalCandidateOwners': external, 'limitations': ['Equal-control whole-boundary matches among active visible BDM venous sources within 2 mm are geometric candidates, not declared interfaces.', 'Hidden, suppressed, other-source, partial-edge and non-boundary connections are excluded.', 'Unmatched mouths are unclassified, not necessarily free ends; this is not mechanical or anatomical closure.']}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: result[k] for k in ('scopedOwnerCount', 'peerMouthCount', 'unmatchedPeerMouthCount', 'externalOwnerCount')}, indent=2))
    for owner, matches in sorted(external.items()):
        print(owner.split('system-')[-1], len(matches))


if __name__ == '__main__':
    main()
