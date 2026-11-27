"""Install the reviewed bilateral source-connection repair; body fit stays provisional.

The source construction report deliberately never authorizes an install. This
adapter requires a later review bound to the emitted catalogue, bank and gates.
It changes only the three existing local BioDigital preview routes, with rollback.
"""
import argparse
import json
import time

from geometry_tools.common import ROOT
from geometry_tools.course_io import sha, write_json
from geometry_tools.install import catalogue_diff, transaction

B = ROOT / '.local-models/geometry-workbench/biodigital-male-circulation-v1'
BASE = B / 'atlas-cardiovascular-review-v32/atlas.json'
OUT = B / 'atlas-cardiovascular-review-v34'
ROUTE = ROOT / 'app/local-ta98-layer-catalogue.ts'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    started = time.monotonic()
    m = json.loads((OUT / 'merge-report.json').read_text())
    review = json.loads((OUT / 'visual-review.json').read_text())
    for path, digest in [(BASE, m['baselineSha256']),
                         (OUT / 'atlas.json', m['catalogueSha256']),
                         (OUT / 'geometry.bin', m['bankSha256']),
                         (ROUTE, m['routeSourceSha256'])]:
        if sha(path) != digest:
            raise ValueError('Prepared install source changed: ' + str(path))
    if not review['accepted_for_local_preview'] or review['anatomical_acceptance']:
        raise ValueError('This adapter requires the explicitly limited local review')
    if not review['remaining_limits'] or review['catalogueSha256'] != m['catalogueSha256']:
        raise ValueError('Missing limits or stale rendered catalogue')
    if review['bankSha256'] != m['bankSha256']:
        raise ValueError('Rendered bank changed')
    for path, digest in review['evidenceBindings'].items():
        if sha(path) != digest:
            raise ValueError('Review evidence changed: ' + path)
    for view in review['views']:
        if sha(view['path']) != view['sha256']:
            raise ValueError('Inspected image changed')
    emitted = OUT / 'emitted-hand-report.json'
    r = json.loads(emitted.read_text())
    audits = json.loads((B / 'hand-emitted-audits-v537/report.json').read_text())
    gate = json.loads((B / 'hand-source-recipient-gate-v536.json').read_text())
    if not r['complete'] or r['skeletonChanged'] or len(r['pairs']) != 90:
        raise ValueError('Incomplete source assembly')
    for evidence in [audits, gate]:
        if not evidence['complete'] or evidence['candidateSha256'] != sha(emitted):
            raise ValueError('Stale or incomplete emitted geometry gate')
    if not audits['allPassed'] or len(audits['rows']) != 34:
        raise ValueError('Source mesh audit failed')
    for row in audits['rows']:
        if not row['passed'] or row['failures'] or sha(row['path']) != row['sha256']:
            raise ValueError('Audited source array changed')
    if any(x['selfPairs'] or x['connectedComponents'] != 1 for x in audits['wholeAssembly']):
        raise ValueError('Whole source network failed')
    if not gate['noContactCountRegressions'] or not gate['remainingWholeHandFitNotAccepted']:
        raise ValueError('Recipient contact scope is not the reviewed limited repair')
    before, after = [json.loads(p.read_text()) for p in [BASE, OUT / 'atlas.json']]
    catalogue_diff(before, after, m['changedParts'])
    for key in ['concepts', 'materials', 'ta98TermBindings', 'defaultHidden', 'sourceCardiacIntegration']:
        if before.get(key) != after.get(key):
            raise ValueError('Unrelated metadata changed: ' + key)
    old = ROUTE.read_bytes()
    needle = b'atlas-cardiovascular-review-v32/atlas.json'
    if old.count(needle) != 3:
        raise ValueError('Unexpected local preview routes')
    if args.apply:
        new = old.replace(needle, b'atlas-cardiovascular-review-v34/atlas.json')
        result = transaction({ROUTE: new}, OUT / 'rollback-route', expected={str(ROUTE): sha(ROUTE)})
        m.update(installed=True, transaction=result, installedRouteSha256=sha(ROUTE),
                 visualReviewSha256=sha(OUT / 'visual-review.json'))
        write_json(OUT / 'merge-report.json', m)
    print(json.dumps(dict(installed=args.apply, changedParts=34, exactSourceJoins=90,
                          sourceSelfPairs=0, skeletonChanged=False,
                          wholeHandBodyFitAccepted=False, elapsedSeconds=time.monotonic()-started)))


if __name__ == '__main__':
    main()
