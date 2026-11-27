"""Retain reviewed source-mouth extensions instead of overwriting them."""
import hashlib
import json
from pathlib import Path


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inherit_or_retain_mouth_ledger(previous, current, review, *, catalogue_sha,
                                  source_report_sha):
    previous,current,review=map(Path,[previous,current,review])
    old=json.loads(previous.read_text())
    if not current.exists() or json.loads(current.read_text())==old:
        current.write_bytes(previous.read_bytes())
        return dict(inherited=True,added=0)
    new=json.loads(current.read_text())
    if not review.exists():
        raise ValueError('Existing extended mouth ledger needs a source-bound review; refuse to overwrite it')
    evidence=json.loads(review.read_text())
    if not(evidence.get('complete') and evidence.get('passed') and
           evidence.get('catalogueSha256')==catalogue_sha and
           evidence.get('sourceReportSha256')==source_report_sha and
           evidence.get('previousLedgerSha256')==_sha(previous) and
           evidence.get('sourcePortLedgerSha256')==_sha(current)):
        raise ValueError('Extended source-mouth ledger review is stale or incomplete')
    count=len(old['interfaces'])
    if (new['interfaces'][:count]!=old['interfaces'] or
        len(new['interfaces'])<=count or
        {k:v for k,v in new.items()if k!='interfaces'}!={k:v for k,v in old.items()if k!='interfaces'}):
        raise ValueError('Reviewed ledger extension must retain the entire previous ledger literally')
    keys=['first','second','firstRing','secondRing']
    added=[{k:r[k]for k in keys}for r in new['interfaces'][count:]]
    if added!=evidence.get('addedInterfaces'):
        raise ValueError('Added source-mouth declarations differ from their reviewed identities')
    return dict(inherited=False,retainedReviewedExtension=True,added=len(added))
