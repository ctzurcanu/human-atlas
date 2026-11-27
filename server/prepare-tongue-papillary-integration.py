"""Exact data-only inputs for integrating three real papillary host controls.

This prepares neither a complete atlas nor an installation approval. A caller
must retain the source-bound male assembly, replace its common finite bank as
one unit and run the normal reviewed-job installer with its own final manifest.
"""
from pathlib import Path
import argparse
import copy
import json
import sys
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'


def main(candidate,controls,audit,review,output):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local stage')
    construction=candidate/'construction.json';c=json.loads(construction.read_text())
    check=json.loads(audit.read_text());ctrl=json.loads(controls.read_text());visual=json.loads(review.read_text())
    if not c['passed']or not check['passed']or check['construction']['sha256']!=digest_file(construction):
        raise ValueError('Failed or stale construction/neighbor replay')
    if visual['constructionSha256']!=digest_file(construction)or not visual['acceptedForPartialFieldPreview']:
        raise ValueError('The actual finite surface/section review is absent')
    wanted={'A05.1.04.013','A05.1.04.014','A05.1.04.015'}
    new=[copy.deepcopy(r)for r in ctrl['termGroups']if r['term']in wanted]
    if len(new)!=3 or {r['term']for r in new}!=wanted:raise ValueError('Require exactly three actual relief controls')
    for row in new:
        if row['partIds']!=['TA98TONGUE:male:mucosa']or not row['surface']['faces']:
            raise ValueError('No additional solid or empty papilla label permitted')
        # Canonical TA98: .013 belongs under tongue; .014/.015 under .013.
        row['canonicalParentTerm']='A05.1.04.001'if row['term']=='A05.1.04.013'else'A05.1.04.013'
    tissue=candidate/'packed/tissue.json';t=json.loads(tissue.read_text())
    if not t['passed']or digest_file(t['bank'])!=t['bank_sha256']:raise ValueError('Changed actual bank')
    output.mkdir(parents=True)
    atomic_json(output/'integration-inputs.json',dict(schema='tongue-papillary-integration-inputs/v1',
        sex='male',sourceTissueManifest=c['source'],
        construction=dict(path=str(construction),sha256=digest_file(construction)),
        actualTissueManifest=dict(path=str(tissue),sha256=digest_file(tissue)),
        bank=dict(path=t['bank'],sha256=t['bank_sha256'],bytes=t['total_bytes']),
        commonBankMustReplaceTogether=True,partRecords=t['parts'],
        unchangedReciprocalCells=['TA98TONGUE:male:aponeurosis','TA98TONGUE:male:intrinsic:superior-longitudinal'],
        mutatedSurfaceHost='TA98TONGUE:male:mucosa',newPhysicalMeshCount=0,newImageCount=0,
        newTermGroups=new,allRemappedSurfaceControls=dict(path=str(controls),sha256=digest_file(controls)),
        sourceFacetAncestry=ctrl['sourceFaceAncestry'],
        neighborAudit=dict(path=str(audit),sha256=digest_file(audit)),
        actualVisualReview=dict(path=str(review),sha256=digest_file(review)),
        fullDorsalDistribution=False,anatomicalAcceptance=False,installed=False,
        installRequires='Build final complete male catalogue from exact selected assembly; preserve canonical terms and existing controls, run reviewed-job install and original-context checks. This input alone does not authorize anatomical completion.',
        remaining=c['remaining'],elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(integration_inputs=str(output/'integration-inputs.json'),new_term_controls=3,
        new_physical_meshes=0,common_bank_bytes=t['total_bytes'],elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for k in ['candidate','controls','audit','review','output']:ap.add_argument('--'+k,type=Path,required=True)
    a=ap.parse_args();main(a.candidate.resolve(),a.controls.resolve(),a.audit.resolve(),a.review.resolve(),a.output.resolve())
