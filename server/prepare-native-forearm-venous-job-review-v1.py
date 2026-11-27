"""Prepare the root's reviewed four-view payload without changing job state."""
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.db import Database
from geometry_tools.jobs import row,artifacts,verify_bindings

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-native-venous-regions-v1/literal-forearm-pack-v1'
JOB='8fb62dcb8565badc97f067310d8ec4eaa16b9e08d52f1db33949434134f5a163'
EXPECTED_REQUESTS={'eaee1f55096474714b02a1435c625590ddeb28e4d276a013e4a775a194892956','9e806c995e1b9d6ba03fb6b0c82e3932a34beb53b68dece92b927402bde3e81e'}


def main():
    gate=json.loads((BASE/'catalogue-candidate-job-gate.json').read_bytes());assert gate['job_id']==JOB and gate['passed']
    review_path=BASE/'root-render-review-v1.json';review=json.loads(review_path.read_bytes());assert review['forearm_accepted_for_local_preview'] is True and review['anatomical_acceptance'] is False
    selected=[r for r in review['requests'] if r['request_id'] in EXPECTED_REQUESTS];assert len(selected)==2
    views=[]
    for request in selected:
        assert request['original_objects_preserved'] is True
        for image in request['views']:
            assert digest_file(image['path'])==image['sha256'];views.append(dict(image,inspected=True,request_id=request['request_id']))
    assert len(views)==4 and sum(bool(v.get('section')) for v in views)==2
    db=Database();job=row(db,JOB);assert job['state']=='passed';artifacts(db,JOB);verify_bindings(job['parameters']['input_files'])
    proposal=gate['result']['proposal'];assert proposal['proposed_sha256']==digest_file(BASE/'proposed-atlas.json')
    spec=json.loads((BASE/'install-spec.json').read_bytes())
    evidence=dict(input_signature=job['input_signature'],reviewer=review['reviewer'],accepted_for_local_preview=True,anatomical_acceptance=False,
        remaining_limits=spec['remaining_limits'],views=views,source_bound_root_review=str(review_path),source_bound_root_review_sha256=digest_file(review_path),
        source_catalogue_sha256=proposal['expected_catalogue_sha256'],proposed_catalogue_sha256=proposal['proposed_sha256'],
        scope='Only the four literal native forearm outer-surface partitions and their unchanged source remainders. No cerebral views, finite wall/lumen or junction acceptance.',
        preparation_only_no_job_state_change=True)
    target=BASE/'job-review-ready.json';atomic_json(target,evidence)
    install=ROOT/f'.local-models/geometry-workbench/jobs/{JOB}/prepared-1/install.json';assert install.exists()
    handoff=dict(job_id=JOB,job_state=job['state'],input_signature=job['input_signature'],review_payload=str(target),review_payload_sha256=digest_file(target),
        prepared_install=str(install),prepared_install_sha256=digest_file(install),proposed_catalogue_sha256=proposal['proposed_sha256'],
        expected_source_catalogue_sha256=proposal['expected_catalogue_sha256'],review_views=4,new_direct_terms=2,
        numerical_geometry_passed=True,review_accepted_by_root_for_local_preview_only=True,root_must_record_job_review=True,root_must_prepare_source_bound_relation_cache=True,
        job_state_unchanged=True,installed=False,anatomical_acceptance=False)
    atomic_json(BASE/'job-ready-handoff.json',handoff);print(json.dumps(handoff))


if __name__=='__main__':main()
