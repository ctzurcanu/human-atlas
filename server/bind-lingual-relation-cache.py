"""Bind only a fully finalized source/resolver cache; never rebuild a forest."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'.local-models/geometry-workbench/circulatory-completion/lingual-controls-relation-cache-v1'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--sex',choices=['male','female'],required=True)
    args=parser.parse_args();began=time.monotonic()
    sha=lambda raw:hashlib.sha256(raw).hexdigest()
    audit_path=FOLDER/'resolver-input-audit.json'
    audit=json.loads(audit_path.read_bytes())
    spec=importlib.util.spec_from_file_location('lingual_parent_snapshot_guard',ROOT/'scripts/structural_relation_cache.py')
    guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
    assert guard.AUDITED_PARENT_ONLY_INPUT_SIGNATURE==audit['parentOnlyInputSignature']==guard.parent_only_input_signature(audit['resolverInputs'])
    source=audit['sources']['prospectiveMale' if args.sex=='male' else 'female']
    assert sha((ROOT/source['path']).read_bytes())==source['sha256']
    path=ROOT/'.local-models/geometry-workbench/dependency-relations'/source['sha256']/(args.sex+'.json')
    original=path.read_bytes();prepared=json.loads(original)
    assert prepared['sourceSha256']==source['sha256']
    assert prepared['resolverInputs']==audit['resolverInputs']
    assert guard.resolver_inputs_match(ROOT,prepared)
    producer=prepared['stagedProducer'];manifest=Path(producer['manifest'])
    assert manifest==FOLDER/(args.sex+'-ordinal-stages')/'manifest.json'
    assert sha(manifest.read_bytes())==producer['manifestSha256']
    m=json.loads(manifest.read_bytes())
    assert m['sourceSha256']==source['sha256'] and m['batchSize']==500
    assert producer['completeOrdinalBatches']==m['batchCount']
    assert len(producer['batchFiles'])==m['batchCount']
    for batch in producer['batchFiles']:
        assert sha(Path(batch['path']).read_bytes())==batch['sha256']
    receipt=FOLDER/(args.sex+'-finalized-cache.json')
    assert not receipt.exists()
    guard.bind_parent_snapshot(ROOT,prepared)
    assert prepared['resolverInputSnapshots'][guard.PARENT_INDEX]==(ROOT/guard.PARENT_INDEX).read_text()
    assert guard.resolver_inputs_match(ROOT,prepared)
    # The exporter already published this complete exact cache. Add only the
    # literal parent snapshot; all numerical edges and source bindings persist.
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix='snapshot-',delete=False) as tmp:
        tmp.write(json.dumps(prepared,separators=(',',':')).encode());temporary=Path(tmp.name)
    os.replace(temporary,path)
    after=json.loads(path.read_bytes())
    assert after['edges']==json.loads(original)['edges']
    assert guard.resolver_inputs_match(ROOT,after)
    result=dict(schema='lingual-controls-source-bound-cache/v1',sex=args.sex,source=source,
                audit_sha256=sha(audit_path.read_bytes()),guard_sha256=sha((ROOT/'scripts/structural_relation_cache.py').read_bytes()),
                parent_only_signature=guard.AUDITED_PARENT_ONLY_INPUT_SIGNATURE,
                cache=str(path.relative_to(ROOT)),cache_before_snapshot_sha256=sha(original),cache_sha256=sha(path.read_bytes()),
                complete_batches=m['batchCount'],ordinal_parts=m['partCount'],edges=len(prepared['edges']),
                resolver_inputs=len(prepared['resolverInputs']),parent_snapshot_sha256=sha(after['resolverInputSnapshots'][guard.PARENT_INDEX].encode()),
                exact_resolver_inputs_match=True,all_batch_payloads_bound=True,forest_or_live_catalogue_modified=False,
                elapsed_seconds=time.monotonic()-began)
    receipt.write_text(json.dumps(result,sort_keys=True)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
