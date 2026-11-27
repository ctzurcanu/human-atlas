"""Bind exact exported resolvers, then verify the standard current forests."""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file

parser = argparse.ArgumentParser()
parser.add_argument('phase', choices=['bind', 'verify'])
args = parser.parse_args()
start = time.monotonic()
folder = ROOT / '.local-models/geometry-workbench/circulatory-completion/current-relationship-finalization-v1'
audit_path = folder/'resolver-input-audit.json'
audit = json.loads(audit_path.read_text())
spec = importlib.util.spec_from_file_location('final_relation_cache', ROOT/'scripts/structural_relation_cache.py')
cache = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cache)
assert cache.AUDITED_PARENT_ONLY_INPUT_SIGNATURE == audit['parentOnlyInputSignature']
assert cache.parent_only_input_signature(audit['resolverInputs']) == audit['parentOnlyInputSignature']
assert all(digest_file(ROOT/path) == expected for path, expected in audit['resolverInputs'].items())
models = []
for sex, expected in audit['sourceCatalogueSha256'].items():
    source = ROOT/audit['sourcePaths'][sex]
    assert digest_file(source) == expected
    cache_path = ROOT/'.local-models/geometry-workbench/dependency-relations'/expected/(sex+'.json')
    prepared = json.loads(cache_path.read_text())
    assert prepared['sourceSha256'] == expected
    assert prepared['resolverInputs'] == audit['resolverInputs']
    assert cache.resolver_inputs_match(ROOT, prepared)
    if args.phase == 'bind':
        cache.bind_parent_snapshot(ROOT, prepared)
        assert cache.PARENT_INDEX in prepared['resolverInputSnapshots']
        atomic_json(cache_path, prepared, overwrite=True)
    row = dict(sex=sex, sourceCatalogue=str(source), sourceCatalogueSha256=expected,
               cache=str(cache_path), cacheSha256=digest_file(cache_path),
               edges=len(prepared['edges']), resolverInputs=len(prepared['resolverInputs']),
               exactResolverInputsMatch=True, parentSnapshotBound=cache.PARENT_INDEX in prepared['resolverInputSnapshots'])
    if args.phase == 'verify':
        catalogue = json.loads(source.read_text())
        parts = {part['id']: part for part in catalogue['parts']}
        assert len(parts) == len(catalogue['parts'])
        output = ROOT/'reports/structural-dependencies/current'
        edge_path, forest_path = output/(sex+'-anatomical-edges.json'), output/(sex+'-forest.json')
        assert edge_path.read_bytes() == cache_path.read_bytes()
        forest = json.loads(forest_path.read_text())
        assert forest['sourceSha256'] == expected
        assert forest['partCount'] == len(parts)
        assert set(forest['partToGroup']) == set(parts)
        assert forest['anatomicalRelationEvidence'] == dict(sha256=digest_file(edge_path), edges=len(prepared['edges']))
        assert all(edge['from'] in parts and edge['to'] in parts and not edge['coMove'] for edge in prepared['edges'])
        assert digest_file(ROOT/forest['spatialIndex']['file']) == forest['spatialIndex']['sha256']
        row.update(forest=str(forest_path), forestSha256=digest_file(forest_path),
                   completeOwnershipLeaves=len(parts), ownershipGroups=forest['assemblyCount'],
                   sourceFrames=forest['sourceFrameCount'], spatialCandidatePairs=forest['spatialIndex']['pairs'],
                   currentEdgeCacheLiteral=True)
    models.append(row)
result = dict(schema='current-relationship-finalization/'+args.phase+'/v1',
              resolverAudit=str(audit_path), resolverAuditSha256=digest_file(audit_path),
              auditedParentOnlyInputSignature=audit['parentOnlyInputSignature'],
              models=models, geometryMastersLiteral=True, geometryAcceptanceClaimed=False,
              elapsedSeconds=time.monotonic()-start)
atomic_json(folder/(args.phase+'-receipt.json'), result)
print(json.dumps(result))
