"""Bind held guide to current source and all declared registered blockers."""
from pathlib import Path
from copy import deepcopy
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file


def main():
    B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=B/'male-left-facial-cervical-routes-palatine-input-v27';assert not O.exists();O.mkdir()
    p=deepcopy(json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes()));binding=D/'source-binding-v19/evidence.json';r=json.loads(binding.read_bytes());assert digest_file(ROOT/r['snapshot'])==r['snapshot_sha256'];inventory=D/'nearby-derivative-context-inventory-v24.json';inv=json.loads(inventory.read_bytes());assert inv['source_catalogue_sha256']==r['snapshot_sha256'];proof=D/'registered-context-v25/evidence.json';gate=json.loads(proof.read_bytes());assert gate['complete']and len(gate['actual_crossings'])==2
    for route in p['routes']:assert digest_file(ROOT/route['guide'])==route['guide_sha256']
    contract=deepcopy(json.loads((ROOT/p['construction_contract']).read_bytes()));contract.update(schema='palatine-all-current-context-contract/v27',source_catalogue=r['snapshot'],source_catalogue_sha256=r['snapshot_sha256'],corrected_routes_only=['ascending-palatine'],literal_passed_routes=['tonsillar','glandular-0','glandular-1'],root_recipient_and_source_calibre_fixed=True,new_actual_crossing_evidence=str(proof.relative_to(ROOT)),new_actual_crossing_evidence_sha256=digest_file(proof),registry_context_inventory=str(inventory.relative_to(ROOT)),registry_context_inventory_sha256=digest_file(inventory),context_registry_revision=inv['registry_revision'],unregistered_frame_or_ownership_ids=sorted(set(x['part_id']for x in inv['unregistered_or_unresolved_candidates'])),unknown_sources_not_counted_passed=True)
    atomic_json(O/'construction-contract.json',contract);atomic_json(ROOT/'reports/ta98-male/facial-palatine-current-context-contract-v27.json',contract)
    p.update(source_catalogue=r['snapshot'],source_catalogue_sha256=r['snapshot_sha256'],literal_context_rebase_proof=str(binding.relative_to(ROOT)),literal_context_rebase_proof_sha256=digest_file(binding),context_registry_revision=inv['registry_revision'],construction_contract=str((O/'construction-contract.json').relative_to(ROOT)),construction_contract_sha256=digest_file(O/'construction-contract.json'));atomic_json(O/'routes.json',p);print(str(O/'routes.json'))


if __name__=='__main__':main()
