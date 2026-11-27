"""Prepare source-bound inventory/job inputs after independent branch replay."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def run(out,side='left'):
    start=time.monotonic();p=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracoacromial-clavicular-branch-v1';proof=json.loads((out/'independent-compiled-replay.json').read_text());assert proof['passed'] and proof['complete']
    invsource=p/'authority-complete-inventory.json';iv=json.loads(invsource.read_text());source=json.loads((p/'source-inputs-and-constraints.json').read_text());lookup={x['part_id']:x for x in source['sources']};base=json.loads((ROOT/'.local-models/geometry-workbench/circulatory-completion/native-iliac-walls-v1/registry-authority-inventory.json').read_text());terms=[]
    for raw in iv['terms']:
        authority={k:v for k,v in raw.items() if k not in ['classification','existing_inventory','boundary_or_entry_patch_established']};existing=raw.get('existing_inventory',{});ids=existing.get('direct_part_ids',[])
        if raw['id']=='A04.4.01.007':ids=['ZA:Subclavius muscle.l','ZA:Subclavius muscle.r']
        meshes=[dict(part_id=pid,geometry_sha256=lookup[pid]['geometry_sha256'],attributes_sha256=lookup[pid]['attributes_sha256']) if pid in lookup else dict(part_id=pid,source_reference_only=True) for pid in ids]
        terms.append(dict(authority=authority,applicable_to_model=True,representation_class=raw['classification'],source_meshes=meshes,missing_distinct_representation=not bool(ids),distinct_representation_verified=False,source_reference_does_not_certify_correct_anatomy=True))
    inventory={k:base[k] for k in ['schema','authority_path','authority_sha256','authority_corrections_path','authority_corrections_sha256']};inventory.update(model='local-ta98',model_sex='male',root_term='A12.2.09.005',additional_context_root='A04.4.01.007',term_count=len(terms),applicable_term_count=len(terms),terms=terms,current_catalogue_sha256=proof['source_catalogue_sha256'],anatomical_acceptance=False)
    ip=out/'registry-authority-inventory.json';dump(ip,inventory);assembly=json.loads((out/'assembly-plan-before-union.json').read_text());assembly['whole_organ']=f'Complete native {side} thoraco-acromial source exterior with one independently constructed clavicular main/muscular sidecourse and an open parent ostium; source parent proximal/distal profiles retained.'
    evidence=[str(out/x) for x in ['packing-proof.json','native-tree-proof.json','compiled-merged-audit.json','diagnostic-connected-blood-classifier-audit.json','independent-compiled-replay.json','assembly-plan-before-union.json']]
    suffix='l'if side=='left'else'r';job=dict(kind='catalogue-candidate',parameters=dict(model='local-ta98',proposed_catalogue=str(out/'proposed-atlas.json'),proposed_catalogue_sha256=H(out/'proposed-atlas.json'),expected_catalogue_sha256=proof['source_catalogue_sha256'],changed_parts=[f'ZA:Thoraco-acromial artery.{suffix}',f'TA98PART:ZA:Thoraco-acromial artery.{suffix}:A12.2.09.008'],require_closed=False,closed_parts=[],inventory=str(ip),inventory_sha256=H(ip),assembly=assembly,evidence=evidence,max_pairs=600000,max_tests=2000000))
    dump(out/'catalogue-job.json',job);dump(out/'job-preparation-proof.json',dict(passed=True,inventory_sha256=H(ip),job_sha256=H(out/'catalogue-job.json'),complete_authority_terms=7,geometry_or_live_changes=False,elapsed_seconds=time.monotonic()-start));print(str(out/'catalogue-job.json'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True);p.add_argument('--side',choices=['left','right'],default='left');a=p.parse_args();run(a.directory.resolve(),a.side)
