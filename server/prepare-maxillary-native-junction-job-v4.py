"""Prepare local job inputs only after root views and packed-context review."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def run(out,reviewpath,pending_final_review=False):
 t=time.monotonic();packing=out/'packing-proof.json';p=json.loads(packing.read_text());review=json.loads(reviewpath.read_text());proposed=out/'proposed-atlas.json'
 assert p['complete']and p['passed']and p['source_master_sha256']==H(CAT)
 assert p['proposed_catalogue_sha256']==H(proposed)
 if pending_final_review:
  assert review.get('accepted_for_raw_local_preview')is True and review['source_catalogue_sha256']==H(CAT)and review['candidate_sha256']==p['producer_cohort_sha256']
  context=json.loads((out/'assembly-plan.json').read_text())['raw_context_evidence']
 else:
  assert review['proposed_catalogue_sha256']==H(proposed)and review['packing_proof_sha256']==H(packing)and review['source_master_sha256']==H(CAT)
  assert review.get('accepted_for_local_preview')is True and review.get('actual_packed_context_review_complete')is True and review.get('remaining_limits_accepted_for_local_preview')is True
  context=review['actual_packed_context_evidence']
 assert context and review.get('views')
 for e in context:assert H(e['path'])==e['sha256']
 for v in review['views']:assert v.get('inspected')is True and H(v['path'])==v['sha256']
 assert p['new_parts']==p['new_concepts']==p['new_TA98_term_count_claim']==0 and len(p['changed_parts'])==6
 inventory=out/'registry-authority-inventory.json';iv=json.loads(inventory.read_text());assert iv['term_count']==len(iv['terms'])==38 and iv['source_catalogue_sha256']==H(CAT)
 assembly=json.loads((out/'assembly-plan.json').read_text());assert assembly['whole_organ']and assembly['adjoining_interfaces']and assembly['finite_wall']is False
 assert assembly['remaining_limits']and assembly['full_anatomical_acceptance']is False
 assembly['final_actual_packed_render_review_pending']=pending_final_review
 evidence=[out/'packing-proof.json',out/'compiled-cohort-audit.json',out/'independent-decoded-bank-replay.json',out/'assembly-plan.json',reviewpath]+[Path(e['path'])for e in context]+[Path(v['path'])for v in review['views']]
 spec=dict(kind='catalogue-candidate',parameters=dict(model='local-ta98',proposed_catalogue=str(proposed),proposed_catalogue_sha256=H(proposed),expected_catalogue_sha256=H(CAT),changed_parts=p['changed_parts'],require_closed=False,closed_parts=[],inventory=str(inventory),inventory_sha256=H(inventory),assembly=assembly,evidence=list(map(str,evidence)),max_pairs=600000,max_tests=3000000))
 job=out/'catalogue-job-spec.json';assert not job.exists();write(job,spec);write(out/'job-preparation-proof.json',dict(complete=True,passed=True,source_catalogue_sha256=H(CAT),proposed_catalogue_sha256=H(proposed),job_sha256=H(job),actual_context_and_raw_root_views_bound=True,final_actual_packed_render_review_pending=pending_final_review,geometry_job_started=False,database_changed=False,installed=False,elapsed_seconds=time.monotonic()-t));print(str(job))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True);p.add_argument('--review',type=Path,required=True);p.add_argument('--pending-final-review',action='store_true');a=p.parse_args();run(a.directory.resolve(),a.review.resolve(),a.pending_final_review)
