"""Explicit evidence-budget continuation; never rebuild the saved candidate."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';candidate=base/'male-left-facial-cervical-finite-network-free-c2-v10';out=candidate/'self-contact-budget-v11';assert not out.exists();out.mkdir()
    plan=json.loads((base/'male-left-facial-cervical-routes-free-c2-v10/routes.json').read_bytes());source=ROOT/plan['source_catalogue'];expected=plan['source_catalogue_sha256'];assert digest_file(source)==expected
    prior=candidate/'wall-audit.json';old=json.loads(prior.read_bytes());assert not old['self_intersections']['complete_within_pair_budget']
    wall=Mesh.load(candidate/'finite-wall.npz');blood=Mesh.load(candidate/'closed-blood-classifier.npz');assert wall.geometry_sha==old['geometry_sha256']
    receipts=[]
    for name,mesh,genus in [('wall',wall,5),('blood',blood,0)]:
        if time.monotonic()-began>43:
            atomic_json(out/'evidence.json',dict(stopped_before_complete_audits=True,reason='43-second internal validation start cutoff',completed=receipts,installed=False,elapsed_seconds=time.monotonic()-began));return
        report=audit(mesh,require_closed=True,expected_components=1,expected_genus=[genus],max_pairs=1000000,max_tests=12000000)
        atomic_json(out/(name+'-audit.json'),report);receipts.append(dict(name=name,geometry_sha256=mesh.geometry_sha,passed=report['passed'],failures=report['failures'],self_intersections=report['self_intersections']))
        if not report['passed']:break
    passed=len(receipts)==2 and all(r['passed']for r in receipts)
    result=dict(schema='saved-facial-complete-self-contact/v11',source_catalogue_sha256=expected,prior_incomplete_receipt=str(prior.relative_to(ROOT)),prior_incomplete_receipt_sha256=digest_file(prior),
                pair_result_budget=1000000,triangle_test_budget=12000000,complete_wall_and_blood_audits_pass=passed,audits=receipts,
                saved_geometry_reconstructed=False,source_calibres_or_anatomy_changed=False,native_context_pending=True,installed=False,blender_queued=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(source)==expected;atomic_json(out/'evidence.json',result)
    print(json.dumps(dict(passed=passed,audits=[dict(name=r['name'],passed=r['passed'],complete=r['self_intersections']['complete_within_pair_budget'],defects=r['self_intersections']['defect_pairs'],tested=r['self_intersections']['candidate_triangle_pairs_tested'],touching=r['self_intersections']['touching_face_pairs'])for r in receipts],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
