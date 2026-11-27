"""Final single bounded construction/context stage for free C2 courses."""
from pathlib import Path
import json,sys,time,runpy,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';route=base/'male-left-facial-cervical-routes-free-c2-v10';candidate=base/'male-left-facial-cervical-finite-network-free-c2-v10'
    stage2=json.loads((route/'stage2-evidence.json').read_bytes());assert stage2['all_free_c2_curvature_and_profile_gates_pass'];expected=stage2['source_catalogue_sha256'];source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(source)==expected
    try:
        sys.argv=['construct-facial-cervical-network.py','--route-version','free-c2-v10','--output-version','free-c2-v10']
        runpy.run_path(str(ROOT/'server/construct-facial-cervical-network.py'),run_name='__main__')
        finite_seconds=time.monotonic()-began
        if finite_seconds>48:
            result=dict(stage=3,maximum_stages=3,finite_network_complete=True,context_not_started_due_internal_stage_budget=True,installed=False,blender_queued=False,elapsed_seconds=finite_seconds)
        else:
            work=min(48.,max(1.,55-finite_seconds))
            sys.argv=['check-facial-cervical-feasibility.py','--route-version','free-c2-v10','--network-version','free-c2-v10','--occupied-volume','--work-seconds',str(work)]
            runpy.run_path(str(ROOT/'server/check-facial-cervical-feasibility.py'),run_name='__main__')
            context=json.loads((candidate/'feasibility-v2/evidence.json').read_bytes())
            result=dict(stage=3,maximum_stages=3,source_catalogue_sha256=expected,finite_network_complete=True,finite_construction_seconds=finite_seconds,
                        complete_guide_segments=context['complete_guide_segments'],wall_guide_intersections=len(context['wall_guide_intersections']),
                        minimum_midpoint_blood_winding=context['minimum_midpoint_blood_winding'],native_context_pairs=context['native_context_branch_pairs'],
                        forbidden_crossings=context['forbidden_crossings'],forbidden_occupied_volume_intersections=context['forbidden_occupied_volume_intersections'],
                        unverified_native_occupied_volumes=context['unverified_native_occupied_volumes'],stopped_before_complete_native_context=context['stopped_before_complete_native_context'],
                        numerical_feasibility_pass=context['numerical_feasibility_pass'],emitted_owner_bank_pending=True,
                        full_pharyngeal_airspace_proof_pending=context['oral_lumen_exclusion_full_proof_pending'],installed=False,blender_queued=False,elapsed_seconds=time.monotonic()-began)
    except Exception as error:
        result=dict(stage=3,maximum_stages=3,source_catalogue_sha256=expected,exact_failure_type=type(error).__name__,exact_failure=str(error),
                    construction_or_context_incomplete=True,installed=False,blender_queued=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(source)==expected;atomic_json(route/'stage3-evidence.json',result)
    print(json.dumps({k:v for k,v in result.items()if k!='unverified_native_occupied_volumes'}))


if __name__=='__main__':main()
