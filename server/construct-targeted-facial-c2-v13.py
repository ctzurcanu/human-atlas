"""Targeted solve + finite geometry, one bounded source-bound stage."""
from pathlib import Path
import sys,time,json,runpy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';out=base/'male-left-facial-cervical-routes-targeted-c2-v13'
    sys.argv=['optimise-free-c2-facial-course-v9.py','--targeted-blockers'];runpy.run_path(str(ROOT/'server/optimise-free-c2-facial-course-v9.py'),run_name='__main__')
    solve=json.loads((out/'stage1-evidence.json').read_bytes())
    if not solve['all_free_c2_curvature_and_profile_gates_pass']:return
    result=dict(stage=1,maximum_stages=3,source_catalogue_sha256=solve['source_catalogue_sha256'],solve_seconds=time.monotonic()-began,installed=False,blender_queued=False)
    if time.monotonic()-began>33:
        result.update(stopped_before_geometry=True,reason='33-second finite-construction start cutoff',elapsed_seconds=time.monotonic()-began);atomic_json(out/'stage1-construction.json',result);return
    try:
        sys.argv=['construct-facial-cervical-network.py','--route-version','targeted-c2-v13','--output-version','targeted-c2-v13','--max-pairs','1000000','--max-tests','12000000']
        runpy.run_path(str(ROOT/'server/construct-facial-cervical-network.py'),run_name='__main__');result['finite_network_complete']=True
    except Exception as error:result.update(finite_network_complete=False,exact_failure_type=type(error).__name__,exact_failure=str(error))
    result['elapsed_seconds']=time.monotonic()-began;atomic_json(out/'stage1-construction.json',result);print(json.dumps(result))


if __name__=='__main__':main()
