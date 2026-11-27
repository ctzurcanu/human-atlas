from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.db import Database
from geometry_tools import jobs
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--stage');args=ap.parse_args();stage=ROOT/'.local-models/geometry-workbench/liver-completion'/(args.stage or'final-install-v266-'+args.sex);result=jobs.execute(Database(),json.loads((stage/'job-spec.json').read_text()));atomic_json(stage/'job-result.json',result);print(json.dumps(dict(job_id=result['job_id'],state=result['state'],passed=result.get('result',{}).get('passed'))),flush=True)
