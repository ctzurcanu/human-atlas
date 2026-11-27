from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.db import Database
from geometry_tools import jobs
stage=ROOT/'.local-models/geometry-workbench/liver-completion/male-liver-completion-install-v201'
spec=json.loads((stage/'job-spec.json').read_text());result=jobs.execute(Database(),spec)
atomic_json(stage/'job-result.json',result)
print(json.dumps(dict(job_id=result['job_id'],state=result['state'],passed=result.get('result',{}).get('passed'))),flush=True)
