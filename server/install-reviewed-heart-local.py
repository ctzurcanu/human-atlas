"""Install a numerically passed, actually inspected local heart proposal."""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.db import Database
from geometry_tools import jobs,install
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--name',required=True);args=ap.parse_args();stage=ROOT/'.local-models/geometry-workbench/heart'/args.stage;r=json.loads((stage/'job-result.json').read_text());db=Database();job=jobs.row(db,r['job_id']);proposal=r['result']['proposal'];views=[dict(name=name,path=str(stage/(name+'.png')),sha256=digest_file(stage/(name+'.png')),inspected=True,section=name=='section')for name in ['surface','section']];evidence=dict(input_signature=job['input_signature'],reviewer='Codex — actual local atlas exterior and coronal section inspection',accepted_for_local_preview=True,anatomical_acceptance=False,remaining_limits=proposal['assembly']['remaining_limits'],views=views);path=stage/'review.json';atomic_json(path,evidence);review=jobs.review(db,job['id'],path);print('Heart review recorded:',review['state'],flush=True);result=install.apply(db,Path(proposal['proposed_catalogue']).parent/'install.json',args.name);atomic_json(stage/'installation.json',result);print(json.dumps(result),flush=True)
