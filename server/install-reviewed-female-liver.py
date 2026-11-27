"""Activate only the passed, visually reviewed local female liver revision."""
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.db import Database
from geometry_tools import install,jobs
base=ROOT/'.local-models/geometry-workbench/liver-completion';stage=base/'female-liver-install-v207';qa=base/'female-transfer-v206/qa'
result=json.loads((stage/'job-result.json').read_text());assert result['passed']and result['result']['passed'];db=Database();job=jobs.row(db,result['job_id'])
views=[dict(name=name,path=str(path),sha256=digest_file(path),section=section,actual_catalogue_geometry=True,inspected=True)for name,path,section in [('female-exterior',qa/'exterior.png',False),('fixed-diaphragm-context',qa/'fixed-diaphragm-side.png',False),('filled-coronal-channel-section',qa/'coronal.png',True)]]
review=dict(input_signature=job['input_signature'],reviewer='Codex: actual fitted female exterior, fixed diaphragm context and live channel section inspected',accepted_for_local_preview=True,anatomical_acceptance=False,views=views,
 remaining_limits=['Same 63/69 represented liver terms as the male: 57 organ-scale plus six illustrative microscopic terms. Six remain unrepresented.','Regional surface, functional and branch extents remain provisional.','The source internal networks and proximal parents are transferred together; native distal female vascular and biliary junctions remain adjacent/unjoined unless separately verified.','Explicit positive uniform scale reuses literal source geometry buffers; renderer applies the same scale once.'])
atomic_json(stage/'review.json',review);checked=jobs.review(db,result['job_id'],stage/'review.json');print('Review:',checked['state'],flush=True)
proposal=Path(result['result']['proposal']['proposed_catalogue']).with_name('install.json');installed=install.apply(db,proposal,'female-liver-v207');atomic_json(stage/'installed.json',installed);print(json.dumps(installed),flush=True)
