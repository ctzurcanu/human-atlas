from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.db import Database
from geometry_tools import install, jobs

base=ROOT/'.local-models/geometry-workbench/liver-completion'
stage=base/'male-liver-section-install-v142'
result=json.loads((stage/'job-result.json').read_text());jid=result['job_id']
db=Database();job=jobs.row(db,jid)
views=[]
for name,path,section in [
    ('unchanged-exterior',base/'delivery-surface-section-v139/front.png',False),
    ('live-coronal-cut',base/'liver-section-volumes-v141/after-coronal.png',True),
    ('axial-cut',base/'liver-section-volumes-v141/after-axial-cut.png',True)]:
    views.append(dict(name=name,path=str(path),sha256=digest_file(path),section=section,
        actual_catalogue_geometry=True,inspected=True))
evidence=dict(input_signature=job['input_signature'],reviewer='Codex: actual source exterior, local coronal/axial cuts and live drag reviewed',
    accepted_for_local_preview=True,anatomical_acceptance=False,views=views,
    remaining_limits=['Render-only classifiers restore the parenchymal cut face; displayed liver surfaces, normals and adjoining conduits remain literal.',
        'The current liver still represents only 14 of 69 TA98 terms provisionally. This section fix adds no TA98 coverage.',
        'Macroscopic regional extents and microscopic features still require their own construction and anatomical review.'])
atomic_json(stage/'review.json',evidence)
review=jobs.review(db,jid,stage/'review.json');print('Review:',review['state'],flush=True)
proposal=Path(result['result']['proposal']['proposed_catalogue']).with_name('install.json')
installed=install.apply(db,proposal,'male-liver-sections-v142')
atomic_json(stage/'installed.json',installed)
print(json.dumps(installed),flush=True)
