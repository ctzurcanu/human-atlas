"""Activate the reviewed male liver locally through the journaled adapter."""
from pathlib import Path
import json, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.db import Database
from geometry_tools import install, jobs

base = ROOT / '.local-models/geometry-workbench/liver-completion'
stage = base / 'male-liver-completion-install-v201'
result = json.loads((stage / 'job-result.json').read_text())
assert result['passed'] and result['result']['passed']
db = Database()
job = jobs.row(db, result['job_id'])
views = []
for name, path, section in [
    ('unchanged-exterior-native-entries', base / 'completion-v200/qa/exterior.png', False),
    ('filled-coronal-section', base / 'completion-v200/qa/coronal.png', True),
    ('functional-territories-coronal', base / 'completion-v200/qa/functional-coronal.png', True),
    ('unchanged-illustrative-microscopic-unit', base / 'completion-v155/qa/microscopic-axial.png', True),
]:
    views.append(dict(name=name, path=str(path), sha256=digest_file(path),
        section=section, actual_catalogue_geometry=True, inspected=True))
coverage = json.loads((base / 'completion-v155/coverage.json').read_text())
review = dict(input_signature=job['input_signature'],
    reviewer='Codex: actual exterior, final filled sections, functional selection and unchanged microscopic detail reviewed',
    accepted_for_local_preview=True, anatomical_acceptance=False, views=views,
    remaining_limits=[
        '63/69 TA98 liver terms represented: 57 organ-scale and six illustrative microscopic terms; not 91.3 percent anatomical completion.',
        'Provisional source-bound surface regions, functional boundaries, ligaments and vascular branch ownership.',
        'Native hilar and caval entry/exit segments are retained; containment excludes only literal source-preserved attachment segments.',
        'Missing: ' + ', '.join(r['name'] for r in coverage['missing']),
    ])
atomic_json(stage / 'review.json', review)
checked = jobs.review(db, result['job_id'], stage / 'review.json')
print('Review:', checked['state'], flush=True)
proposal = Path(result['result']['proposal']['proposed_catalogue']).with_name('install.json')
installed = install.apply(db, proposal, 'male-liver-completion-v201')
atomic_json(stage / 'installed.json', installed)
print(json.dumps(installed), flush=True)
