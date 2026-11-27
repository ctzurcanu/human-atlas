"""Save incremental literal-cell quality evidence for a pulmonary candidate."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_json, digest_file
ap=argparse.ArgumentParser();ap.add_argument('--catalogue',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args()
path=ROOT/args.catalogue;a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False);started=time.monotonic()
report=dict(catalogue_sha256=digest_file(path),cells=[],complete=False)
for pid in a['pulmonarySegmentIntegration']['partIds']:
 m=reader.read_part(a,parts[pid],frame='assembled-metres',catalogue_path=path)
 row=dict(id=pid,audit=audit(m,require_closed=True,max_tests=4000000))
 if not row['audit']['passed']:
  c,s=pair_contacts(m,m,0.,selfcheck=True,max_tests=4000000)
  row['bad_faces']=c['pairs'][np.isin(c['relation'],['crossing','coplanar-overlap'])].tolist()
 report['cells'].append(row);report['elapsed_seconds']=time.monotonic()-started;atomic_json(out/'report.json',report,overwrite=True)
 print(pid, 'PASS'if row['audit']['passed']else row['audit'].get('failures'), 'bad',len(row.get('bad_faces',[])),flush=True)
report.update(complete=True,passed=all(r['audit']['passed']for r in report['cells']),elapsed_seconds=time.monotonic()-started);atomic_json(out/'report.json',report,overwrite=True)
