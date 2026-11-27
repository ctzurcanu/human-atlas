"""Independent per-cell checks with durable bounded-stage evidence."""
from pathlib import Path
import argparse,sys,time,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
ap=argparse.ArgumentParser();ap.add_argument('--source',default='chart-coats-v59');ap.add_argument('--start',type=int,default=1);ap.add_argument('--stop',type=int,default=9);args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart'/args.source;record=json.loads((b/'construction.json').read_text());rows=[]
for i in range(args.start,min(args.stop,len(record['names']))):
 path=b/('cell-'+str(i)+'.npz');m=Mesh.load(path);proof=audit(m,require_closed=True,max_pairs=1300000,max_tests=14000000);atomic_json(b/('cell-'+str(i)+'-audit.json'),proof);row=dict(index=i,name=record['names'][i],source_sha256=digest_file(path),passed=proof['passed'],failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs'],components=len(proof['topology']['components']));rows.append(row);print(row,flush=True)
atomic_json(b/('checks-'+str(args.start)+'-'+str(args.stop)+'.json'),dict(passed=all(r['passed']for r in rows),cells=rows,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False))
