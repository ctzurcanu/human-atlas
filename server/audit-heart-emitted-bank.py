"""Bounded, resumable checks of the literal shared cardiac render bank."""
from pathlib import Path
from types import SimpleNamespace
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.tissue import decode_packed,validate_cells
from geometry_tools.model import AtlasReader
ap=argparse.ArgumentParser();ap.add_argument('--bank',required=True);ap.add_argument('--start',type=int,default=0);ap.add_argument('--stop',type=int);ap.add_argument('--report-name');args=ap.parse_args()
start=time.monotonic();root=ROOT/'.local-models/geometry-workbench/heart'/args.bank;e=json.loads((root/'emission.json').read_text());bank=Path(e['bank']);assert digest_file(bank)==e['bank_sha256'];packed=SimpleNamespace(data=bank.read_bytes(),cells=e['parts']);meshes=decode_packed(packed,e['frame'],np.asarray(e['origin']));
if e.get('catalogue'):
 catalogue=Path(e['catalogue']);a=json.loads(catalogue.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();meshes={name:reader.read_part(a,parts[pid],frame=e['frame'],catalogue_path=catalogue)for name,pid in e['part_ids'].items()}
names=list(meshes);selected=names[args.start:args.stop];checks=validate_cells({k:meshes[k]for k in selected},max_pairs=1300000,max_tests=14000000,checkpoint_directory=root/'checkpoints');atomic_json(root/(args.report_name or f'checks-{args.start}-{args.stop}.json'),dict(checks,source_bank_sha256=e['bank_sha256'],elapsed_seconds=time.monotonic()-start));print(json.dumps(dict(passed=checks['passed'],failures=checks['failures'],cells={k:dict(passed=v['passed'],failures=v['failures'],self_crossings=v['self_intersections'].get('defect_pairs'))for k,v in checks['cells'].items()},pairs=len(checks['pair_occupied_volume_checks']),elapsed_seconds=time.monotonic()-start)),flush=True)
