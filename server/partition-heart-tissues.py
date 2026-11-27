"""Build complementary finite tissues with the bounded native PLC engine."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.native_labelled_tetrahedra import partition_constraints_native
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='tissues-v31');ap.add_argument('--source',default='tissue-fields-v30');args=ap.parse_args()
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';source=base/args.source;out=base/args.stage;out.mkdir(exist_ok=False);z=np.load(source/'fields.npz');spec=json.loads((source/'cells.json').read_text());origin=z['origin'];result=partition_constraints_native(z['positions']-origin,z['tetrahedra'],z['fields'],spec['constraints'],z['zero_side_priority'],out/'native',seconds=90.,maximum_tetrahedra=1000000);arrays={'positions':result['positions']+origin};counts={}
for i,name in enumerate(spec['names']):arrays['cell_'+str(i)]=result['cells'][i];counts[name]=len(result['cells'][i])
atomic_npz(out/'cells.npz',**arrays);report=dict(constructed=True,source_field_sha256=digest_file(source/'fields.npz'),positions=len(result['positions']),cell_triangles=counts,cells=spec['names'],native_evidence=result['evidence'],elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
