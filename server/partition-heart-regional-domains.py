"""One reciprocal shared boundary complex for actual cardiac tissue domains."""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.native_labelled_tetrahedra import partition_constraints_native
ap=argparse.ArgumentParser();ap.add_argument('--source',default='regional-domains-v77');ap.add_argument('--stage',default='regional-tissues-v78');ap.add_argument('--conform-membranous-interface',action='store_true');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/args.source;out=b/args.stage;out.mkdir(exist_ok=False);z=np.load(source/'fields.npz');s=json.loads((source/'cells.json').read_text());r=partition_constraints_native(z['positions']-z['origin'],z['tetrahedra'],z['fields'],s['constraints'],z['zero_side_priority'],out/'native',seconds=60.,maximum_tetrahedra=1000000,tetrahedron_domains=z['domains'],cell_domains=z['cell_domains'],surface_refinement_planes=[z['fields'].shape[1]-1]if args.conform_membranous_interface else());arrays={'positions':r['positions']+z['origin']};counts={}
for i,name in enumerate(s['names']):arrays['cell_'+str(i)]=r['cells'][i];counts[name]=len(r['cells'][i])
atomic_npz(out/'cells.npz',**arrays);atomic_json(out/'construction.json',dict(constructed=True,source_sha256=digest_file(source/'fields.npz'),names=s['names'],cell_triangles=counts,native=r['evidence'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False));print('Shared cardiac regional tissues',counts,flush=True)
