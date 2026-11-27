from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.planar_sliver_flip import flip_planar_slivers
ap=argparse.ArgumentParser();ap.add_argument('--source',default='whole-wall-v18');ap.add_argument('--stage',default='whole-wall-v23');ap.add_argument('--plane-tolerance',type=float,default=2e-8);args=ap.parse_args();assert 0<args.plane_tolerance<=.000002
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/args.source/'whole-wall.npz';out=b/args.stage;out.mkdir(exist_ok=False);m=Mesh.load(source);f,changes=flip_planar_slivers(m.vertices,m.faces,minimum_altitude_m=.00001,plane_tolerance_m=args.plane_tolerance,maximum_flips=200,seconds=8.);mesh=Mesh(m.vertices,f,m.frame);mesh.save(out/'whole-wall.npz');atomic_json(out/'diagonals.json',changes);proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof);report=dict(passed=proof['passed'],source_sha256=digest_file(source),controls_literal=True,flips=len(changes['flips']),remaining_slivers=changes['remaining_below_threshold'],maximum_patch_plane_residual_m=max((q['maximum_plane_residual_m']for q in changes['flips']),default=0),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
