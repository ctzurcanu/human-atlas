"""Audit the actual Float32 reference buffers before installing a sidecar."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args()
stage=ROOT/'.local-models/geometry-workbench/heart'/args.stage
meta=json.loads((stage/'construction.json').read_text());started=time.monotonic();rows={}
for term,g in meta['geometries'].items():
    m=Mesh.load(stage/g['file']);m=Mesh(m.vertices.astype('float32').astype(float),m.faces,'heart-native-local-metres')
    evidence=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000)
    evidence['body_extent_mm']=(np.ptp(m.vertices,axis=0)*.7658069963266287*1000).tolist()
    rows[term]=evidence
    atomic_json(stage/'float32-audit.json',dict(complete=False,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True)
    print(term,'pass',evidence.get('passed'), 'failures',evidence.get('failures'),'components',len(evidence['topology']['components']),'extent mm',np.round(evidence['body_extent_mm'],2).tolist(),flush=True)
atomic_json(stage/'float32-audit.json',dict(complete=True,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True)
print('Audit completed',len(rows),'elapsed',round(time.monotonic()-started,2),'s',flush=True)
