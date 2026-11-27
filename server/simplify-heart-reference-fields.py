"""Reduce rejected field candidates, then perform complete buffer checks."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();base=ROOT/'.local-models/geometry-workbench/heart';src=base/'typed-field-audited-v660';out=base/args.stage;out.mkdir(exist_ok=False);started=time.monotonic();meta=copy.deepcopy(json.loads((src/'construction.json').read_text()));previous=json.loads((src/'float32-audit.json').read_text())['geometries'];rows={}
for term,g in meta['geometries'].items():
    m=Mesh.load(src/g['file']);evidence=previous[term]
    if not evidence['passed']:
        tolerance=g['fieldConditioning']['samplePitchMetres']*.25;s=md.Manifold(md.Mesh64(np.array(m.vertices,order='C'),np.array(m.faces,dtype='uint64',order='C'))).as_original().simplify(tolerance);b=s.to_mesh64();m=Mesh(np.asarray(b.vert_properties[:,:3]).astype('float32').astype(float),np.asarray(b.tri_verts,dtype=int),'heart-native-local-metres');evidence=audit(m,require_closed=True,max_pairs=500000,max_tests=2000000);g['fieldConditioning']['simplificationToleranceMetres']=tolerance
    m.save(out/g['file']);g.update(vertices=len(m.vertices),triangles=len(m.faces));rows[term]=evidence;atomic_json(out/'float32-audit.json',dict(complete=False,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True);print(term,evidence['passed'],evidence['failures'],len(m.faces),'triangles',flush=True)
atomic_json(out/'construction.json',meta);atomic_json(out/'float32-audit.json',dict(complete=True,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True);failed=[t for t,e in rows.items()if not e['passed']];print('Complete',len(rows),'failed',failed,'elapsed',round(time.monotonic()-started,2),flush=True)
if failed:raise SystemExit(2)
