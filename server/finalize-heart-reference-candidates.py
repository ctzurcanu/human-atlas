"""Select passing buffers; regularize only new self-touching reference solids.

The maximum dilation is 10 micrometres in the common cardiac chart. Installed
organ/skeleton geometry and accepted shared pericardial strata are immutable.
"""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();base=ROOT/'.local-models/geometry-workbench/heart';out=base/args.stage;out.mkdir(exist_ok=False);started=time.monotonic()
source=base/'typed-finite-stored-v656';meta=copy.deepcopy(json.loads((source/'construction.json').read_text()));previous=[base/'typed-finite-final-v654',base/'typed-packed-final-v658'];accepted=[(p,json.loads((p/'float32-audit.json').read_text())['geometries'])for p in previous];rows={}
for term,g in meta['geometries'].items():
    chosen=next(((p,d[term])for p,d in accepted if term in d and d[term]['passed']),None)
    if chosen:
        p,evidence=chosen;m=Mesh.load(p/g['file']);m=Mesh(m.vertices.astype('float32').astype(float),m.faces,m.frame);g['acceptedSourceStage']=p.name
    else:
        assert g['type']=='finite-mesh'and not term.startswith('A12.1.08.'),'Never dilate an empty space or a shared pericardial stratum'
        m=Mesh.load(source/g['file']);s=md.Manifold(md.Mesh64(np.array(m.vertices,order='C'),np.array(m.faces,dtype='uint64',order='C')))
        s=s.minkowski_sum(md.Manifold.sphere(.00001,4)).as_original().simplify(.000001);b=s.to_mesh64();m=Mesh(np.asarray(b.vert_properties[:,:3]).astype('float32').astype(float),np.asarray(b.tri_verts,dtype=int),'heart-native-local-metres');evidence=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000)
        g['referenceRegularization']=dict(method='finite morphological dilation removes subprecision self-contacts',maximumDilationMetres=.00001,originalAnatomyChanged=False)
    m.save(out/g['file']);g.update(vertices=len(m.vertices),triangles=len(m.faces));rows[term]=evidence
    atomic_json(out/'float32-audit.json',dict(complete=False,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True)
    print(term,evidence['passed'],evidence['failures'],len(evidence['topology']['components']),'components',flush=True)
atomic_json(out/'construction.json',meta);atomic_json(out/'float32-audit.json',dict(complete=True,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True)
failed=[term for term,e in rows.items()if not e['passed']];print('Complete',len(rows),'failed',failed,'elapsed',round(time.monotonic()-started,2),flush=True)
if failed:raise SystemExit(2)
