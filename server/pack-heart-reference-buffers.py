"""Conform new reference meshes to their literal Float32 coordinate bank."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args()
base=ROOT/'.local-models/geometry-workbench/heart';src=base/args.source;out=base/args.stage;out.mkdir(exist_ok=False);started=time.monotonic();meta=copy.deepcopy(json.loads((src/'construction.json').read_text()))
for term,g in meta['geometries'].items():
    m=Mesh.load(src/g['file'])
    # All mutable controls belong to newly constructed reference geometry.
    # A declared 1-micrometre construction grid removes subprecision Boolean
    # slivers; it never rounds or replaces an installed anatomical source.
    positions=m.vertices if term.startswith('A12.1.08.')else np.round(m.vertices/.000001)*.000001
    v,remap=np.unique(positions.astype('float32'),axis=0,return_inverse=True);f=remap[m.faces];q=v[f].astype(float);keep=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0;f=f[keep]
    first=np.argmin(f,axis=1);cyclic=np.take_along_axis(f,(first[:,None]+np.arange(3))%3,axis=1);sign=np.where(cyclic[:,1]<cyclic[:,2],1,-1);canonical,inv=np.unique(np.sort(f,axis=1),axis=0,return_inverse=True);coefficient=np.bincount(inv,weights=sign,minlength=len(canonical));selected=coefficient!=0;emitted=canonical[selected].copy();emitted[coefficient[selected]<0]=emitted[coefficient[selected]<0][:,[0,2,1]]
    fixed=Mesh(v.astype(float),emitted,m.frame);fixed.save(out/g['file']);g.update(vertices=len(fixed.vertices),triangles=len(fixed.faces),float32Conformity=dict(constructionGridMetres=0 if term.startswith('A12.1.08.')else .000001,maximumPositionErrorMetres=float(np.linalg.norm(v[remap]-m.vertices,axis=1).max()),removedDegenerateFaces=int((~keep).sum()),normalizedDuplicateFaces=int(len(f)-len(emitted)),originalAnatomyChanged=False,requiresFullAudit=True))
    print(term,g['triangles'],'stored triangles',flush=True)
meta.update(source_stage=args.source,elapsed_seconds=time.monotonic()-started)
atomic_json(out/'construction.json',meta)
