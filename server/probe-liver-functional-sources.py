"""Bounded source audit for the eight native functional territories."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.audit import topology

start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion'
out=base/'functional-source-probe-v144';out.mkdir(exist_ok=False)
rows=[]
for number,label in [(0,'VI'),(1,'V'),(2,'III'),(3,'IV'),(4,'II'),(7,'VII'),(8,'VIII'),(9,'I')]:
    mesh=Mesh.load(base/f'common-field-v13/{number:03d}.npz')
    try:
        solid,info=make_solid(mesh.vertices,mesh.faces)
        rows.append(dict(label=label,passed=True,volume=float(solid.volume()),info=info))
    except Exception as exc:
        rows.append(dict(label=label,passed=False,error=str(exc),topology=topology(mesh)))
    print(json.dumps(rows[-1]),flush=True)
atomic_json(out/'evidence.json',dict(parts=rows,elapsedSeconds=time.monotonic()-start))
