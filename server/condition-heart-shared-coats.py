"""Propagate bounded core rounding to all shared lining interfaces."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.bounded_float32_rounding import round_diagnosed_crossings
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'chart-coats-v73';out=b/'chart-coats-v74';out.mkdir(exist_ok=False);m=Mesh.load(source/'cell-0.npz');a=Mesh.load(b/'native-pv-mouths-v69/whole-wall.npz');o=np.array([0.,1.375,.03125]);proof=json.loads((source/'core-audit.json').read_text());pairs=proof['self_intersections']['examples'];assert len(pairs)==proof['self_intersections']['defect_pairs'];dtype=np.dtype((np.void,12));stored=(m.vertices-o).astype('<f4').view(dtype).ravel();original=(a.vertices-o).astype('<f4').view(dtype).ravel();protected=np.flatnonzero(np.isin(stored,original));result=round_diagnosed_crossings(m.vertices,m.faces,pairs,origin=o,protected=protected,seconds=20.,maximum_moves=32);atomic_json(out/'rounding.json',result['evidence']);z=np.load(source/'cells.npz');p=z['positions'].copy();lookup={tuple(v):i for i,v in enumerate(p)};changed=np.flatnonzero(np.any(result['positions']!=m.vertices,axis=1))
for i in changed:p[lookup[tuple(m.vertices[i])]]=result['positions'][i]
arrays={'positions':p};metadata=json.loads((source/'construction.json').read_text());names=metadata['names'];rows=[]
for i,name in enumerate(names):
 f=z['cell_'+str(i)];arrays['cell_'+str(i)]=f;used,inv=np.unique(f,return_inverse=True);mesh=Mesh(p[used],inv.reshape(-1,3),m.frame);mesh.save(out/('cell-'+str(i)+'.npz'));rows.append(dict(name=name,geometry_sha256=mesh.geometry_sha,triangles=len(f)))
atomic_npz(out/'cells.npz',**arrays);core=Mesh.load(out/'cell-0.npz');final=audit(core,require_closed=True,max_pairs=3000000,max_tests=22000000);atomic_json(out/'core-audit.json',final);report=dict(constructed=True,passed=final['passed']and result['evidence']['passed'],source_sha256=digest_file(source/'cells.npz'),source_geometry_sha256=a.geometry_sha,names=names,cells=rows,source_exterior_literal=True,shared_interfaces_updated_together=True,uncoated_surface_fraction=metadata['uncoated_surface_area_mm2']/metadata['source_surface_area_mm2'],minimum_live_thickness_m=metadata['minimum_live_thickness_m'],rounding=result['evidence'],core_audit_passed=final['passed'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
