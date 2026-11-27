"""Finite peritoneal coats within the existing diagnostic liver envelope.

The bare area remains fibrous capsule. All fractions are illustrative layer
allocation, not measured human histological thicknesses.
"""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import AtlasReader
from geometry_tools.harmonic_wall import tetrahedralize
from geometry_tools.implicit_cells import Arrangement
from geometry_tools.open3d_backend import library
from geometry_tools.tissue import pack
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces

start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'covering-coats-v148';out.mkdir(exist_ok=False)
a=json.loads((base/'liver-surface-controls-v143/with-landmark.json').read_text());part=next(p for p in a['parts']if p['id']=='TA98PART:ZA:Liver:capsule');m=AtlasReader().read_part(a,part,frame='model:local-ta98:assembled-metres')
outer=m.faces[:part['tissueDrawGroups'][0]['count']//3];inner=np.concatenate([m.faces[g['start']//3:(g['start']+g['count'])//3]for g in part['tissueDrawGroups']if g.get('neighbor')])
p,t,*_,te=tetrahedralize(m.vertices,m.faces);atomic_json(out/'tetrahedral-evidence.json',te)
o3d=library()
def distances(f):
    scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor(m.vertices.astype(np.float32)),o3d.core.Tensor(f.astype(np.uint32)))
    return scene.compute_distance(o3d.core.Tensor(p.astype(np.float32)),nthreads=1).numpy().astype(float)
do,di=distances(outer),distances(inner)
# Literal constrained boundaries receive their exact potentials.
lookup={tuple(q):i for i,q in enumerate(p)}
for i in np.unique(outer):do[lookup[tuple(m.vertices[i])]]=0.
for i in np.unique(inner):di[lookup[tuple(m.vertices[i])]]=0.
bare=next(g for g in a['liverFeatureIntegration']['termGroups']if g['term']=='A05.8.01.008')['surface']['faces'];covered=np.setdiff1d(np.arange(len(outer)),bare)
db,dc=distances(outer[bare]),distances(outer[covered]);bare_field=dc-db
values=np.c_[.8*do-.2*di,.55*do-.45*di,bare_field]
def label(x):return'capsule'if x[2]>=0 or x[1]>=0 else'subserosa'if x[0]>=0 else'serosa'
center=np.array([0.,1.25,0.]);arr=Arrangement(p-center,values,label);v,surfaces,cells=arr.build(t);v+=center
mapped,_,conform=conform_quantized_surfaces(v,surfaces,center);surfaces=mapped
atomic_npz(out/'arrangement-source.npz',positions=v,**{'surface_'+str(i):f for i,f in enumerate(surfaces.values())})
ids={name:'TA98LIVER:male:coat:'+name for name in cells}
result=pack(v,surfaces,cells,ids,m.frame,out/'packed',origin=center,max_pairs=2000000,max_tests=16000000,checkpoint_directory=out/'checkpoints')
atomic_json(out/'evidence.json',dict(passed=result['passed'],coats=list(cells),existingEnvelopePreserved=True,bareAreaExcludedFromPeritonealCoats=True,
    sourceCapsuleGeometrySha256=m.geometry_sha,allocationFractions=[.2,.45],thicknessIllustrative=True,conform=conform,elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],bytes=result['total_bytes'],elapsedSeconds=time.monotonic()-start)),flush=True)
