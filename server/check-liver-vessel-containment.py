"""Measure actual intrahepatic branch walls against the repaired liver."""
from pathlib import Path
import sys,json,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.open3d_backend import library
from geometry_tools.common import atomic_json

start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'vessel-containment-v157';out.mkdir(exist_ok=False)
current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(current.read_text());r=AtlasReader();parts={p['id']:p for p in a['parts']}
solids=[]
for name in ['right','left','quadrate','caudate','capsule']:
    m=Mesh.load(base/f'closed-liver-cells-v102/{name}.npz');s,_=make_solid(m.vertices,m.faces);solids.append(s)
whole=md.Manifold.batch_boolean(solids,md.OpType.Add);data=whole.to_mesh64();v=np.asarray(data.vert_properties)[:,:3];f=np.asarray(data.tri_verts)
body=Mesh(v,f,'model:local-ta98:assembled-metres');body.save(out/'liver-envelope.npz')
o3d=library();scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor(v.astype(np.float32)),o3d.core.Tensor(f.astype(np.uint32)))
terms=['A12.2.12.'+n for n in ['032','033','034','036','037','038']]+['A12.3.12.'+n for n in ['003','004','007','009','011']]+['A12.3.09.'+n for n in ['006','007','008']]+['A05.8.01.'+n for n in ['063','064','066','067','069']]
rows=[]
for term in terms:
    g=next(g for g in a['hepaticConnectionIntegration']['termGroups']if g['term']==term)
    for pid in g['parts']:
        p=parts[pid];m=r.read_part(a,p,frame=body.frame);ids=np.unique(m.faces);points=m.vertices[ids]
        signed=scene.compute_signed_distance(o3d.core.Tensor(points.astype(np.float32)),nthreads=1,nsamples=7).numpy().astype(float)
        tri_centers=m.triangles.mean(1);centroid_distance=scene.compute_signed_distance(o3d.core.Tensor(tri_centers.astype(np.float32)),nthreads=1,nsamples=7).numpy().astype(float)
        worst=int(signed.argmax());witness=points[worst];lookup=np.full(len(m.vertices),np.nan);lookup[ids]=signed
        record=dict(term=term,id=pid,name=p['name'],vertices=len(ids),outside_vertices=int(sum(signed>3e-7)),
            outside_triangle_centers=int(sum(centroid_distance>3e-7)),wholly_outside_faces=int(np.all(lookup[m.faces]>3e-7,axis=1).sum()),
            maximum_outside_m=float(max(0,signed.max())),furthest_point_m=witness.tolist(),geometrySha256=m.geometry_sha)
        rows.append(record);print(json.dumps(record),flush=True)
atomic_json(out/'evidence.json',dict(sourceCatalogueSha256=__import__('hashlib').sha256(current.read_bytes()).hexdigest(),
    measuredIntrahepaticBranches=rows,containerGeometrySha256=body.geometry_sha,
    expectedExtrahepaticNotTested=['caval trunk','hepatic portal trunk','hepatic artery proper and common hepatic artery','lower common hepatic duct'],
    coordinateGuardMetres=3e-7,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start))
