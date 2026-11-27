from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.open3d_backend import library
base=ROOT/'.local-models/geometry-workbench/liver-completion';a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());r=AtlasReader();parts={p['id']:p for p in a['parts']};o3d=library();body=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor(body.vertices.astype(np.float32)),o3d.core.Tensor(body.faces.astype(np.uint32)))
old=json.loads((base/'vessel-containment-v157/evidence.json').read_text())['measuredIntrahepaticBranches'];lookups={}
for role in ['artery','portal','hepatic_vein','biliary']:
    source=r.read_part(a,parts['TA98HEPATIC:male:'+role],frame=body.frame);candidate=Mesh.load(base/f'hepatic-contained-course-v158/{role}.npz');lookups[role]={tuple(p):q for p,q in zip(source.vertices,candidate.vertices)}
for row in old:
    if not row['outside_vertices']:continue
    p=parts[row['id']];m=r.read_part(a,p,frame=body.frame);role=p['hepaticConnectionRole'];used=np.unique(m.faces);source=m.vertices[used];v=np.array([lookups[role][tuple(q)]for q in source]);signed=scene.compute_signed_distance(o3d.core.Tensor(v.astype(np.float32)),nthreads=1,nsamples=7).numpy();i=int(signed.argmax());point=v[i];hit=scene.compute_closest_points(o3d.core.Tensor(point[None].astype(np.float32)),nthreads=1)
    print(json.dumps(dict(term=row['term'],point=point.tolist(),sourcePoint=source[i].tolist(),outside=float(signed[i]),nearest=hit['points'].numpy().tolist()[0],normal=hit['primitive_normals'].numpy().tolist()[0])),flush=True)
