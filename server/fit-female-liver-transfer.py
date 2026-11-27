"""Bounded, orientation-preserving registration to the native female liver."""
from pathlib import Path
import json,sys,time
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.fit import transform_mesh

started=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion/female-transfer-v202'
plan=json.loads((base/'plan.json').read_text());mp=Path(plan['sourceCatalogue']);fp=Path(plan['targetCatalogue'])
assert digest_file(mp)==plan['sourceCatalogueSha256'] and digest_file(fp)==plan['targetCatalogueSha256']
m=json.loads(mp.read_text());f=json.loads(fp.read_text());reader=AtlasReader()
host=next(p for p in m['parts']if p['id']=='TA98PART:ZA:Liver:capsule')
full=reader.read_part(m,host,frame='model:local-ta98:assembled-metres')
outer=Mesh(full.vertices,full.faces[:host['tissueDrawGroups'][0]['count']//3],full.frame)
target=reader.read_part(f,next(p for p in f['parts']if p['id']==plan['targetReference']),frame='model:local-female-ta98:assembled-metres')
origin=target.bounds.mean(0)
tm=o3d.t.geometry.TriangleMesh(o3d.core.Tensor((target.vertices-origin).astype('float32')),o3d.core.Tensor(target.faces.astype('uint32')))
scene=o3d.t.geometry.RaycastingScene();scene.add_triangles(tm)
tri=outer.triangles;areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
rng=np.random.default_rng(202);ids=rng.choice(len(tri),size=4096,p=areas/areas.sum());w=rng.random((4096,2));w[w.sum(1)>1]=1-w[w.sum(1)>1]
x=tri[ids,0]+w[:,0,None]*(tri[ids,1]-tri[ids,0])+w[:,1,None]*(tri[ids,2]-tri[ids,0])
source_size=np.ptp(outer.bounds,axis=0);target_size=np.ptp(target.bounds,axis=0)
scale=float(source_size@target_size/(source_size@source_size));translation=origin-scale*outer.bounds.mean(0);initial=(scale,translation.copy());iterations=[]
for iteration in range(8):
    q=scale*x+translation
    closest=scene.compute_closest_points(o3d.core.Tensor((q-origin).astype('float32')))['points'].numpy().astype(float)+origin
    distance=np.linalg.norm(q-closest,axis=1);keep=distance<=np.quantile(distance,.9)
    xx=x[keep];yy=closest[keep];cx=xx.mean(0);cy=yy.mean(0)
    # Keep the three native extent measurements as scale constraints. A
    # one-sided surface objective alone can shrink into a concave target.
    proposed=initial[0];assert .62<=proposed<=.72
    dt=cy-proposed*cx;assert np.linalg.norm((dt+proposed*outer.bounds.mean(0))-origin)<.015
    iterations.append(dict(iteration=iteration,scale=proposed,rms_m=float(np.sqrt(np.mean(distance**2))),p95_m=float(np.quantile(distance,.95))))
    scale,translation=proposed,dt
matrix=np.eye(4);matrix[:3,:3]*=scale;matrix[:3,3]=translation
candidate=transform_mesh(outer,matrix,'model:local-female-ta98:assembled-metres');candidate.save(base/'fitted-exterior.npz');target.save(base/'native-female-reference.npz')
actual=scale*x+translation;distance=scene.compute_distance(o3d.core.Tensor((actual-origin).astype('float32'))).numpy()
record=dict(scale=scale,translation_m=translation.tolist(),matrix=matrix.tolist(),rotationDegrees=0,scaleBounds=[.62,.72],iterations=iterations,
    sourceCatalogueSha256=digest_file(mp),targetCatalogueSha256=digest_file(fp),sourcePartId=host['id'],targetPartId=plan['targetReference'],
    fitRmsMetres=float(np.sqrt(np.mean(distance**2))),fitP95Metres=float(np.quantile(distance,.95)),scaleFromThreeNativeExtentMeasurements=True,sourceSurfaceGeometrySha256=outer.geometry_sha,targetGeometrySha256=target.geometry_sha,
    fixedSkeleton=True,allMovingPartsUseOneTransform=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-started)
atomic_json(base/'fit.json',record);print(json.dumps({k:record[k]for k in ['scale','translation_m','fitRmsMetres','fitP95Metres','elapsedSeconds']}))
