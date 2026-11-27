"""Area/edge-weighted binary chart cut; removes ray-classifier speckle."""
from pathlib import Path
import argparse,sys,time,json
from collections import deque
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
ap=argparse.ArgumentParser();ap.add_argument('--source',default='whole-wall-v25');ap.add_argument('--charts',default='wall-regions-v27');ap.add_argument('--stage',default='lining-charts-v63');args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False);m=Mesh.load(b/args.source/'whole-wall.npz');z=np.load(b/args.charts/'charts.npz');f=m.faces;q=m.triangles;normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);length=np.linalg.norm(normal,axis=1);normal/=length[:,None];area=length*.5*1e6;visibility=np.clip(z['external_visibility'],.005,.995);edges=np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]];face_id=np.tile(np.arange(len(f)),3);keys,inv,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True);assert np.all(counts==2);order=np.argsort(inv,kind='stable').reshape(-1,2);pairs=face_id[order];edge_length=np.linalg.norm(m.vertices[keys[:,1]]-m.vertices[keys[:,0]],axis=1)*1000;weight=.6*edge_length*np.maximum(.05,(1+np.einsum('ij,ij->i',normal[pairs[:,0]],normal[pairs[:,1]]))*.5);source=len(f);sink=source+1;ids=np.arange(len(f));u=np.r_[np.full(len(f),source),ids,pairs[:,0],pairs[:,1]];v=np.r_[ids,np.full(len(f),sink),pairs[:,1],pairs[:,0]];capacity=np.r_[-np.log(visibility)*area,-np.log(1-visibility)*area,weight,weight];capacity=np.maximum(1,np.rint(capacity*10000)).astype('int64');graph=coo_matrix((capacity,(u,v)),shape=(len(f)+2,len(f)+2)).tocsr();result=maximum_flow(graph,source,sink);residual=(graph-result.flow).tocsr();reached=np.zeros(len(f)+2,bool);reached[source]=True;queue=deque([source])
while queue:
 i=queue.popleft()
 for j,c in zip(residual.indices[residual.indptr[i]:residual.indptr[i+1]],residual.data[residual.indptr[i]:residual.indptr[i+1]]):
  if c>0 and not reached[j]:reached[j]=True;queue.append(int(j))
inner=reached[:len(f)];assert inner.any()and not inner.all();atomic_npz(out/'charts.npz',endocardial_faces=inner,face_owner=z['face_owner'],external_visibility=z['external_visibility']);report=dict(source_geometry_sha256=m.geometry_sha,source_chart_sha256=digest_file(b/args.charts/'charts.npz'),source_inner_faces=int(z['endocardial_faces'].sum()),smoothed_inner_faces=int(inner.sum()),changed_faces=int(np.sum(inner!=z['endocardial_faces'])),regularization_width_m=.0006,coordinates_unchanged=True,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False);atomic_json(out/'classification.json',report);print(json.dumps(report),flush=True)
