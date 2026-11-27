"""Intrinsic source-shape measurements; no anatomical partition is inferred."""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz
start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';src=h/'apparatus-source-v203';out=h/'chordal-shape-fields-v204';out.mkdir(exist_ok=False);rows=[]
for index in [1,5,8,9,10]:
    source=Mesh.load(src/f'part-{index}.npz');p,inv=np.unique(source.vertices,axis=0,return_inverse=True);f=inv[source.faces]
    edges=np.unique(np.sort(np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]],axis=1),axis=0);length=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)
    graph=coo_matrix((np.r_[length,length],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(p),len(p))).tocsr()
    tri=p[f];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2;weight=np.bincount(f.ravel(),weights=np.repeat(area/3,3),minlength=len(p));eigen=np.full((2,len(p),3),np.nan);radii=[.0015,.0025]
    for first in range(0,len(p),96):
        if time.monotonic()-start>25:raise TimeoutError('Saved shape-field stage has a 25-second native limit')
        distances=dijkstra(graph,directed=False,indices=np.arange(first,min(first+96,len(p))),limit=max(radii))
        for offset,d in enumerate(distances):
            for ri,radius in enumerate(radii):
                ids=np.flatnonzero(d<=radius)
                if len(ids)<4:continue
                xyz=p[ids];w=weight[ids];w=w/w.sum();center=(xyz*w[:,None]).sum(0);q=xyz-center;cov=np.einsum('i,ij,ik->jk',w,q,q);eigen[ri,first+offset]=np.linalg.eigvalsh(cov)
    ratios=eigen[:,:,1]/eigen[:,:,2];atomic_npz(out/f'part-{index}-fields.npz',positions=p,faces=f,source_vertex_map=inv,eigenvalues=eigen,transverse_to_longitudinal_ratio=ratios)
    row=dict(source=index,geometry_sha256=source.geometry_sha,exact_coordinate_welding_only=True,radii_m=radii,ratio_quantiles=[np.nanquantile(v,[0,.1,.25,.5,.75,.9,1]).tolist()for v in ratios],measured_vertices=len(p),complete=bool(np.isfinite(ratios).all()))
    rows.append(row);print(json.dumps(row),flush=True)
atomic_json(out/'measurements.json',dict(measured=True,geometry_modified=False,source=src.as_posix(),parts=rows,anatomical_partition_established=False,elapsed_seconds=time.monotonic()-start))
