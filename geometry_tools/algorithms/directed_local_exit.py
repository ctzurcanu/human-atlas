"""Candidate planes from verified local interior probes and actual tube exits.

This is the directed-exterior construction used by source_control_clearance,
with local oriented probes for open source envelopes. It establishes neither
closed-solid containment nor anatomical fit; actual Float64 contacts remain a
mandatory gate. All inputs, ray distance and active support are explicit.
"""
import numpy as np
from ..open3d_backend import library
from .feature_normals import FeatureNormals,classified_oriented_distance

def directed_local_exit_planes(index,points,active,direction,*,maximum_ray_m,probe_m=.000001):
    points=np.asarray(points,float);active=np.asarray(active,int);direction=np.asarray(direction,float)
    if points.ndim!=2 or points.shape[1]!=3 or direction.shape!=(3,)or abs(np.linalg.norm(direction)-1)>1e-12 or not 0<maximum_ray_m<=.01 or not 0<probe_m<=.00001:
        raise ValueError('Explicit finite points, unit direction and bounded local exit required')
    if not len(active):return np.empty((0,3)),dict(active=0,verifiedLocalProbes=0)
    feature=FeatureNormals(index.mesh);h=index.closest_points(points[active])[0];normal,rim,_=feature.at_masked(h)
    probes=h['points']-normal*probe_m;ph=index.closest_points(probes)[0];pn,pr,_=feature.at_masked(ph);d,k=classified_oriented_distance(probes,ph,pn,pr)
    if not np.all(k&(d<0)):raise ValueError('A directed local exit has an unverified interior probe')
    rays=np.c_[(probes-index.origin)/index.scale,np.broadcast_to(direction,probes.shape)].astype(np.float32)
    hit=index.scene.cast_rays(library().core.Tensor(rays),nthreads=2);distance=hit['t_hit'].numpy()*index.scale;ids=hit['primitive_ids'].numpy()
    if not np.isfinite(distance).all()or np.any(distance<0)or np.any(distance>maximum_ray_m):raise ValueError('Local exit exceeds the declared finite ray bound')
    triangles=index.mesh.triangles[ids];n=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
    if np.any(n@direction<=1e-6):raise ValueError('Directed ray does not leave the locally oriented envelope')
    ray_points=probes+distance[:,None]*direction
    uv=hit['primitive_uvs'].numpy().astype(float);bary=np.c_[1-uv.sum(1),uv]
    exits=np.einsum('ni,nij->nj',bary,triangles)
    error=float(np.linalg.norm(ray_points-exits,axis=1).max())
    if error>1e-7:raise ValueError('Ray display projection exceeds the explicit candidate-plane replay bound')
    return exits,dict(active=len(active),verifiedLocalProbes=len(active),maximumRayTravel_m=float(distance.max()),maximumAllowedRay_m=maximum_ray_m,probe_m=probe_m,direction=direction.tolist(),sourceFaceIds=ids.tolist(),barycentric=bary.tolist(),maximumRayProjectionReplayError_m=error,closedSolidContainmentEstablished=False,actualTriangleContactsRequired=True)
