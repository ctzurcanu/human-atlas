"""Bound an inward finite-wall field on an explicitly supplied outer sheet.

Incident-plane cones prevent corner offsets from crossing their own source
facets. Opposite-wall rays limit thickness; strict mesh/prism audits remain
required. No anatomical thickness or lumen identity is inferred here.
"""
import numpy as np
from scipy.optimize import linprog

def inward_field(vertices,faces,desired,*,ray_scene=True):
    v=np.asarray(vertices,float);f=np.asarray(faces,int);tri=v[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);length=np.linalg.norm(fn,axis=1)
    if np.any(length==0):raise ValueError('Repair zero-area source facets before constructing a wall')
    unit=fn/length[:,None];sums=np.zeros_like(v);incident=[[]for _ in v]
    for j in range(3):np.add.at(sums,f[:,j],fn)
    for i,t in enumerate(f):
        for k in t:incident[k].append(i)
    directions=np.empty_like(v);margins=[];solves=0
    for i,ids in enumerate(incident):
        n=sums[i]/np.linalg.norm(sums[i]);normals=unit[ids];margin=float(np.min(normals@n))
        if margin<=.05:
            q=linprog([0,0,0,-1],A_ub=np.c_[-normals,np.ones(len(normals))],b_ub=np.zeros(len(normals)),bounds=[(-1,1)]*3+[(0,1)],method='highs')
            if not q.success or q.x[3]<=1e-7:raise ValueError('No strictly inward source normal cone at control '+str(i))
            n=q.x[:3]/np.linalg.norm(q.x[:3]);margin=float(np.min(normals@n));solves+=1
        directions[i]=n;margins.append(margin)
    thickness=np.full(len(v),float(desired));hits=0
    if ray_scene:
        from geometry_tools.open3d_backend import library
        o3d=library();origin=v.mean(0);scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((v-origin).astype('<f4')),o3d.core.Tensor(f.astype('<u4')))
        epsilon=desired*.001;rays=np.c_[v-origin-epsilon*directions,-directions].astype('<f4');q=scene.cast_rays(o3d.core.Tensor(rays),nthreads=1);distance=q['t_hit'].numpy().astype(float)+epsilon;which=q['primitive_ids'].numpy()
        for i,(d,face)in enumerate(zip(distance,which)):
            if np.isfinite(d)and int(face)not in incident[i]and d>0:
                thickness[i]=min(thickness[i],.2*d);hits+=1
    inner=v-thickness[:,None]*directions
    return inner,dict(normal_cone_solves=solves,minimum_incident_plane_margin=min(margins),opposite_wall_ray_bounds=hits,minimum_thickness_m=float(thickness.min()),maximum_thickness_m=float(thickness.max()),thickness_is_construction_parameter=True,requires_emitted_geometry_audit=True)
