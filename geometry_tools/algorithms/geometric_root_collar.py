"""Shared physical sections for a declared, root-connected hollow stem.

The caller identifies the source component and its distal plane before using
this module. Returning distal branches are never selected by axis position.
No cap triangles or new anatomical structures are constructed here.
"""
import numpy as np
from .annular_pipe_envelopes import boundary_loops, face_components, projected_area
from .selective_plane_refinement import refine_selected_plane


def subdivide_geometric_collar(vertices, faces, root_region, *, normal, origin,
                               start, end, levels=(0., .25, .5, .75, 1.),
                               tolerance_m=2e-12):
    p=np.asarray(vertices,float);f=np.asarray(faces,int)
    region=np.asarray(root_region,bool);n=np.asarray(normal,float)
    o=np.asarray(origin,float);levels=np.asarray(levels,float)
    if region.shape!=(len(f),) or not region.any():
        raise ValueError('An explicit nonempty source root region is required')
    if n.shape!=(3,) or o.shape!=(3,) or not np.isfinite(np.r_[n,o,start,end]).all() or np.linalg.norm(n)==0 or end<=start:
        raise ValueError('A finite ordered physical collar is required')
    if levels.ndim!=1 or len(levels)<2 or levels[0]!=0 or levels[-1]!=1 or np.any(np.diff(levels)<=0):
        raise ValueError('Ordered physical stations must include zero and one')
    n=n/np.linalg.norm(n);parents=np.arange(len(f));bary=np.broadcast_to(np.eye(3),(len(f),3,3)).copy()
    initial_rings=boundary_loops(f[region])
    if face_components(f[region])[0]!=1 or len(initial_rings)!=2:
        raise ValueError('The declared stem must be connected with two outlet profiles')
    if any(np.max(abs((p[r]-o)@n-end))>tolerance_m for r in initial_rings):
        raise ValueError('Inner and outer outlets must share the distal physical plane')
    steps=[]
    for level in levels[:-1]:
        coordinate=float(start+level*(end-start))
        z=refine_selected_plane(p,f,normal=n,coordinate=coordinate,origin=o,
            active_faces=region,plane_tolerance_m=tolerance_m)
        parent=z['parent_faces'];bary=np.einsum('nij,njk->nik',z['corner_barycentrics'],bary[parent])
        parents=parents[parent];region=region[parent];p=z['positions'];f=z['faces'];steps.append(z['evidence'])
    physical=(p-o)@n;station=np.ones(len(p));owned=np.unique(f[region])
    station[owned]=np.clip((physical[owned]-start)/(end-start),0.,1.)
    profiles=[]
    for level in levels:
        coordinate=float(start+level*(end-start))
        mask=region & np.all(physical[f]<=coordinate+tolerance_m,axis=1)
        partial=f[mask];rings=boundary_loops(partial)
        edges=np.unique(np.sort(np.vstack([partial[:,[0,1]],partial[:,[1,2]],partial[:,[2,0]]]),axis=1),axis=0)
        chi=len(np.unique(partial))-len(edges)+len(partial)
        if face_components(partial)[0]!=1 or len(rings)!=2 or chi!=0:
            raise ValueError('A collar station crosses a branch or an unrelated component')
        if any(np.max(abs(physical[r]-coordinate))>tolerance_m for r in rings):
            raise ValueError('A station boundary is not a common physical section')
        for ring in rings:station[ring]=float(level)
        profiles.append(dict(station=float(level),coordinate_m=coordinate,rings=rings,
            areas_m2=[projected_area(p[r])for r in rings]))
    return dict(vertices=p,faces=f,station=station,root_region=region,
        parent_faces=parents,corner_barycentrics=bary,profiles=profiles,
        refinement_steps=steps,normal=n,origin=o)


def planar_open_wall_volume(vertices, faces, rings, *, origin, tolerance_m=1e-10):
    """Oriented wall volume with an analytic closure on one shared plane.

    Oppositely oriented outer/inner cycles subtract. No fluid-mouth disk or
    triangulated closure is returned, displayed, or installed.
    """
    p=np.asarray(vertices,float);f=np.asarray(faces,int);o=np.asarray(origin,float)
    if len(rings)!=2:raise ValueError('One paired annular outlet is required')
    actual=boundary_loops(f)
    if {frozenset(r)for r in actual}!={frozenset(r)for r in rings}:
        raise ValueError('Declared profiles must be the actual wall boundary')
    allp=p[np.concatenate(actual)]-o;_,_,vh=np.linalg.svd(allp-allp.mean(0),full_matrices=False)
    if np.max(abs((allp-allp.mean(0))@vh[-1]))>tolerance_m:
        raise ValueError('An analytic closure requires coplanar inner and outer profiles')
    vectors=np.array([np.cross(p[r]-o,p[np.roll(r,-1)]-o).sum(0)*.5 for r in actual])
    if np.dot(vectors[0],vectors[1])>=0:
        raise ValueError('The inner and outer boundary orientations must oppose')
    cap_area_vector=-vectors.sum(0);tri=p[f]-o
    side=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6
    cap=np.dot(allp[0],cap_area_vector)/3
    return float(side+cap)
