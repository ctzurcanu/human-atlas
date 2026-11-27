"""Literal uncapped cross-section of a declared closed fluid compartment.

This helper selects the simple contour owning a query point. It explicitly
rejects hole-bearing tissue sections; callers must retain their separate
contour ownership rather than filling a lumen with an inferred disk.
"""
import numpy as np
from .scalar_surface_clip import clip_scalar_surface
from .annular_pipe_envelopes import boundary_loops
from geometry_tools.audit import topology


def polygon_contains(poly,point):
    a=np.asarray(poly,float)-point;b=np.roll(a,-1,axis=0)
    angle=np.arctan2(a[:,0]*b[:,1]-a[:,1]*b[:,0],np.einsum('ij,ij->i',a,b)).sum()
    return bool(abs(angle)>np.pi)


def polygon_clearance(poly,point):
    a=np.asarray(poly,float);b=np.roll(a,-1,axis=0);edge=b-a;den=np.einsum('ij,ij->i',edge,edge)
    # A literal repeated point contributes a zero-length segment to the same
    # planar polygon. Omit only exactly zero projected edges, never a short
    # anatomical edge or a near-coincident vertex under a distance tolerance.
    keep=den>0
    if not np.any(keep):raise ValueError('A section has no nonzero projected edges')
    a,edge,den=a[keep],edge[keep],den[keep]
    t=np.clip(np.einsum('ij,ij->i',point-a,edge)/den,0,1)
    return float(np.linalg.norm(a+t[:,None]*edge-point,axis=1).min())


def simple_compartment_section(mesh,point,normal,*,classification_tolerance_m=2e-13):
    point=np.asarray(point,float);normal=np.asarray(normal,float)
    if point.shape!=(3,)or normal.shape!=(3,)or not np.isfinite(np.r_[point,normal]).all()or np.linalg.norm(normal)==0:
        raise ValueError('A finite point and nonzero physical plane normal are required')
    if not topology(mesh)['closed_oriented_manifold']:
        raise ValueError('An open surface cannot define the requested fluid compartment')
    if not np.isfinite(classification_tolerance_m)or classification_tolerance_m<0:
        raise ValueError('A finite nonnegative plane-classification bound is required')
    normal=normal/np.linalg.norm(normal);axis=np.eye(3)[np.argmin(abs(normal))]
    u=np.cross(normal,axis);u/=np.linalg.norm(u);v=np.cross(normal,u);basis=np.c_[u,v]
    scalar=(mesh.vertices-point)@normal
    near=abs(scalar)<=classification_tolerance_m
    if np.any(np.all(near[mesh.faces],axis=1)):
        raise ValueError('The requested section coincides with a source face')
    # Apply the caller's declared plane-classification bound to the field,
    # never to source coordinates. A roundoff-sized opposite sign otherwise
    # produces an edge fraction rounded to one and an unrepresentable sliver.
    # This neither welds vertices nor broadens the existing error bound.
    scalar[near]=0.
    cut=clip_scalar_surface(mesh.vertices,mesh.faces,scalar,0.);loops=boundary_loops(cut['negative_faces']);rows=[]
    for ring in loops:
        world=cut['positions'][ring];xy=(world-point)@basis;cross=xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1);area=.5*cross.sum()
        if area==0:raise ValueError('A section contour has zero projected area')
        centroid=((xy+np.roll(xy,-1,axis=0))*cross[:,None]).sum(0)/(6*area)
        rows.append(dict(world_points=world,plane_points=xy,area_m2=abs(area),signed_area_m2=area,centroid_plane=centroid))
    enclosing=[i for i,row in enumerate(rows)if polygon_contains(row['plane_points'],np.zeros(2))]
    if not enclosing:raise ValueError('No section component contains the query point')
    owner=min(enclosing,key=lambda i:rows[i]['area_m2']);row=rows[owner]
    nested=[i for i,other in enumerate(rows)if i!=owner and polygon_contains(row['plane_points'],other['plane_points'][0])]
    if len(enclosing)>1 or nested:
        raise ValueError('Hole-bearing sections require explicit multi-contour ownership')
    return {**row,'point_m':point,'normal':normal,'basis':basis,
        'centroid_m':point+basis@row['centroid_plane'],
        'centroid_inside_contour':polygon_contains(row['plane_points'],row['centroid_plane']),
        'query_clearance_m':polygon_clearance(row['plane_points'],np.zeros(2)),
        'centroid_clearance_m':polygon_clearance(row['plane_points'],row['centroid_plane']),
        'other_disjoint_contours':len(rows)-1,'caps_created':0,'source_geometry_changed':False,
        'plane_classification_tolerance_m':float(classification_tolerance_m),
        'controls_classified_on_plane':int(np.count_nonzero(near)),
        'anatomical_acceptance':False}
