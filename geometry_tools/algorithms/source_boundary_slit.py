"""Zip an explicitly classified narrow source crack without adding a cap.

This is not automatic hole filling: callers must inventory the defect, identify
the two corresponding shores and exclude a genuine adjoining vessel mouth.
"""
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics
from .tube_continuity import boundary_loops
from .surface_normals import transport_vertex_words

def _cycle(points):
    p=[tuple(x) for x in points]
    return min(tuple(p[k:]+p[:k]) for k in range(len(p)))

def zip_source_slit(mesh,first_arc,second_arc,*,maximum_motion_m,
                    maximum_singular_strain=.04,maximum_facet_turn_degrees=3.):
    a,b=[np.asarray(x,int) for x in [first_arc,second_arc]]
    if a.ndim!=1 or a.shape!=b.shape or len(a)<4:
        raise ValueError('Two explicit corresponding boundary shores are required')
    if not 0<maximum_motion_m<=1e-4:
        raise ValueError('Explicit source crack displacement at most 0.1 mm required')
    v=mesh.vertices;f=mesh.faces
    if not np.array_equal(v[a[[0,-1]]],v[b[[0,-1]]]):
        raise ValueError('Both source crack endpoints must already be identical')
    gap=np.linalg.norm(v[b]-v[a],axis=1)
    span=np.linalg.norm(v[a[-1]]-v[a[0]])
    if gap.max()>maximum_motion_m or gap.max()>=.02*span:
        raise ValueError('Declared shores are not one sufficiently narrow source crack')
    physical,alias=np.unique(v,axis=0,return_inverse=True)
    # Only for reading the original intrinsic boundary, keep the declared raw
    # second-shore controls distinct at an already collapsed pinch.
    logical=alias.copy();positions=physical.tolist()
    for raw in b[1:-1]:
        if raw in a:raise ValueError('Interior shore controls must be distinct source controls')
        if alias[raw] in alias[a]:
            logical[raw]=len(positions);positions.append(v[raw].tolist())
    loops=boundary_loops(logical[f]);target=np.r_[logical[a],logical[b][-2:0:-1]]
    matches=[i for i,loop in enumerate(loops) if set(loop)==set(target)]
    if len(matches)!=1:raise ValueError('The two shores must bound exactly one complete source crack')
    loop=loops[matches[0]]
    if not any(np.array_equal(loop,np.roll(order,k)) for order in [target,target[::-1]] for k in range(len(target))):
        raise ValueError('Source crack winding is inconsistent with the declared shores')
    positions=np.asarray(positions)
    genuine=[_cycle(positions[loop]) for i,loop in enumerate(loops) if i!=matches[0]]
    q=v.copy();moves=[]
    for first,second in zip(a[1:-1],b[1:-1]):
        if np.array_equal(v[first],v[second]):continue
        ids=np.flatnonzero(np.all(v==v[second],axis=1))
        if np.any(np.isin(ids,a)):raise ValueError('Moving alias belongs to the fixed shore')
        q[ids]=v[first]
        moves.append(dict(rawIds=ids.tolist(),fixedSourceControl=int(first),
            original_m=v[second].tolist(),target_m=v[first].tolist(),distance_m=float(np.linalg.norm(v[second]-v[first]))))
    if not moves:raise ValueError('No source crack correction required')
    triangles=v[f];E=(triangles[:,1:]-triangles[:,:1]).transpose(0,2,1)
    _,s,V=np.linalg.svd(E,full_matrices=False)
    if np.any(s<=0):raise ValueError('Source has degenerate facets')
    R=np.swapaxes(V,1,2)/s[:,None,:]
    new=(q[f][:,1:]-q[f][:,:1]).transpose(0,2,1)
    strain=np.abs(np.linalg.svd(new@R,compute_uv=False)-1).max(1)
    if np.any(strain>maximum_singular_strain):raise ValueError('Source facet metric budget exceeded')
    normals,transport=transport_vertex_words(v,q,f,mesh.attributes['normals'],maximum_facet_rotation_degrees=maximum_facet_turn_degrees)
    attrs=dict(mesh.attributes);attrs['normals']=normals
    candidate=Mesh(q,f,mesh.frame,attrs,mesh.source)
    d=diagnostics(candidate,include_pairs=True)
    if d['nonmanifold_edges'] or d['nonmanifold_vertices'] or d['self_intersection_pairs'] or not d['orientable']:
        raise ValueError('Zipped whole source must be manifold, orientable and self-disjoint')
    points,ids=np.unique(q,axis=0,return_inverse=True);newloops=boundary_loops(ids[f])
    if sorted(_cycle(points[loop]) for loop in newloops)!=sorted(genuine):
        raise ValueError('A genuine source mouth changed')
    return candidate,dict(complete=True,firstArc=a.tolist(),secondArc=b.tolist(),moves=moves,
        maximumMotion_m=max(x['distance_m'] for x in moves),maximumSingularStrain=float(strain.max()),
        maximumSingularStrainAllowed=maximum_singular_strain,normalTransport=transport,
        maximumFacetTurnAllowed_degrees=maximum_facet_turn_degrees,
        sourceCoordinatesElsewhereLiteral=True,allSourceFacesRetained=True,
        genuineSourceMouthsRetained=len(genuine),removedDefectiveSlitBoundaries=1,
        facesCreated=0,capsCreated=0,anatomicalAcceptance=False)
