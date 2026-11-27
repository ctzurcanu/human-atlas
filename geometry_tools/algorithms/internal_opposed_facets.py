"""Cancel only opposed source facets lying inside a closed regular envelope.

This is a source-envelope union operation. It never treats opposed external
walls as redundant and does not infer blood-wall thickness or clinical flow.
All surviving source coordinates, winding and appearance words remain literal.
"""
import numpy as np
from geometry_tools.contact import pair_contacts
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.spatial import BVH

def internal_opposed_facets(mesh):
    import igl
    physical,alias=np.unique(mesh.vertices,axis=0,return_inverse=True)
    f=alias[mesh.faces]
    _,group,count=np.unique(np.sort(f,axis=1),axis=0,return_inverse=True,return_counts=True)
    pairs=[]
    for k in np.flatnonzero(count==2):
        ids=np.flatnonzero(group==k);a,b=ids
        if any(np.array_equal(f[a],np.roll(f[b][::-1],j)) for j in range(3)):
            pairs.append(ids.tolist())
    if not pairs:raise ValueError('No exactly opposed paired source facets')
    removed=np.asarray(pairs,int).ravel();keep=np.setdiff1d(np.arange(len(f)),removed)
    candidate=Mesh(mesh.vertices,mesh.faces[keep],mesh.frame,mesh.attributes,mesh.source)
    cf=f[keep];edges=cf[:,[[0,1],[1,2],[2,0]]].reshape(-1,2)
    ue,inv,inc=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    if not np.all(inc==2):raise ValueError('Retained carrier must be closed with exactly two faces per edge')
    signs=np.where(edges[:,0]<edges[:,1],1,-1)
    if np.any(np.bincount(inv,weights=signs)!=0):raise ValueError('Retained envelope winding is inconsistent')
    diag=diagnostics(candidate,include_pairs=True)
    if diag['nonmanifold_edges'] or diag['nonmanifold_vertices'] or diag['self_intersection_pairs'] or not diag['orientable']:
        raise ValueError('Retained closed envelope is not regular and self-disjoint')
    old_triangles=mesh.triangles
    triangles=old_triangles[np.asarray(pairs)[:,0]]
    normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
    centers=triangles.mean(1);tree=BVH(candidate)
    distance=tree.closest_points(centers)[0]['distance_m']
    if np.any(distance<1e-7):raise ValueError('An opposed facet touches the retained exterior at its center')
    epsilon=np.minimum(distance*.1,1e-5)
    points=np.concatenate([centers,centers+normals*epsilon[:,None],centers-normals*epsilon[:,None]])
    origin=mesh.vertices.mean(0)
    winding=igl.winding_number(candidate.vertices-origin,candidate.faces,points-origin)
    # On a removed source face, the original winding query is undefined.
    # Compare original/candidate values only at the two off-face probes.
    old_winding=igl.winding_number(mesh.vertices-origin,mesh.faces,points[len(centers):]-origin)
    if np.max(np.abs(np.abs(winding)-1))>1e-9:raise ValueError('Opposed facets are not inside the retained solid envelope')
    if np.max(np.abs(winding[len(centers):]-old_winding))>1e-9:raise ValueError('Source oriented volume changed')
    deleted=Mesh(mesh.vertices,mesh.faces[removed],mesh.frame)
    arrays,contact=pair_contacts(deleted,candidate,0.,max_tests=250000)
    if not contact['complete_within_pair_budget'] or contact['crossing_face_pairs'] or contact['coplanar_overlap_face_pairs']:
        raise ValueError('Removed facets cross or overlap the retained exterior')
    # Opposite exact triangles cancel as an oriented chain, not only in volume.
    for a,b in pairs:
        assert any(np.array_equal(old_triangles[a],np.roll(old_triangles[b][::-1],j,axis=0)) for j in range(3))
    return keep,dict(complete=True,removedSourceFaceIds=removed.tolist(),retainedSourceFaceIds=keep.tolist(),opposedSourcePairs=pairs,
        everyRemovedPairCancelsAsAnOrientedChain=True,retainedEnvelopeClosed=True,retainedEnvelopeConsistentlyOriented=True,
        retainedEnvelopeNonmanifoldEdges=0,retainedEnvelopeNonmanifoldVertices=0,retainedEnvelopeSelfPairs=0,
        minimumRemovedCentroidDistanceToExterior_m=float(distance.min()),
        maximumInteriorWindingError=float(np.abs(np.abs(winding)-1).max()),
        maximumSourceWindingDifference=float(np.abs(winding[len(centers):]-old_winding).max()),
        windingProbePoints_m=points.tolist(),retainedEnvelopeWinding=winding.tolist(),
        removedVsExteriorContact=contact,sourceCoordinatesAndRetainedFieldsLiteral=True,
        anatomicalAcceptance=False,finiteWallOrClinicalPatencyEstablished=False)

def replay_proof_matches(first,second):
    """Native winding summation may differ by a few final floating-point bits."""
    numerical={'maximumInteriorWindingError','maximumSourceWindingDifference','retainedEnvelopeWinding'}
    assert {k:v for k,v in first.items() if k not in numerical}=={k:v for k,v in second.items() if k not in numerical}
    for key in numerical:np.testing.assert_allclose(first[key],second[key],rtol=0,atol=1e-12)
