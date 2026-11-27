"""Disjoint rigid closures of literal profiles and unrecovered source charts.

No branch identity is inferred. Every actual control, exact coordinate alias,
annular section and unrecovered chart retains one owner and its complete shape.
"""
import numpy as np
from .source_tube_topology import recover_regular_source_annuli

def rigid_source_profile_groups(mesh,*,maximum_core_patch_m=None):
    rec=recover_regular_source_annuli(mesh.vertices,mesh.faces,np.arange(len(mesh.faces)//2*2).reshape(-1,2))
    alias=np.asarray(rec['raw_to_physical']);parent=np.arange(len(rec['physical_points']))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def merge(ids):
        ids=np.unique(ids);r=root(int(ids[0]))
        for i in ids:parent[root(int(i))]=r
    if maximum_core_patch_m is None:
        for ring in rec['rings']:merge(ring)
        for faces in rec['unrecovered_source_face_components']:merge(alias[np.unique(mesh.faces[faces])])
        physical=np.array([root(i)for i in alias]);_,group=np.unique(physical,return_inverse=True)
    else:
        from .source_profile_axis import source_profile_graph,complete_source_sections
        from .source_control_clearance import source_control_basis
        guide=source_profile_graph(mesh.vertices,mesh.faces,triangle_pairing=np.arange(len(mesh.faces)//2*2).reshape(-1,2))
        sections=complete_source_sections(guide,alias)
        W,_,_,_=source_control_basis(mesh.vertices,mesh.faces,sections,alias,core_patch_diameter_m=maximum_core_patch_m)
        assert np.all((W==0)|(W==1))and np.all(W.sum(1)==1)
        group=W.argmax(1)
    groups=[np.flatnonzero(group==i)for i in range(int(group.max())+1)]
    edges=np.sort(alias[np.concatenate([mesh.faces[:,[0,1]],mesh.faces[:,[1,2]],mesh.faces[:,[2,0]]])],axis=1)
    edges,count=np.unique(edges,axis=0,return_counts=True);boundary_physical=np.unique(edges[count==1])
    boundary_raw=np.flatnonzero(np.isin(alias,boundary_physical));fixed=np.unique(group[boundary_raw])
    adjacency=np.zeros((len(groups),len(groups)))
    for face in group[mesh.faces]:
        for i,j in zip(face,np.roll(face,-1)):
            if i!=j:adjacency[i,j]=adjacency[j,i]=1
    assert np.isfinite(mesh.vertices).all() and all(len(g)for g in groups)
    return dict(topology=rec,groups=groups,group_of_vertex=group,adjacency=adjacency,
                boundary_raw=boundary_raw,boundary_groups=fixed,
                sourceFaceOwnershipComplete=True,anatomicalIdentityInferred=False)
