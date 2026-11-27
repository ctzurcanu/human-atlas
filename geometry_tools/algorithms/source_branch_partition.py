"""Compact disjoint source-face partitions; no reconstructed or overlapping faces."""
import numpy as np
from geometry_tools.model import Mesh

def partition_source_faces(mesh,groups):
    owner=np.full(len(mesh.faces),-1,int);outputs=[]
    for k,ids in enumerate(groups):
        ids=np.asarray(ids,int)
        if not len(ids)or ids.min()<0 or ids.max()>=len(owner)or len(np.unique(ids))!=len(ids)or np.any(owner[ids]!=-1):raise ValueError('Distinct nonempty original source faces required')
        owner[ids]=k
    all_groups=[np.flatnonzero(owner==-1),*[np.asarray(g,int)for g in groups]]
    if not len(all_groups[0]):raise ValueError('The parent vessel must retain a nonempty trunk')
    for ids in all_groups:
        used,inv=np.unique(mesh.faces[ids],return_inverse=True)
        attrs={}
        for key,value in mesh.attributes.items():
            if len(value)!=len(mesh.vertices):raise ValueError('Explicit per-vertex source attributes required')
            attrs[key]=value[used].copy()
        candidate=Mesh(mesh.vertices[used],inv.reshape(-1,3),mesh.frame,attrs,mesh.source)
        assert np.array_equal(candidate.triangles,mesh.triangles[ids])
        outputs.append((candidate,ids,used))
    assert np.array_equal(np.sort(np.concatenate([x[1]for x in outputs])),np.arange(len(mesh.faces)))
    return outputs
