"""Orient a tetrahedral complex from its actual boundary and shared faces.

    Metric determinants alone cannot orient a flat interface transition:
    cancellation can report either sign for an exactly coplanar tetrahedron.
    Such a transition is returned explicitly, never called occupied tissue.
"""
from collections import deque
from fractions import Fraction
import numpy as np


def _parity(faces):
    f=np.asarray(faces,int)
    n=(f[:,0]>f[:,1]).astype(int)+(f[:,0]>f[:,2]).astype(int)+(f[:,1]>f[:,2]).astype(int)
    return 1-2*(n%2)


def _exact_determinant(q):
    p=[[Fraction(float(x))for x in row]for row in q]
    a,b,c=[[x-y for x,y in zip(row,p[0])]for row in p[1:]]
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


def orient_complex(points,tetrahedra,boundary):
    p=np.asarray(points,float);t=np.asarray(tetrahedra,int).copy();boundary=np.asarray(boundary,int)
    faces=np.concatenate([t[:,r]for r in ((1,2,3),(0,3,2),(0,1,3),(0,2,1))])
    owners=np.tile(np.arange(len(t)),4);parity=_parity(faces)
    keys,inverse,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
    if np.any(count>2):
        raise ValueError('The tetrahedral complex has a nonmanifold face')
    order=np.argsort(inverse,kind='stable');start=np.r_[0,np.cumsum(count)[:-1]]
    boundary_keys={tuple(sorted(f)):s for f,s in zip(boundary,_parity(boundary))}
    singles=order[start[count==1]]
    if {tuple(sorted(faces[i]))for i in singles}!=set(boundary_keys):
        raise ValueError('The tetrahedral complex changed the literal boundary facets')
    directions=np.zeros(len(t),np.int8);queue=deque();adjacency=[[]for _ in t]
    def assign(i,value):
        if directions[i] and directions[i]!=value:
            raise ValueError('The tetrahedral complex has inconsistent oriented adjacency')
        if not directions[i]:
            directions[i]=value;queue.append(int(i))
    for i in singles:
        assign(owners[i],int(boundary_keys[tuple(sorted(faces[i]))]*parity[i]))
    for index in start[count==2]:
        a,b=order[index:index+2];i,j=owners[a],owners[b];factor=int(-parity[a]*parity[b])
        adjacency[i].append((int(j),factor));adjacency[j].append((int(i),factor))
    while queue:
        i=queue.popleft()
        for j,factor in adjacency[i]:
            assign(j,int(directions[i])*factor)
    if np.any(directions==0):
        raise ValueError('A tetrahedral component is disconnected from the actual boundary')
    changed=directions<0;t[changed]=t[changed][:,[0,2,1,3]]
    q=p[t];a,b,c=(q[:,k]-q[:,0]for k in (1,2,3))
    det=np.einsum('ij,ij->i',a,np.cross(b,c))
    error=128*np.finfo(float).eps*np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)*np.linalg.norm(c,axis=1)
    if np.any(det < -error):
        i=int(np.flatnonzero(det < -error)[0]);failure=ValueError('A nonzero construction tetrahedron is geometrically inverted')
        failure.details=dict(tetrahedron=i,controls=t[i].tolist(),points_m=q[i].tolist(),determinant_m3=float(det[i]),
                             exact_predicate=False,positions_modified=False)
        raise failure
    flat=[];exact=0
    for i in np.flatnonzero(np.abs(det)<=error):
        value=_exact_determinant(q[i]);exact+=1
        if value<0:
            failure=ValueError('An exact nonzero construction tetrahedron is geometrically inverted')
            failure.details=dict(tetrahedron=int(i),controls=t[i].tolist(),points_m=q[i].tolist(),
                determinant_m3=float(value),exact_predicate=True,positions_modified=False)
            raise failure
        if value==0:
            flat.append(int(i))
    return t,dict(orientations_changed=int(changed.sum()),exact_determinants_checked=exact,
        flat_interface_glue_tetrahedra=flat,nonzero_tetrahedra_positive=True,
        literal_boundary_facets_preserved=True,all_shared_faces_oppositely_oriented=True,
        flat_transitions_require_homogeneous_tissue_ownership=True,positions_modified=False)
