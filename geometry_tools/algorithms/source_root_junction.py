"""Join declared original branch collars with an explicit majority root chart.

Actual root incidence is required before the bounded internal chart correction.
This is not a nearest-mesh weld or a clinical attachment classifier.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
from geometry_tools.model import Mesh
from .graft_collar_union import union_graft_collars


def root_area(v):
    q=v-v.mean(0)
    return float(np.linalg.norm(np.cross(q,np.roll(q,-1,axis=0)).sum(0)*.5))


def union_source_junction(mesh,cells,roots,canonical,*,branch_count,maximum_correction_m,
                      maximum_root_area_relative_change):
    roots=[np.asarray(x,int)for x in roots];canonical=np.asarray(canonical,int)
    if branch_count not in (3,4)or len(cells)!=branch_count or len(roots)!=branch_count or any(len(x)!=6 for x in roots):
        raise ValueError('Explicit three/four-branch six-control annular roots required')
    reference=mesh.vertices[canonical];rootset=frozenset(map(tuple,reference))
    if sum(frozenset(map(tuple,mesh.vertices[x]))==rootset for x in roots)<2:
        raise ValueError('Two literal complete roots must establish the source chart')
    v=mesh.vertices.tolist();attrs={k:x.tolist()for k,x in mesh.attributes.items()}
    faces=mesh.faces.copy();edits=[];areas=[]
    for ci,root in enumerate(roots):
        w=mesh.vertices[root];distance=np.linalg.norm(w[:,None]-reference[None],axis=2)
        first,last=linear_sum_assignment(distance)
        if np.max(distance[first,last])>maximum_correction_m:
            raise ValueError('Source-root correction exceeds explicit construction bound')
        match=canonical[last];shared=distance[first,last]==0
        if np.count_nonzero(shared)<2:
            raise ValueError('An actual shared source edge is required')
        common={tuple(x)for x in w[shared]}
        if not any(tuple(a)in common and tuple(b)in common
                   for a,b in zip(w,np.roll(w,-1,axis=0))):
            raise ValueError('Shared source points do not form a literal root edge')
        if not any(np.array_equal(np.roll(match,k),canonical)or
                   np.array_equal(np.roll(match[::-1],k),canonical)for k in range(6)):
            raise ValueError('Source-root correspondence changes cyclic order')
        area=root_area(w);new_area=root_area(mesh.vertices[match])
        relative=abs(new_area/area-1)
        if relative>maximum_root_area_relative_change:
            raise ValueError('Root chart correction exceeds declared area bound')
        areas.append(dict(cell=ci,sourceArea_m2=area,canonicalArea_m2=new_area,
                          relativeAreaChange=relative,exactSharedControls=int(shared.sum())))
        for old,new,delta in zip(root,match,distance[first,last]):
            if delta==0:continue
            # A root coordinate may have multiple source corner-field aliases.
            # Remap every literal alias used by this cell, retaining each field.
            local=faces[cells[ci]].copy();used=np.unique(mesh.faces[cells[ci]])
            source_aliases=used[np.all(mesh.vertices[used]==mesh.vertices[old],axis=1)]
            if not len(source_aliases):
                raise ValueError('Declared root coordinate is absent from the source cell')
            for original in source_aliases:
                alias=len(v);v.append(mesh.vertices[new].tolist())
                for key in attrs:attrs[key].append(mesh.attributes[key][original].tolist())
                local[local==original]=alias
                edits.append(dict(cell=ci,oldRawId=int(original),canonicalRawId=int(new),
                                  newAliasRawId=alias,displacement_m=float(delta)))
            faces[cells[ci]]=local
    working=Mesh(np.array(v),faces,mesh.frame,
                 {k:np.array(x,dtype=mesh.attributes[k].dtype)for k,x in attrs.items()},mesh.source)
    result,proof=union_graft_collars(working,cells,[canonical.tolist()],junction_count=1)
    proof.update(rootSeamEdits=edits,rootAreas=areas,
                 maximumRootSeamCorrection_m=max((x['displacement_m']for x in edits),default=0),
                 originalCoordinateAndFieldPrefixLiteral=True,
                 sourceRootIncidenceRequired=True,clinicalAttachmentEstablished=False)
    return result,proof
