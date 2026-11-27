"""Resolve redundant same-oriented draw faces with conflicting normal words.

Geometry and every retained source field remain literal. An existing variant
is selected by agreement with unique area-weighted physical face normals.
Other attribute conflicts, opposed sheets and inverted variants are retained.
This changes ambiguous duplicate shading; it is not an appearance identity.
"""
import numpy as np


def resolve_normal_variants(mesh):
    _,physical=np.unique(mesh.vertices,axis=0,return_inverse=True)
    faces=physical[mesh.faces]
    unique,first,group,counts=np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_inverse=True,return_counts=True)
    triangles=mesh.triangles
    fn=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    summed=np.zeros((int(physical.max())+1,3))
    for corner in range(3):np.add.at(summed,faces[first,corner],fn[first])
    lengths=np.linalg.norm(summed,axis=1)
    reference=np.divide(summed,lengths[:,None],out=np.zeros_like(summed),where=lengths[:,None]>0)
    normal=np.asarray(mesh.attributes['normals'],float)*np.asarray(mesh.attributes['normalSigns'])[:,None]
    normal_length=np.linalg.norm(normal,axis=1)
    if not np.isfinite(normal).all() or np.any(normal_length==0):
        raise ValueError('Finite nonzero source normal words required')
    normal/=normal_length[:,None]
    keep=np.ones(len(faces),bool);resolved=[];retained=[]
    for gi in np.flatnonzero(counts>1):
        ids=np.flatnonzero(group==gi);orders=[];base=int(ids[0]);reason=None
        for face in ids:
            match=next((np.roll(np.arange(3),k) for k in range(3) if np.array_equal(faces[face,np.roll(np.arange(3),k)],faces[base])),None)
            if match is None:reason='Opposed winding retained';break
            orders.append(match)
        if reason is None:
            for face,order in zip(ids,orders):
                x,y=mesh.faces[base],mesh.faces[face,order]
                if any(not np.array_equal(v[x],v[y]) for key,v in mesh.attributes.items() if key not in ['normals','normalSigns']):reason='Other corner field conflict retained';break
                face_length=np.linalg.norm(fn[face])
                if face_length==0 or np.any(normal[y]@(fn[face]/face_length)<=0):reason='Degenerate or inward source normal variant retained';break
            if np.any(lengths[faces[base]]==0):reason='Unresolved physical normal star retained'
        if reason is not None:
            retained.append(dict(sourceFaces=ids.tolist(),reason=reason));continue
        scores=[]
        for face,order in zip(ids,orders):
            dots=np.einsum('ij,ij->i',normal[mesh.faces[face,order]],reference[faces[base]])
            scores.append((float(dots.min()),float(dots.sum()),-int(face)))
        selected=int(max(range(len(ids)),key=lambda i:scores[i]));chosen=int(ids[selected])
        removed=[int(face) for face in ids if face!=chosen]
        keep[removed]=False
        max_angle=max(float(np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',normal[mesh.faces[chosen,orders[selected]]],normal[mesh.faces[face,order]]),-1,1))).max()) for face,order in zip(ids,orders))
        resolved.append(dict(sourceFaces=ids.tolist(),retainedSourceFace=chosen,removedSourceFaces=removed,sourceCornerOrders=[x.tolist() for x in orders],agreementScores=[list(s) for s in scores],maximumSourceNormalVariantAngleDegrees=max_angle))
    kept=np.flatnonzero(keep)
    assert len(kept)==len(np.unique(kept))
    assert set(map(tuple,np.sort(faces[kept],axis=1)))==set(map(tuple,unique))
    return kept,dict(complete=True,retainedSourceFaceIds=kept.tolist(),resolvedGroups=resolved,retainedConflicts=retained,duplicateCopiesRemoved=int((~keep).sum()),sourceSurfaceSupportUnchanged=True,positionsAndAllRetainedSourceFieldsLiteral=True,normalWordsReconstructed=False,appearanceIdenticalClaim=False,opposedSheetsMerged=False,anatomicalAcceptance=False)
