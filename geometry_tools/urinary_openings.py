"""Bind an empty anatomical opening to its actual finite mucosal surround.

The curve is metadata over an already verified shared wall. This creates no
solid, disk, overlay mesh, texture or duplicate tissue surface.
"""
import numpy as np
from .model import Mesh
from .spatial import BVH


def bind_opening_boundary(name, curve, supports, plan_sha, *, maximum_offset_m=1e-7):
    curve=np.asarray(curve,float)
    if curve.ndim!=2 or curve.shape[1]!=3 or len(curve)<3 or not np.isfinite(curve).all():
        raise ValueError('An explicit finite closed clinical rim is required')
    fractions=np.arange(4)/4
    samples=(curve[:,None]+fractions[None,:,None]*(np.roll(curve,-1,axis=0)-curve)[:,None]).reshape(-1,3)
    distance=np.full(len(samples),np.inf);points=np.full_like(samples,np.nan)
    faces=np.full(len(samples),-1,int);owners=np.full(len(samples),'',object);bary=np.full((len(samples),3),np.nan)
    for pid,(mesh,record) in supports.items():
        groups=[np.arange(g['start']//3,(g['start']+g['count'])//3)
            for g in record['tissueDrawGroups'] if not g.get('neighbor')]
        free=np.concatenate(groups)if groups else np.empty(0,int)
        if not len(free):raise ValueError('A named ostial surround must touch the physical urinary lining: '+pid)
        query,_=BVH(Mesh(mesh.vertices,mesh.faces[free],mesh.frame)).closest_points(samples)
        take=query['distance_m']<distance
        distance[take]=query['distance_m'][take];points[take]=query['points'][take]
        faces[take]=free[query['face_ids'][take]];owners[take]=pid;bary[take]=query['barycentric'][take]
    if not supports or not np.isfinite(distance).all() or distance.max()>maximum_offset_m:
        raise ValueError('The named mucosal tissue does not contain the complete clinical opening curve')
    return dict(name=name,representation='wall-rim',partIds=list(supports),
        points=points.tolist(),samplePartIds=owners.tolist(),faceIds=faces.tolist(),barycentric=bary.tolist(),
        sourcePlanSha256=plan_sha,maximumSourceOffsetM=float(distance.max()),
        sourceCurveControls=len(curve),samplesPerSourceEdge=4)


def shared_lumen_rim(first, first_record, second, second_record):
    """Extract the real lumen edge of one reciprocal annular tissue join.

    Coordinates are literal, including corner-normal duplicates. Both mucosal
    free surfaces must meet the same single closed loop on the shared annulus.
    No nearest-point tolerance, invented disk or new geometry is introduced.
    """
    from collections import Counter, defaultdict
    if first.frame != second.frame:
        raise ValueError('Opening peers must occupy the same declared frame')
    records = (first_record, second_record)
    if not all(isinstance(r.get('id'), str) for r in records):
        raise ValueError('Both tissue records need explicit identities')
    meshes = (first, second)
    points, inverse = np.unique(np.concatenate([m.vertices for m in meshes]), axis=0, return_inverse=True)
    offsets = np.cumsum([0] + [len(m.vertices) for m in meshes])
    seams = []
    rims = []
    def boundary(faces):
        count = Counter(tuple(sorted(edge)) for face in faces
                        for edge in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])))
        if any(n > 2 for n in count.values()):
            raise ValueError('Nonmanifold tissue patch cannot establish a lumen rim')
        return {edge for edge, n in count.items() if n == 1}
    def signature(mesh, faces, reverse=False, normals=False):
        result = []
        words = mesh.attributes.get('normals')
        if normals:
            words = words.astype(np.int64)
            if 'normalSigns' in mesh.attributes:
                words = words * mesh.attributes['normalSigns'][:, None]
        for face in faces:
            if reverse:
                face = face[[0, 2, 1]]
            corners = [(tuple(mesh.vertices[v]), tuple((-1 if reverse else 1) * words[v]))
                       if normals else tuple(mesh.vertices[v]) for v in face]
            result.append(min(tuple(corners[k:] + corners[:k]) for k in range(3)))
        return Counter(result)
    for i, mesh in enumerate(meshes):
        groups = records[i].get('tissueDrawGroups', [])
        for group in groups:
            if (type(group.get('start')) is not int or type(group.get('count')) is not int
                or group['start'] < 0 or group['count'] <= 0 or group['start'] % 3 or group['count'] % 3
                or group['start'] + group['count'] > len(mesh.faces) * 3):
                raise ValueError('Invalid tissue draw range')
        seam = [k for g in groups if g.get('neighbor') == records[1-i]['id']
                for k in range(g['start']//3, (g['start']+g['count'])//3)]
        free = [k for g in groups if not g.get('neighbor')
                for k in range(g['start']//3, (g['start']+g['count'])//3)]
        if not seam or not free:
            raise ValueError('A lumen rim needs a declared join and real free mucosa on both peers')
        welded = inverse[offsets[i]:offsets[i+1]][mesh.faces]
        rims.append(boundary(welded[seam]) & boundary(welded[free]))
        seams.append(np.asarray(seam, int))
    if signature(first, first.faces[seams[0]]) != signature(second, second.faces[seams[1]], True):
        raise ValueError('Opening peers do not have literal opposed annular triangles')
    normal_presence = ['normals' in m.attributes for m in meshes]
    if any(normal_presence) and not all(normal_presence):
        raise ValueError('Only one opening peer has a render-normal field')
    if all(normal_presence) and signature(first, first.faces[seams[0]], normals=True) != signature(second, second.faces[seams[1]], True, True):
        raise ValueError('Opening peers do not have literal opposed normal words')
    if not rims[0] or rims[0] != rims[1]:
        raise ValueError('Both mucosal surfaces must meet the identical annular lumen edge')
    adjacency = defaultdict(list)
    for u, v in rims[0]:
        adjacency[u].append(v); adjacency[v].append(u)
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise ValueError('The lumen edge must have degree two at every point')
    walk = [min(adjacency)]
    previous = None
    while True:
        nxt = min(v for v in adjacency[walk[-1]] if v != previous)
        if nxt == walk[0]:
            break
        if nxt in walk:
            raise ValueError('The lumen edge revisits a point before closing')
        previous = walk[-1]; walk.append(nxt)
    if len(walk) != len(adjacency):
        raise ValueError('More than one lumen edge exists; inner and outer surfaces are ambiguous')
    curve = points[walk]; center = curve.mean(0)
    vector_area = .5 * np.cross(curve-center, np.roll(curve, -1, axis=0)-center).sum(0)
    area = float(np.linalg.norm(vector_area))
    if not np.isfinite(area) or area <= 0:
        raise ValueError('Collapsed or non-finite lumen reference')
    return dict(representation='existing-shared-lumen-rim', first_part=records[0]['id'], second_part=records[1]['id'],
                points=curve.tolist(), centroid_m=center.tolist(), edge_count=len(walk),
                shared_tissue_triangles=len(seams[0]), opposed_normal_words_verified=all(normal_presence),
                projected_vector_area_m2=area,
                plane_residual_m=float(np.max(abs((curve-center) @ (vector_area/area)))),
                new_geometry_bytes=0, anatomical_acceptance=False)
