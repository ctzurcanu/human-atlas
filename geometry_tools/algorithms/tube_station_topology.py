"""Recover serial section rings from an unbranched triangulated tube strip.

Only source connectivity supplies the sections. No nearest-neighbour joining,
anatomical radius, branch identity, inner surface or clinical claim is inferred.
Disconnected source cylinders are handled separately; real forks are rejected.
"""
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

class StationTopologyError(ValueError):pass

def recover_station_rings(positions,faces):
    v=np.asarray(positions,float);f=np.asarray(faces,int)
    if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3:raise StationTopologyError('Triangle tube arrays required')
    edges=np.vstack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);graph=coo_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),len(v))).tocsr();_,labels=connected_components(graph,directed=False);used=np.unique(f);adj=defaultdict(set)
    for x,y in edges:adj[int(x)].add(int(y));adj[int(y)].add(int(x))
    result=[]
    for label in sorted(set(labels[used])):
        face_ids=np.flatnonzero(labels[f[:,0]]==label);ff=f[face_ids];mouths=boundary_loops(ff)
        if len(mouths)!=2 or len(mouths[0])!=len(mouths[1]):raise StationTopologyError('Each source cylinder needs two equal-count open mouths')
        first=mouths[0];count=len(first);groups=[first];seen=set(map(int,first));all_ids=set(map(int,ff.ravel()))
        while seen!=all_ids:
            frontier=set().union(*(adj[int(x)]for x in groups[-1]))-seen
            if len(frontier)!=count:raise StationTopologyError('Source connectivity does not define equal-size serial section rings')
            # Same-station source edges must form one simple cycle.
            cycle=[min(frontier)];previous=None
            while True:
                nxt=(adj[cycle[-1]]&frontier)-({previous}if previous is not None else set())
                if len(nxt)!=(2 if previous is None else 1):raise StationTopologyError('Section is not a degree-two source cycle')
                nxt=min(nxt)
                if nxt==cycle[0]:break
                if nxt in cycle:raise StationTopologyError('Section cycle self-repeats')
                previous=cycle[-1];cycle.append(nxt)
            if len(cycle)!=count:raise StationTopologyError('Section has disconnected source cycles')
            groups.append(np.asarray(cycle,int));seen|=frontier
        if set(map(int,groups[-1]))!=set(map(int,mouths[1])):raise StationTopologyError('Final section does not match the other source mouth')
        station=np.full(len(v),-1,int)
        for j,ring in enumerate(groups):station[ring]=j
        sf=station[ff]
        if not np.all(sf.max(1)-sf.min(1)==1):raise StationTopologyError('A source face skips a section or caps a mouth')
        result.append(dict(face_ids=face_ids,rings=tuple(groups),station_of_vertex=station,source_component=int(label),source_mouths=mouths))
    return result


def align_station_correspondence(positions, faces, rings):
    """Recover a consistent phase using the literal source strip triangles.

    recover_station_rings intentionally returns independent cycles; their
    index minima do not define corresponding material points. For each new
    interval, enumerate cyclic phases/directions, require the original face
    set to replay exactly through its quads (either diagonal per quad), then
    choose the shortest longitudinal source-edge assignment. No vertex moves
    or proximity welds occur.
    """
    v = np.asarray(positions, float); f = np.asarray(faces, int)
    rings = [np.asarray(q, int) for q in rings]
    station = np.full(len(v), -1, int)
    for j, ring in enumerate(rings): station[ring] = j
    if np.any(station[f] < 0) or not np.all(np.ptp(station[f], axis=1) == 1):
        raise StationTopologyError('Every source face must join consecutive recovered rings')
    result = [rings[0]]; rows = []
    canon = lambda t: tuple(sorted(map(int, t)))
    for interval, right in enumerate(rings[1:]):
        left = result[-1]; source_ids = np.flatnonzero(station[f].min(1) == interval)
        source_keys = {canon(t) for t in f[source_ids]}; choices = []
        for reverse in (False, True):
            for shift in range(len(right)):
                ordered = np.roll(right[::-1] if reverse else right, shift)
                replay = set(); diagonals = []
                for k in range(len(left)):
                    j = (k + 1) % len(left)
                    a, b, c, d = map(int, (left[k], left[j], ordered[j], ordered[k]))
                    first = {canon((a, c, b)), canon((a, d, c))}
                    second = {canon((a, d, b)), canon((b, d, c))}
                    if first <= source_keys: replay |= first; diagonals.append(0)
                    elif second <= source_keys: replay |= second; diagonals.append(1)
                    else: break
                if len(diagonals) == len(left) and replay == source_keys:
                    cost = float(np.sum((v[left] - v[ordered]) ** 2))
                    choices.append((cost, reverse, shift, ordered, diagonals))
        if not choices:
            raise StationTopologyError('No phase correspondence replays the literal source triangles')
        choices.sort(key=lambda q: (q[0], q[1], q[2])); cost, reverse, shift, ordered, diagonals = choices[0]
        ambiguity = len(choices) > 1 and abs(choices[1][0] - cost) <= 256 * np.finfo(float).eps * max(cost, np.finfo(float).tiny)
        if ambiguity:
            raise StationTopologyError('Source longitudinal phase is geometrically ambiguous')
        result.append(ordered)
        rows.append(dict(sourceInterval=interval, reversed=reverse, cyclicShift=shift,
            longitudinalEdgeCostMetersSquared=cost, diagonalPerQuad=diagonals,
            sourceFaceIds=source_ids.tolist(), sourceFacesReplayedExactly=True))
    return tuple(result), rows
