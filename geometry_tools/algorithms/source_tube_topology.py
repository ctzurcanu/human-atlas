"""Recover source annuli without imposing circular or planar tube profiles.

The caller supplies the exporter's actual triangle pairing. This is topology
evidence, not inferred anatomy, a centerline, a finite wall or a lumen proof.
"""
from collections import Counter, defaultdict
import numpy as np


def _edge(a, b):
    return tuple(sorted((int(a), int(b))))


def _cycles(edges):
    graph = defaultdict(set)
    for a, b in edges:
        graph[int(a)].add(int(b)); graph[int(b)].add(int(a))
    left = set(graph); cycles = []
    while left:
        start = min(left); component = set(); stack = [start]
        while stack:
            v = stack.pop()
            if v in component:
                continue
            component.add(v); stack.extend(graph[v] - component)
        left -= component
        if not all(len(graph[v]) == 2 for v in component):
            continue
        order = [start]; previous = None
        while len(order) < len(component):
            choices = graph[order[-1]] - set(order)
            if not choices:
                raise ValueError('Source cycle does not have a simple boundary')
            previous = min(choices); order.append(previous)
        if order[0] not in graph[order[-1]]:
            raise ValueError('Unclosed source cycle')
        cycles.append(order)
    return cycles


def _annulus(faces, first_ring, second_ring):
    vertices = set(map(int, faces.ravel()))
    if len(vertices) != len(first_ring) + len(second_ring):
        return False
    edge_faces = defaultdict(list); directed = Counter()
    for i, tri in enumerate(faces):
        for a, b in zip(tri, np.roll(tri, -1)):
            edge_faces[_edge(a, b)].append(i); directed[int(a), int(b)] += 1
    if any(len(v) not in (1, 2) for v in edge_faces.values()):
        return False
    if len(vertices) - len(edge_faces) + len(faces) != 0:
        return False
    if any(len(ids) == 2 and (directed[e] != 1 or directed[e[::-1]] != 1)
           for e, ids in edge_faces.items()):
        return False
    adjacency = defaultdict(set)
    for ids in edge_faces.values():
        if len(ids) == 2:
            a, b = ids; adjacency[a].add(b); adjacency[b].add(a)
    seen = set(); stack = [0]
    while stack:
        i = stack.pop()
        if i not in seen:
            seen.add(i); stack.extend(adjacency[i] - seen)
    if len(seen) != len(faces):
        return False
    rim = [e for e, ids in edge_faces.items() if len(ids) == 1]
    cycles = _cycles(rim)
    return len(cycles) == 2 and {frozenset(c) for c in cycles} == {
        frozenset(first_ring), frozenset(second_ring)}


def recover_source_annuli(vertices, faces, triangle_pairing, *, ring_size=6):
    """Return source point aliases, original face IDs and certified annular strips.

    No coordinates/attributes move. Exact-equal positions are used only for
    incidence. Non-circular/warped profiles are retained; fork collars and
    unpaired faces stay explicit instead of being invented as regular rings.
    """
    v = np.asarray(vertices, float); raw = np.asarray(faces)
    pairing = np.asarray(triangle_pairing)
    if v.ndim != 2 or v.shape[1] != 3 or not np.isfinite(v).all():
        raise ValueError('Finite metre coordinates required')
    if raw.ndim != 2 or raw.shape[1] != 3 or not np.issubdtype(raw.dtype, np.integer):
        raise ValueError('Original integer triangle indices required')
    if pairing.ndim != 2 or pairing.shape[1] != 2 or not np.issubdtype(pairing.dtype, np.integer):
        raise ValueError('Explicit original source triangle pairs required')
    if ring_size < 3 or raw.min() < 0 or raw.max() >= len(v) or pairing.min() < 0 or pairing.max() >= len(raw):
        raise ValueError('Source index outside its bank')
    if len(np.unique(pairing)) != pairing.size:
        raise ValueError('A source face cannot belong to multiple exporter quads')
    physical, representative, inverse = np.unique(v, axis=0, return_index=True, return_inverse=True)
    f = inverse[raw]; constraints = defaultdict(list); quads = []; invalid = []
    for qid, ids in enumerate(pairing):
        ff = f[ids]; counts = Counter(_edge(a, b) for tri in ff for a, b in zip(tri, np.roll(tri, -1)))
        rim = [e for e, n in counts.items() if n == 1]; cycles = _cycles(rim)
        if len(np.unique(ff)) != 4 or len(cycles) != 1 or len(cycles[0]) != 4:
            invalid.append(qid); continue
        cycle = cycles[0]; edges = [_edge(a, b) for a, b in zip(cycle, np.roll(cycle, -1))]
        for i, e in enumerate(edges):
            constraints[e].extend([(edges[(i + 1) % 4], 1), (edges[(i + 2) % 4], 0)])
        quads.append((qid, edges))
    colors = {}; classes = []
    for start in sorted(constraints):
        if start in colors:
            continue
        colors[start] = 0; todo = [start]; group = [[], []]
        while todo:
            e = todo.pop(); group[colors[e]].append(e)
            for other, parity in constraints[e]:
                expected = colors[e] ^ parity
                if other not in colors:
                    colors[other] = expected; todo.append(other)
                elif colors[other] != expected:
                    raise ValueError('Source quad direction parity is inconsistent at a junction')
        classes.append(group)
    rings = []
    for group in classes:
        choices = [[c for c in _cycles(es) if len(c) == ring_size] for es in group]
        if choices[0] and choices[1]:
            raise ValueError('Both source quad directions contain the requested cycle size')
        rings.extend(choices[0] or choices[1])
    owner = {}
    for i, ring in enumerate(rings):
        for point in ring:
            if point in owner:
                raise ValueError('Ambiguous overlapping source cross-section cycles')
            owner[point] = i
    links = defaultdict(list)
    for qid, _ in quads:
        ids = {owner.get(int(x), -1) for x in f[pairing[qid]].ravel()}
        if -1 not in ids and len(ids) == 2:
            links[tuple(sorted(ids))].extend(map(int, pairing[qid]))
    strips = []; used_faces = set()
    for (a, b), ids in links.items():
        if len(ids) != 2 * ring_size or not _annulus(f[ids], rings[a], rings[b]):
            continue
        if used_faces.intersection(ids):
            raise ValueError('Recovered annuli duplicate original source faces')
        used_faces.update(ids); strips.append(dict(rings=[a, b], source_face_ids=sorted(ids)))
    return dict(schema='source-quad-annular-topology/v1', ring_size=ring_size,
                physical_points=physical, representative_raw_ids=representative,
                raw_to_physical=inverse, rings=rings, strips=strips,
                invalid_quad_pair_ids=invalid,
                unassigned_source_face_ids=sorted(set(range(len(raw))) - used_faces),
                exact_position_incidence_only=True, geometry_changed=False,
                planar_or_circular_profile_assumed=False, anatomical_acceptance=False)


def recover_regular_source_annuli(vertices, faces, triangle_pairing, *, ring_size=6):
    """Isolate source-grid singularities before recovering regular strips.

    Forks can invalidate a global quad-direction coloring. Only literal
    exporter quads whose incident vertices have regular quad valence enter
    the strip recovery: four quads internally, two at a surface boundary.
    Every excluded/unrecovered face remains owned, unchanged core data.
    These components are candidates for further junction review, not inferred
    anatomical forks. No tolerance weld, surface motion or face removal occurs.
    """
    v=np.asarray(vertices,float);raw=np.asarray(faces);pairing=np.asarray(triangle_pairing)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all():
        raise ValueError('Finite metre coordinates required')
    if raw.ndim!=2 or raw.shape[1]!=3 or not np.issubdtype(raw.dtype,np.integer) or not len(raw):
        raise ValueError('Nonempty original integer triangles required')
    if pairing.ndim!=2 or pairing.shape[1]!=2 or not np.issubdtype(pairing.dtype,np.integer) or not len(pairing):
        raise ValueError('Nonempty explicit original source triangle pairs required')
    if raw.min()<0 or raw.max()>=len(v) or pairing.min()<0 or pairing.max()>=len(raw):
        raise ValueError('Source index outside its bank')
    if np.unique(pairing).size!=pairing.size:
        raise ValueError('A source face cannot belong to multiple exporter quads')
    physical,representative,inverse=np.unique(v,axis=0,return_index=True,return_inverse=True)
    f=inverse[raw];edge_owners=defaultdict(list)
    for i,tri in enumerate(f):
        for a,b in zip(tri,np.roll(tri,-1)):edge_owners[_edge(a,b)].append(i)
    if any(len(ids)>2 for ids in edge_owners.values()):
        raise ValueError('Nonmanifold source cannot define regular junction collars')
    boundary=set(x for e,ids in edge_owners.items() if len(ids)==1 for x in e)
    quad_vertices=[];invalid=[];valence=np.zeros(len(physical),int)
    for i,ids in enumerate(pairing):
        ff=f[ids];counts=Counter(_edge(a,b) for tri in ff for a,b in zip(tri,np.roll(tri,-1)))
        cycles=_cycles([e for e,n in counts.items() if n==1]);points=np.unique(ff)
        if len(points)!=4 or len(cycles)!=1 or len(cycles[0])!=4:
            invalid.append(i);quad_vertices.append(points);continue
        quad_vertices.append(points);valence[points]+=1
    irregular=set(i for i,n in enumerate(valence) if n!=(2 if i in boundary else 4))
    for i in invalid:irregular.update(map(int,quad_vertices[i]))
    unpaired=sorted(set(range(len(raw)))-set(map(int,pairing.ravel())))
    if unpaired:irregular.update(map(int,f[unpaired].ravel()))
    invalid_set=set(invalid)
    regular=[i for i,points in enumerate(quad_vertices) if i not in invalid_set and not irregular.intersection(map(int,points))]
    if regular:
        result=recover_source_annuli(v,raw,pairing[regular],ring_size=ring_size)
    else:
        result=dict(ring_size=ring_size,physical_points=physical,representative_raw_ids=representative,
                    raw_to_physical=inverse,rings=[],strips=[],unassigned_source_face_ids=list(range(len(raw))))
    owned=[i for strip in result['strips'] for i in strip['source_face_ids']]
    if len(owned)!=len(set(owned)):
        raise ValueError('Regular strips duplicate actual source faces')
    core=set(range(len(raw)))-set(owned);adj=defaultdict(set)
    for ids in edge_owners.values():
        if len(ids)==2 and ids[0] in core and ids[1] in core:
            adj[ids[0]].add(ids[1]);adj[ids[1]].add(ids[0])
    unseen=set(core);components=[]
    while unseen:
        todo=[min(unseen)];component=set()
        while todo:
            i=todo.pop()
            if i not in unseen:continue
            unseen.remove(i);component.add(i);todo.extend(adj[i])
        components.append(sorted(component))
    assert sorted(owned+list(core))==list(range(len(raw)))
    result.update(schema='junction-separated-source-annuli/v1',
                  original_pairing=pairing.copy(),regular_quad_pair_ids=regular,
                  excluded_quad_pair_ids=sorted(set(range(len(pairing)))-set(regular)),
                  invalid_quad_pair_ids=invalid,irregular_physical_vertex_ids=sorted(irregular),
                  physical_quad_valence=valence,unassigned_source_face_ids=sorted(core),
                  unrecovered_source_face_components=components,
                  source_face_ownership_complete=True,core_anatomical_identity_established=False,
                  geometry_changed=False,exact_position_incidence_only=True,
                  planar_or_circular_profile_assumed=False,anatomical_acceptance=False)
    return result
