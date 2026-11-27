"""Recover only certified regular source annuli, retaining a nonmanifold core.

Unlike the strict whole-carrier utility this accepts known singularities only as
explicit unreconstructed core. It never calls those singularities regular.
"""
from collections import Counter,defaultdict
import numpy as np
from .source_tube_topology import _edge,_cycles,recover_source_annuli

def recover_nonsingular_source_annuli(vertices, faces, triangle_pairing, *, ring_size=6):
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
    singular_edges=[e for e,ids in edge_owners.items()if len(ids)>2]
    singular_vertices={i for e in singular_edges for i in e}
    boundary=set(x for e,ids in edge_owners.items() if len(ids)==1 for x in e)
    quad_vertices=[];invalid=[];valence=np.zeros(len(physical),int)
    for i,ids in enumerate(pairing):
        ff=f[ids];counts=Counter(_edge(a,b) for tri in ff for a,b in zip(tri,np.roll(tri,-1)))
        cycles=_cycles([e for e,n in counts.items() if n==1]);points=np.unique(ff)
        if len(points)!=4 or len(cycles)!=1 or len(cycles[0])!=4:
            invalid.append(i);quad_vertices.append(points);continue
        quad_vertices.append(points);valence[points]+=1
    irregular=set(i for i,n in enumerate(valence) if n!=(2 if i in boundary else 4))
    irregular.update(singular_vertices)
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
    result.update(schema='nonmanifold-core-separated-source-annuli/v1',
                  singular_source_edges=singular_edges,singular_source_vertices=sorted(singular_vertices),
                  singular_core_not_reconstructed=True,
                  original_pairing=pairing.copy(),regular_quad_pair_ids=regular,
                  excluded_quad_pair_ids=sorted(set(range(len(pairing)))-set(regular)),
                  invalid_quad_pair_ids=invalid,irregular_physical_vertex_ids=sorted(irregular),
                  physical_quad_valence=valence,unassigned_source_face_ids=sorted(core),
                  unrecovered_source_face_components=components,
                  source_face_ownership_complete=True,core_anatomical_identity_established=False,
                  geometry_changed=False,exact_position_incidence_only=True,
                  planar_or_circular_profile_assumed=False,anatomical_acceptance=False)
    return result
