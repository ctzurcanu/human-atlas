"""Recover literal collars of an explicitly declared closed source root.

Incidence is exact-coordinate only. This certifies source-cell ownership,
never clinical vessel identity or an anastomosis inferred from proximity.
"""
from collections import defaultdict,Counter
import numpy as np
from .source_tube_topology import _annulus
from .tube_continuity import boundary_loops


def recover_threeway_cells(vertices,faces,root_ids):
    v=np.asarray(vertices,float);f=np.asarray(faces,int)
    physical,rep,inv=np.unique(v,axis=0,return_index=True,return_inverse=True)
    root=np.asarray(root_ids,int);r=inv[root];rootset=set(map(int,r))
    if len(root)!=6 or len(rootset)!=6:raise ValueError('One ordered six-control source root required')
    ff=inv[f];star=np.flatnonzero(np.isin(ff,r).any(1))
    if len(star)!=36:raise ValueError('Closed three-way root does not own exactly36literal collar faces')
    incidence=Counter(tuple(sorted(map(int,[a,b])))for t in ff for a,b in zip(t,np.roll(t,-1)))
    if any(incidence[tuple(sorted(map(int,[a,b])))]!=3 for a,b in zip(r,np.roll(r,-1))):
        raise ValueError('Every declared root edge must have exactly three source faces')
    graph=defaultdict(set)
    for fi in star:
        outside=set(map(int,ff[fi]))-rootset
        if len(outside)not in [1,2]:raise ValueError('A source root star contains a cap or non-collar triangle')
        for p in outside:graph[p].update(outside-{p})
    unseen=set(graph);components=[]
    while unseen:
        todo=[min(unseen)];seen=set()
        while todo:
            p=todo.pop()
            if p not in seen:seen.add(p);todo.extend(graph[p]-seen)
        unseen-=seen;components.append(seen)
    if len(components)!=3 or any(len(c)!=6 for c in components):
        raise ValueError('Source star does not separate into three distinct six-control outer profiles')
    cells=[]
    for port in components:
        ids=[int(fi)for fi in star if (set(map(int,ff[fi]))-rootset)<=port]
        loops=boundary_loops(ff[ids]);outer=next((p for p in loops if set(p)==port),None)
        if len(ids)!=12 or outer is None or not _annulus(ff[ids],r,outer):
            raise ValueError('Recovered source cell is not an oriented twelve-triangle annulus')
        cells.append(dict(sourceFaceIds=ids,rootRawIds=root.tolist(),outerRawIds=rep[outer].tolist()))
    if len({i for cell in cells for i in cell['sourceFaceIds']})!=36:
        raise ValueError('Recovered source cells duplicate source faces')
    return dict(complete=True,cells=cells,rootRawIds=root.tolist(),allSourceCoordinatesLiteral=True,
        sourceFacesPartitionedExactly=True,anatomicalAcceptance=False)
