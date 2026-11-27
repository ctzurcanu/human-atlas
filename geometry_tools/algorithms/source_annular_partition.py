"""Cut source face ownership at existing annular edges without moving tissue."""
import numpy as np


def annular_face_partition(faces, raw_to_physical, ring, *, root_face):
    f=np.asarray(faces,int);aliases=np.asarray(raw_to_physical,int);ring=np.asarray(ring,int)
    if f.ndim!=2 or f.shape[1]!=3 or aliases.ndim!=1 or not 0<=root_face<len(f)or ring.ndim!=1 or len(ring)<3 or len(set(ring))!=len(ring):
        raise ValueError('Valid source triangles, exact aliases, simple ring and root face required')
    physical=aliases[f];edges={}
    for i,tri in enumerate(physical):
        if len(set(tri))!=3:raise ValueError('Degenerate source face cannot be partitioned')
        for a,b in zip(tri,np.roll(tri,-1)):edges.setdefault(tuple(sorted((int(a),int(b)))),[]).append(i)
    if any(len(v)>2 for v in edges.values()):raise ValueError('Source nonmanifold edges need repair before partition')
    cut={tuple(sorted((int(a),int(b))))for a,b in zip(ring,np.roll(ring,-1))}
    if any(len(edges.get(e,[]))!=2 for e in cut):raise ValueError('Every seam edge must be a literal internal source edge')
    def components(excluded):
        adj=[[]for _ in f]
        for e,owners in edges.items():
            if len(owners)==2 and e not in excluded:
                a,b=owners;adj[a].append(b);adj[b].append(a)
        labels=np.full(len(f),-1,int);n=0
        for seed in range(len(f)):
            if labels[seed]>=0:continue
            stack=[seed];labels[seed]=n
            while stack:
                for j in adj[stack.pop()]:
                    if labels[j]<0:labels[j]=n;stack.append(j)
            n+=1
        return labels,n
    _,before=components(set());labels,after=components(cut)
    if before!=1 or after!=2:raise ValueError('A branch collar must split one connected surface into exactly two regions')
    root=labels[root_face];trunk=np.flatnonzero(labels==root);distal=np.flatnonzero(labels!=root)
    shared=set(physical[trunk].ravel())&set(physical[distal].ravel())
    if shared!=set(ring):raise ValueError('Shared boundary differs from the declared complete annulus')
    return trunk,distal,dict(complete=True,sourceComponents=before,partitionComponents=after,seamPhysicalVertices=ring.tolist(),seamEdges=len(cut),sourceFacesDuplicated=0,newFaces=0,coordinatesMoved=False,anatomicalBoundaryAccepted=False)
