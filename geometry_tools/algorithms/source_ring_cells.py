"""Recover declared manifold annuli at simple nonmanifold source rings.

The local star is separated at nonmanifold edges. Remote vascular cycles may
reconnect two branches; global connected-component identity is not a local
branch-incidence test. No geometric movement or proximity welding is used.
"""
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .source_root_sections import recover_root_sections

def explicit_ring_cells(mesh,*,arity):
    assert isinstance(arity,int)and arity>=3
    vertices,rep,inv=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True);f=inv[mesh.faces]
    e=np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1);edges,which,count=np.unique(e,axis=0,return_inverse=True,return_counts=True);order=np.argsort(which);at=np.r_[0,np.cumsum(count)]
    pairs=np.asarray([order[at[i]:at[i+1]]//3 for i in np.flatnonzero(count==2)],int)
    graph=coo_matrix((np.ones(2*len(pairs)),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(f),len(f))).tocsr()
    adj=defaultdict(set);edge_index={tuple(map(int,edge)):i for i,edge in enumerate(edges)}
    for a,b in edges[count>2]:adj[int(a)].add(int(b));adj[int(b)].add(int(a))
    unseen=set(adj);contracts=[];deferred=[]
    while unseen:
        root=set();stack=[min(unseen)]
        while stack:
            p=stack.pop()
            if p in root:continue
            root.add(p);stack+=list(adj[p]-root)
        unseen-=root;ids=np.asarray([edge_index[tuple(sorted((a,b)))]for a in root for b in adj[a]if a<b],int)
        if len(root)<3 or len(ids)!=len(root)or not all(len(adj[p])==2 for p in root)or not np.all(count[ids]==arity):
            deferred.append(dict(vertices=len(root),edges=len(ids),reason='Not one simple ring of the declared sheet incidence'));continue
        near=np.flatnonzero(np.any(np.isin(f,list(root)),axis=1));_,labels=connected_components(graph[near][:,near],directed=False);lookup=np.full(len(f),-1,int);lookup[near]=labels
        owners=[set(lookup[order[at[i]:at[i+1]]//3])for i in ids]
        if len(owners[0])!=arity or any(owner!=owners[0]for owner in owners):
            deferred.append(dict(vertices=len(root),edges=len(ids),reason='Local root edges do not meet the same distinct source collars'));continue
        rawroot=rep[np.asarray(sorted(root))];cells=[];failed=None
        for label in sorted(owners[0]):
            members=near[labels==label]
            try:sections=recover_root_sections(mesh.vertices,mesh.faces[members],rawroot,maximum_layers=1)
            except (ValueError,AssertionError)as error:failed=str(error);break
            if len(sections)!=1:failed='One complete source annulus is not available';break
            section=sections[0];cells.append(dict(sourceFaceIds=members[np.asarray(section['sourceFaceIds'])].tolist(),distalPortRawIds=section['rawIds'],certifiedManifoldAnnulus=section['certifiedManifoldAnnulus']))
        if failed:deferred.append(dict(vertices=len(root),edges=len(ids),reason=failed));continue
        assert len(cells)==arity
        contracts.append(dict(rootRawIds=rawroot.tolist(),cells=cells,center_m=vertices[list(root)].mean(0).tolist(),sourceSheetIncidence=arity,portsGeometricallyDefined=True,clinicalAttachmentEstablished=False))
    return dict(complete=True,contracts=contracts,deferredComponents=deferred,exactCoordinateWeldingOnly=True,sourceGeometryChanged=False,clinicalAcceptance=False)
