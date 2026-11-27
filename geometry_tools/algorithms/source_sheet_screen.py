"""Screen whole disconnected carriers before proposing an internal mouth weld.

Close opposed boundary cycles may be the two sides of an existing source
envelope. Complete component extent and normal evidence is required; this
screen does not establish histological layers or authorize a weld.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def screen_source_components(mesh, *, maximum_separation_m=.0002,
                             minimum_opposed_fraction=.95):
    if not 0<maximum_separation_m<.001 or not 0<minimum_opposed_fraction<=1:
        raise ValueError('Explicit submillimetre separation and opposed fraction required')
    v,rep,inv=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True)
    f=inv[mesh.faces]
    if not len(f): raise ValueError('A complete nonempty source carrier required')
    edges=np.unique(np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
    graph=coo_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),len(v))).tocsr()
    used=np.unique(f); count,local=connected_components(graph[used][:,used],directed=False)
    labels=np.full(len(v),-1,int); labels[used]=local
    result=dict(complete=True,sourceGeometrySha256=mesh.geometry_sha,components=count,
                facesPerComponent=[int(np.sum(labels[f[:,0]]==i)) for i in range(count)],
                completeComponentComparison=[],pairedOpposedSourceSheets=False,
                automaticWeldAllowed=False,clinicalWallThicknessEstablished=False,
                maximumSeparation_m=maximum_separation_m,minimumOpposedFraction=minimum_opposed_fraction)
    if count!=2 or not {'normals','normalSigns'}<=mesh.attributes.keys():
        result['classification']='unresolved-component-role';return result
    normals=mesh.attributes['normals'].astype(float)*mesh.attributes['normalSigns'][:,None]
    length=np.linalg.norm(normals,axis=1)
    if np.any(length==0):
        result['classification']='unresolved-zero-source-normals';return result
    normals/=length[:,None]
    for i,j in [(0,1),(1,0)]:
        first=np.flatnonzero(labels==i);second=np.flatnonzero(labels==j)
        distance,near=cKDTree(v[second]).query(v[first])
        cosine=np.einsum('ij,ij->i',normals[rep[first]],normals[rep[second[near]]])
        result['completeComponentComparison'].append(dict(maximumDistance_m=float(distance.max()),
            medianDistance_m=float(np.median(distance)),medianNormalCosine=float(np.median(cosine)),
            opposedFraction=float(np.mean(cosine<-.9))))
    result['pairedOpposedSourceSheets']=all(r['maximumDistance_m']<=maximum_separation_m and
        r['opposedFraction']>=minimum_opposed_fraction for r in result['completeComponentComparison'])
    result['classification']='paired-opposed-source-sheets-do-not-weld' if result['pairedOpposedSourceSheets'] else 'unresolved-component-role'
    return result
