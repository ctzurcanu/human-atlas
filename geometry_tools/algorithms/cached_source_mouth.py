"""Per-operation boundary index for the unchanged complete-mouth predicate.

Mesh arrays are immutable. Strong references keep each identity stable for
this checker lifetime; nothing is shared between catalogue revisions/jobs.
"""
import numpy as np

from .tube_continuity import boundary_loops


class SourceMouthChecks:
    def __init__(self):
        self._meshes={}
        self.geometry_builds=0

    def _boundary(self,mesh,ids):
        ids=np.asarray(ids,int)
        if ids.ndim!=1 or len(ids)<3 or np.any(ids<0)or np.any(ids>=len(mesh.vertices)):
            raise ValueError('Every declared mouth control must exist')
        key=id(mesh)
        if key not in self._meshes:
            unique,inv=np.unique(mesh.vertices,axis=0,return_inverse=True)
            loops={}
            for loop in boundary_loops(inv[mesh.faces]):
                loops.setdefault(frozenset(loop),[]).append(loop)
            self._meshes[key]=(mesh,unique,inv,loops)
            self.geometry_builds+=1
        retained,unique,inv,loops=self._meshes[key]
        assert retained is mesh
        target=inv[ids];actual=loops.get(frozenset(target),[])
        if len(np.unique(target))!=len(target)or len(actual)!=1:
            raise ValueError('Declared controls are not one complete external source mouth')
        ring=unique[actual[0]];center=ring.mean(0)
        normal=np.cross(ring-center,np.roll(ring,-1,axis=0)-center).sum(0)
        if np.linalg.norm(normal)<=1e-20:raise ValueError('Zero-area source mouth')
        return mesh.vertices[ids],ring,normal

    def check(self,first,second,first_ids,second_ids,*,tolerance_m=1e-6,area_tolerance=.01):
        if first.frame!=second.frame:raise ValueError('Source mouths must use the same explicit coordinate frame')
        p,actual_p,na=self._boundary(first,first_ids);q,actual_q,nb=self._boundary(second,second_ids)
        curve=None
        if p.shape!=q.shape:
            from .closed_curve_correspondence import closed_curve_correspondence
            curve=closed_curve_correspondence(actual_p,actual_q,opposite=True)
            gap=curve['maximumGap_m']
        else:gap=float(np.linalg.norm(p-q,axis=1).max())
        area_error=float(abs(np.linalg.norm(nb)/np.linalg.norm(na)-1))
        cosine=float(na@nb/(np.linalg.norm(na)*np.linalg.norm(nb)))
        result=dict(passed=bool(gap<=tolerance_m and area_error<=area_tolerance and cosine<-.99),
            maximumGap_m=gap,relativeAreaError=area_error,normalCosine=cosine,
            completeSourceBoundaries=True,geometricCorrespondenceOnly=True,
            wholeLumenPatencyEstablished=False,anatomicalAttachmentEstablished=False)
        if curve is not None:result['curveCorrespondence']=curve
        return result
