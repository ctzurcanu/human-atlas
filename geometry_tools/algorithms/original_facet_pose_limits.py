"""Original source shape plus bounded actual facet-normal turns."""
import numpy as np
from scipy import sparse
from scipy.spatial.transform import Rotation
from .original_singular_pose_rows import OriginalSingularPoseRows,skew,right_rotation_jacobian
from .rigid_section_twists import replay_section_twists


class OriginalFacetPoseLimits(OriginalSingularPoseRows):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        t=self.current.triangles
        normals=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        self.original_normals=normals/np.linalg.norm(normals,axis=1)[:,None]

    def orientation_values(self,parameters):
        q,_=replay_section_twists(self.current.vertices,self.groups,self.centers,self.levers,parameters)
        t=q[self.current.faces];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        normal/=np.linalg.norm(normal,axis=1)[:,None]
        return np.clip(np.einsum('ij,ij->i',normal,self.original_normals),-1,1)

    def orientation_violations(self,parameters,*,maximum_degrees=10.):
        angle=np.degrees(np.arccos(self.orientation_values(parameters)))
        bad=np.flatnonzero(angle>maximum_degrees+1e-9)
        return bad,dict(complete=True,allOriginalFacetsChecked=len(angle),
            violatingOriginalFacets=len(bad),maximumFacetTurn_deg=float(angle.max(initial=0)),
            maximumAllowed_deg=maximum_degrees,exactProperPoseReplay=True)

    def orientation_gradients(self,parameters,face_ids):
        x=np.asarray(parameters,float);ids=np.asarray(face_ids,int)
        q,_=replay_section_twists(self.current.vertices,self.groups,self.centers,self.levers,x)
        faces=self.current.faces[ids];t=q[faces];a=t[:,1]-t[:,0];b=t[:,2]-t[:,0]
        normal=np.cross(a,b);area=np.linalg.norm(normal,axis=1);normal/=area[:,None]
        cosine=np.einsum('ij,ij->i',normal,self.original_normals[ids])
        w=(self.original_normals[ids]-cosine[:,None]*normal)/area[:,None]
        g1=np.cross(b,w);g2=np.cross(w,a)
        point_gradient=np.stack([-g1-g2,g1,g2],axis=1)
        phi=x[self.count:]/self.levers[:,None]
        rotations=Rotation.from_rotvec(phi).as_matrix();jacobian=right_rotation_jacobian(phi)
        rr=[];cc=[];vv=[]
        for corner in range(3):
            groups=self.groups[faces[:,corner]];live=groups>=0;g=groups[live]
            offset=self.current.vertices[faces[live,corner]]-self.centers[g]
            J=-rotations[g]@skew(offset)@jacobian[g]/self.levers[g,None,None]
            translation=point_gradient[live,corner];rotation=np.einsum('ni,nij->nj',translation,J)
            rows=np.repeat(np.flatnonzero(live),3)
            rr.extend([rows,rows]);cc.extend([(g[:,None]*3+np.arange(3)).ravel(),(3*self.count+g[:,None]*3+np.arange(3)).ravel()])
            vv.extend([translation.ravel(),rotation.ravel()])
        matrix=sparse.coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(len(ids),6*self.count)).tocsr()
        matrix.eliminate_zeros()
        return cosine,matrix

    def local_orientation_rows(self,parameters,face_ids,*,maximum_degrees=10.):
        cosine,G=self.orientation_gradients(parameters,face_ids)
        # A stricter numerical target leaves the physical ten-degree bound
        # intact, rather than accepting a solver-residual overshoot.
        target=np.cos(np.radians(maximum_degrees-1e-4))
        norm=np.sqrt(np.asarray(G.power(2).sum(axis=1)).ravel());live=norm>1e-12
        if np.any(~live&(cosine<np.cos(np.radians(maximum_degrees)))):
            raise ValueError('Unresolved orientation derivative at violating source facet')
        multiplier=np.maximum(1.,norm[live]/10.)/norm[live]
        rhs=(target-cosine[live]+G[live]@np.asarray(parameters).ravel())*multiplier
        return sparse.diags(multiplier)@G[live],rhs,dict(complete=True,sourceFaceIds=np.asarray(face_ids).tolist(),
            originalFacetTurnBound_deg=maximum_degrees,numericalTarget_deg=maximum_degrees-1e-4,
            exactReplayRequired=True,localLinearizationOnly=True)


def original_facet_limits(reference,current,groups,centers,levers,*,baseline=None):
    """Rebase pose coordinates without resetting the original shape/turn budgets."""
    model=OriginalFacetPoseLimits(reference,current,groups,centers,levers)
    if baseline is not None:
        original=OriginalFacetPoseLimits(reference,baseline,groups,centers,levers)
        model.budget=original.budget.copy()
        model.original_normals=original.original_normals.copy()
    return model
