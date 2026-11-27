"""Actual source-facet metric/normal inequalities with analytic derivatives."""
import numpy as np


class PatchMetricConstraints:
    def __init__(self, vertices, faces, moving, *, scale_m, strain, rotation_degrees):
        self.vertices=np.asarray(vertices,float);self.faces=np.asarray(faces,int);self.moving=np.asarray(moving,int)
        self.scale=float(scale_m);self.strain=float(strain);self.cos=np.cos(np.radians(rotation_degrees))
        if not 0<self.scale or not 0<self.strain<1 or len(np.unique(self.moving))!=len(self.moving):
            raise ValueError('Distinct controls, positive scale and finite strain domain required')
        x=self.vertices[self.faces];E=np.swapaxes(x[:,1:]-x[:,:1],1,2)
        u,s,_=np.linalg.svd(E,full_matrices=False)
        if np.any(s[:,-1]<=1e-14):raise ValueError('Nondegenerate source facets required')
        T=np.linalg.pinv(E)@u
        self.C=np.concatenate([-T.sum(1)[:,None,:],T],axis=1)
        n=np.cross(E[:,:,0],E[:,:,1]);self.normal=n/np.linalg.norm(n,axis=1)[:,None]
        self.mapping=np.full(len(self.vertices),-1,int);self.mapping[self.moving]=np.arange(len(self.moving))

    def construct(self,z):
        q=self.vertices.copy();q[self.moving]+=np.asarray(z).reshape(-1,3)*self.scale;return q

    def evaluate(self,z):
        q=self.construct(z);tri=q[self.faces];J=np.einsum('fvi,fvj->fij',tri-tri[:,:1],self.C)
        gram=np.swapaxes(J,1,2)@J;values,vectors=np.linalg.eigh(gram)
        e1=tri[:,1]-tri[:,0];e2=tri[:,2]-tri[:,0];n=np.cross(e1,e2);length=np.linalg.norm(n,axis=1)
        if np.any(length<=1e-24):raise ValueError('Collapsed proposed facet')
        nh=n/length[:,None];align=np.einsum('ij,ij->i',nh,self.normal)
        count=len(self.faces);width=3*len(self.moving);g=np.zeros((count,2,3,3))
        for k in range(2):
            v=vectors[:,:,k];jv=np.einsum('fij,fj->fi',J,v);cv=np.einsum('fvj,fj->fv',self.C,v)
            g[:,k]=2*cv[:,:,None]*jv[:,None,:]*self.scale
        normal_adjoint=(self.normal-align[:,None]*nh)/length[:,None]
        grad1=np.cross(e2,normal_adjoint);grad2=np.cross(normal_adjoint,e1)
        gn=np.stack([-grad1-grad2,grad1,grad2],axis=1)*self.scale
        jac_values=np.zeros((count,2,width));jac_normal=np.zeros((count,width))
        for corner in range(3):
            ids=self.mapping[self.faces[:,corner]]
            valid=np.flatnonzero(ids>=0)
            for axis in range(3):
                cols=ids[valid]*3+axis
                for k in range(2):jac_values[valid,k,cols]+=g[valid,k,corner,axis]
                jac_normal[valid,cols]+=gn[valid,corner,axis]
        constraints=np.r_[(values-(1-self.strain)**2).ravel(),((1+self.strain)**2-values).ravel(),align-self.cos]
        jac=np.concatenate([jac_values.reshape(-1,width),-jac_values.reshape(-1,width),jac_normal],axis=0)
        return constraints,jac
