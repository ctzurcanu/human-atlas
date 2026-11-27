"""Regular Freudenthal tetrahedral cages with fixed, literal outer boundaries.

Positive element Jacobians and a fixed convex boundary provide the continuous
piecewise affine embedding. Source triangles still must conform to tetrahedral
faces before a mapped surface can inherit that embedding's geometric proof.
"""
from itertools import permutations
import numpy as np


class RegularTetrahedralCage:
    def __init__(self,minimum,step,shape):
        self.minimum=np.asarray(minimum,float);self.step=float(step);self.shape=np.asarray(shape,int)
        if self.minimum.shape!=(3,)or self.shape.shape!=(3,)or np.any(self.shape<2)or self.step<=0 \
                or not np.isfinite(np.r_[self.minimum,self.step]).all():
            raise ValueError('Finite minimum, positive spacing and three node counts >= 2 required')
        self.maximum=self.minimum+self.step*(self.shape-1)
        indices=np.indices(tuple(self.shape)).reshape(3,-1).T
        self.nodes=self.minimum+indices*self.step
        self.boundary=np.any((indices==0)|(indices==self.shape-1),axis=1)
        self.permutations=np.array(list(permutations(range(3))),int)
        self._permutation_codes=self.permutations@np.array([9,3,1])
        cube=np.indices(tuple(self.shape-1)).reshape(3,-1).T
        tetra=[]
        for c in cube:
            for order in self.permutations:
                corners=np.array([c,c+np.eye(3,dtype=int)[order[0]],
                    c+np.eye(3,dtype=int)[order[:2]].sum(0),c+1])
                tetra.append(self.node_ids(corners))
        self.tetrahedra=np.asarray(tetra,int)
        source=self.nodes[self.tetrahedra]
        self.source_edges=np.transpose(source[:,1:]-source[:,:1],(0,2,1))
        self.inverse_source=np.linalg.inv(self.source_edges)

    def node_ids(self,index):
        q=np.asarray(index,int)
        return (q[...,0]*self.shape[1]+q[...,1])*self.shape[2]+q[...,2]

    def conformance_planes(self,stage):
        """Planes of the regular cage, grouped into finite saved stages."""
        result=[]
        if stage=='box':
            for axis in range(3):
                for coordinate in [0.,(self.shape[axis]-1)*self.step]:
                    result.append((np.eye(3)[axis],float(coordinate),'overlap'))
        elif stage=='grid':
            for axis in range(3):
                for k in range(1,self.shape[axis]-1):
                    result.append((np.eye(3)[axis],float(k*self.step),'inside'))
        elif stage in ['diagxy','diagxz','diagyz']:
            a,b={'diagxy':(0,1),'diagxz':(0,2),'diagyz':(1,2)}[stage]
            for k in range(-(self.shape[b]-2),self.shape[a]-1):
                result.append((np.eye(3)[a]-np.eye(3)[b],float(k*self.step),'inside'))
        else:raise ValueError('Unknown saved cage-conformance stage')
        return result

    def active_faces(self,positions,faces,scope,*,tolerance_m=2e-12):
        q=np.asarray(positions,float)[np.asarray(faces,int)];lo=q.min(1);hi=q.max(1)
        if scope=='overlap':return np.all((hi>=self.minimum-tolerance_m)&(lo<=self.maximum+tolerance_m),axis=1)
        if scope=='inside':return np.all((lo>=self.minimum-tolerance_m)&(hi<=self.maximum+tolerance_m),axis=1)
        raise ValueError('Unknown cage facet scope')

    def triangle_conformance(self,positions,faces,*,tolerance_m=2e-11):
        q=np.asarray(positions,float)[np.asarray(faces,int)];location=self.locate(q.mean(1))
        ids=np.flatnonzero(location['inside']);t=location['tetrahedron'][ids]
        source=self.nodes[self.tetrahedra[t]]
        coordinates=np.einsum('tij,tkj->tki',self.inverse_source[t],q[ids]-source[:,:1])
        bary=np.concatenate([1-coordinates.sum(2,keepdims=True),coordinates],axis=2)
        violation=np.maximum(-bary.min((1,2)),0)*self.step
        bad=ids[violation>tolerance_m]
        return dict(passed=len(bad)==0,triangles_inside_cage=len(ids),nonconforming_triangles=len(bad),
            examples=bad[:20].tolist(),maximum_violation_m=float(violation.max(initial=0)),
            declared_tolerance_m=tolerance_m,storage_and_anatomical_checks_pending=True)

    def locate(self,points):
        p=np.asarray(points,float)
        if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite 3D points required')
        scaled=(p-self.minimum)/self.step
        inside=np.all((scaled>=0)&(scaled<=self.shape-1),axis=1)
        cube=np.clip(np.floor(scaled).astype(int),0,self.shape-2)
        fraction=np.clip(scaled-cube,0,1)
        order=np.argsort(-fraction,axis=1,kind='stable')
        sorted_fraction=np.take_along_axis(fraction,order,axis=1)
        weights=np.c_[1-sorted_fraction[:,0],sorted_fraction[:,0]-sorted_fraction[:,1],
                      sorted_fraction[:,1]-sorted_fraction[:,2],sorted_fraction[:,2]]
        first=np.eye(3,dtype=int)[order[:,0]]
        second=first+np.eye(3,dtype=int)[order[:,1]]
        indices=np.stack([cube,cube+first,cube+second,cube+1],axis=1)
        nodes=self.node_ids(indices);weights[~inside]=0
        codes=order@np.array([9,3,1]);permutation=np.argmax(codes[:,None]==self._permutation_codes[None],axis=1)
        cell=(cube[:,0]*(self.shape[1]-1)+cube[:,1])*(self.shape[2]-1)+cube[:,2]
        tetrahedron=cell*6+permutation;tetrahedron[~inside]=-1
        return dict(node_ids=nodes,weights=weights,inside=inside,tetrahedron=tetrahedron)

    def jacobians(self,displacement):
        d=np.asarray(displacement,float)
        if d.shape!=self.nodes.shape or not np.isfinite(d).all():raise ValueError('One finite displacement per cage node required')
        q=(self.nodes+d)[self.tetrahedra]
        edges=np.transpose(q[:,1:]-q[:,:1],(0,2,1))
        return np.einsum('tij,tjk->tik',edges,self.inverse_source)

    def check(self,displacement,*,minimum_jacobian=0.):
        d=np.asarray(displacement,float);j=self.jacobians(d);det=np.linalg.det(j)
        boundary_literal=bool(np.all(d[self.boundary]==0))
        return dict(passed=boundary_literal and bool(np.all(det>minimum_jacobian)),
            boundary_literal=boundary_literal,minimum_jacobian=float(det.min()),
            declared_minimum_jacobian=float(minimum_jacobian),tetrahedra=len(j),
            source_surface_conformance_required=True,anatomical_acceptance=False)

    def map(self,points,displacement,*,minimum_jacobian=0.):
        proof=self.check(displacement,minimum_jacobian=minimum_jacobian)
        if not proof['passed']:raise ValueError('Cage boundary or positive-element gate failed')
        p=np.asarray(points,float);location=self.locate(p)
        return p+np.einsum('nk,nki->ni',location['weights'],np.asarray(displacement)[location['node_ids']])

    def point_jacobians(self,points,displacement):
        location=self.locate(points);out=np.broadcast_to(np.eye(3),(len(points),3,3)).copy()
        take=location['inside'];out[take]=self.jacobians(displacement)[location['tetrahedron'][take]]
        return out
