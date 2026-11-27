"""Replay Float32 atlas encoding while fitting; keep original physical limits."""
import numpy as np
from scipy import sparse
from .packed_world_positions import encode_world_positions,decode_source_positions
from .rigid_section_twists import replay_section_twists

class PackedPoseLimits:
    def __init__(self,current,sources,parts,models,groups,centers,levers,planes=None,*,global_turn_degrees=10.,step_models=None,step_turn_degrees=10.):
        self.current=current;self.sources=sources;self.parts=parts;self.models=models
        self.groups=groups;self.centers=centers;self.levers=levers;self.planes=planes
        self.count=len(centers);self.lookup={}
        self.global_turn_degrees=global_turn_degrees;self.step_models=step_models;self.step_turn_degrees=step_turn_degrees
        for actor,(mesh,group) in enumerate(zip(current,groups)):
            for vertex,g in enumerate(group):
                if g>=0:self.lookup[(int(g),mesh.vertices[vertex].tobytes())]=(actor,vertex)
        self.plane_indices=[]
        if planes is not None:
            for points,owner in [(planes.a,planes.ga),(planes.b,planes.gb)]:
                self.plane_indices.append([self.lookup[(int(g),p.tobytes())] if g>=0 else None for p,g in zip(points,owner)])

    def vertices(self,parameters):
        exact=[];packed=[]
        for mesh,source,groups in zip(self.current,self.sources,self.groups):
            q,_=replay_section_twists(mesh.vertices,groups,self.centers,self.levers,parameters)
            result=mesh.vertices.copy();seen=np.zeros(len(result),bool)
            for projection in source['activeProjections']:
                ids=np.asarray(projection['carrierControlIds'],int);part=self.parts[projection['partId']]
                value=mesh.vertices[ids].copy();moved=np.any(q[ids]!=value,axis=1)
                value[moved]=decode_source_positions(encode_world_positions(q[ids[moved]],part),part)
                np.testing.assert_array_equal(result[ids[seen[ids]]],value[seen[ids]])
                result[ids]=value;seen[ids]=True
            exact.append(q);packed.append(result)
        return exact,packed

    def rows(self,parameters):
        exact,packed=self.vertices(parameters);blocks=[];bounds=[];details=[];violations=0
        for actor,(mesh,model,q,p) in enumerate(zip(self.current,self.models,exact,packed)):
            t=p[mesh.faces]
            singular=np.linalg.svd((t[:,1:]-t[:,:1]).transpose(0,2,1)@model.chart,compute_uv=False)
            strain=np.abs(singular-1).max(1);budget=model.budget.copy()
            bad=np.flatnonzero(strain-budget>1e-8);violations+=len(bad)
            if len(bad):
                reserve=np.maximum(strain-np.abs(model.values(parameters)-1).max(1),0)+5e-7
                try:
                    model.budget=budget.copy();model.budget[bad]-=reserve[bad]
                    A,b,_=model.local_rows(parameters,bad);blocks.append(A);bounds.append(b)
                finally:model.budget=budget
            normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
            orientations=[(model,self.global_turn_degrees)]
            if self.step_models is not None:orientations.append((self.step_models[actor],self.step_turn_degrees))
            turn_count=0
            for orientation,limit in orientations:
                angle=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',normal,orientation.original_normals),-1,1)))
                bad_turn=np.flatnonzero(angle>limit+1e-9);violations+=len(bad_turn);turn_count+=len(bad_turn)
                if len(bad_turn):
                    exact_angle=np.degrees(np.arccos(orientation.orientation_values(parameters)))
                    reserve=float(np.maximum(angle[bad_turn]-exact_angle[bad_turn],0).max())+.0002
                    A,b,_=orientation.local_orientation_rows(parameters,bad_turn,maximum_degrees=limit-reserve)
                    blocks.append(A);bounds.append(b)
            if len(bad) or turn_count:details.append(dict(actor=actor,packedShapeViolations=len(bad),packedTurnViolations=turn_count))
        plane_bad=np.array([],int)
        if self.planes is not None:
            values=[]
            for original,lookup in zip([self.planes.a,self.planes.b],self.plane_indices):
                value=original.copy()
                for j,key in enumerate(lookup):
                    if key is not None:value[j]=packed[key[0]][key[1]]
                values.append(value)
            actual=np.einsum('ij,ij->i',self.planes.axis,values[0]-values[1])
            minimum=self.planes.minimum.copy();plane_bad=np.flatnonzero(minimum-actual>1e-10);violations+=len(plane_bad)
            if len(plane_bad):
                proper=self.planes.values(parameters)
                try:
                    self.planes.minimum=minimum.copy()
                    self.planes.minimum[plane_bad]+=np.maximum(proper[plane_bad]-actual[plane_bad],0)+2e-10
                    A,b,_=self.planes.local_rows(parameters,plane_bad);blocks.append(A);bounds.append(b)
                finally:self.planes.minimum=minimum
        A=sparse.vstack(blocks,format='csr') if blocks else sparse.csr_matrix((0,6*self.count))
        b=np.concatenate(bounds) if bounds else np.array([])
        return A,b,dict(complete=True,violatingExactConstraints=violations,packedFacetFailures=details,
            packedPlaneFailures=len(plane_bad),actualAtlasFloat32EncodingReplayed=True,originalPhysicalBudgetsUnchanged=True,
            maximumPackingError_m=max(float(np.linalg.norm(a-b,axis=1).max()) for a,b in zip(exact,packed)))
