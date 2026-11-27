"""Reject new global crossings during bounded fixed-topology deformations.

The caller supplies a complete independently audited source defect list.
An immutable source BVH with a twice-displacement margin covers every pair
that can intersect after either triangle moves. Every affected pair is tested
on the actual proposed coordinates, not just the original defect worklist.
Topology, normals, occupied volume and anatomy remain separate obligations.
"""
from dataclasses import dataclass
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH,box_distance_squared
from geometry_tools.distance import triangle_pairs

BAD=['crossing','coplanar-overlap']


class CrossingBudgetExceeded(RuntimeError):pass


@dataclass(frozen=True)
class CrossingProposal:
    owner:object
    revision:int
    vertices:np.ndarray
    values:np.ndarray
    remaining:frozenset


@dataclass(frozen=True)
class FaceChangeProposal:
    owner:object
    revision:int
    ids:np.ndarray
    faces:np.ndarray
    active:np.ndarray
    remaining:frozenset


class LocalCrossingGuard:
    def __init__(self,positions,faces,complete_defect_face_pairs,*,maximum_displacement_m,bvh=None,maximum_pair_tests=2000000):
        self.positions=np.asarray(positions,float).copy();self.faces=np.asarray(faces,dtype=np.int64).copy()
        self.maximum=float(maximum_displacement_m)
        if not np.isfinite(self.maximum)or self.maximum<=0:raise ValueError('Finite positive source displacement bound required')
        if not isinstance(maximum_pair_tests,int)or maximum_pair_tests<1:raise ValueError('Finite positive pair-test budget required')
        source=Mesh(self.positions,self.faces,'bounded-crossing-source')
        if bvh is not None and(not np.array_equal(bvh.mesh.vertices,self.positions)or not np.array_equal(bvh.mesh.faces,self.faces)):raise ValueError('BVH must contain the literal source vertices and faces')
        self.tree=bvh or BVH(source);self.source=self.tree.mesh.vertices
        self.source_triangles=self.tree.triangles
        self.low=self.source_triangles.min(1);self.high=self.source_triangles.max(1)
        pairs=np.asarray(complete_defect_face_pairs,dtype=np.int64).reshape(-1,2)
        if len(pairs)and(pairs.min()<0 or pairs.max()>=len(self.faces)or np.any(pairs[:,0]==pairs[:,1])):raise ValueError('Distinct existing defect face pairs required')
        self.remaining=frozenset(tuple(map(int,p))for p in np.unique(np.sort(pairs,axis=1),axis=0))
        if len(pairs)and not np.isin(triangle_pairs(self.source_triangles[pairs[:,0]],self.source_triangles[pairs[:,1]])['relation'],BAD).all():raise ValueError('Every supplied source defect must replay on the actual triangles')
        self.faces.flags.writeable=False;self.active=np.ones(len(self.faces),bool);self.topology_changed=False
        self.initial=len(self.remaining);self.revision=0;self.owner=object();self.cache={}
        order=np.argsort(self.faces.ravel(),kind='stable');counts=np.bincount(self.faces.ravel(),minlength=len(self.positions))
        self.starts=np.r_[0,np.cumsum(counts)];self.incident=order//3
        self.maximum_pair_tests=maximum_pair_tests;self.tests=0;self.proposals=0;self.rejected_new=0;self.last_new_pairs=[]

    def pairs_for(self,vertices):
        affected=np.unique(np.concatenate([self.incident[self.starts[v]:self.starts[v+1]]for v in vertices]))
        return self.pairs_for_faces(affected)

    def pairs_for_faces(self,affected):
        key=tuple(map(int,affected));cached=self.cache.get(key)
        if cached is not None:return cached
        margin=float(np.nextafter(2*self.maximum,np.inf));found=[]
        for face in affected:
            near=self.tree.query_bounds(self.low[face],self.high[face],distance_m=margin)
            near=near[near!=face]
            found.extend((min(int(face),int(other)),max(int(face),int(other)))for other in near)
        pairs=np.asarray(sorted(set(found)),dtype=np.int64).reshape(-1,2)
        self.cache[key]=pairs;return pairs

    def check_move(self,vertices,values):
        if self.topology_changed:raise ValueError('Rebuild the vertex incidence index before moving a retriangulated surface')
        ids=np.asarray(vertices,dtype=np.int64);values=np.asarray(values,float)
        if ids.ndim!=1 or not 1<=len(ids)<=32 or len(np.unique(ids))!=len(ids)or ids.min()<0 or ids.max()>=len(self.positions)or values.shape!=(len(ids),3)or not np.isfinite(values).all():raise ValueError('One to 32 distinct source vertices and finite proposed positions required')
        if np.linalg.norm(values-self.source[ids],axis=1).max()>self.maximum:raise ValueError('Proposal exceeds the declared immutable-BVH displacement bound')
        pairs=self.pairs_for(ids)
        if self.tests+len(pairs)>self.maximum_pair_tests:raise CrossingBudgetExceeded('Global crossing guard pair-test budget exhausted')
        triangles=self.positions[self.faces[pairs]].copy()
        for vertex,value in zip(ids,values):triangles[self.faces[pairs]==vertex]=value
        overlap=box_distance_squared(triangles[:,0].min(1),triangles[:,0].max(1),triangles[:,1].min(1),triangles[:,1].max(1))<=0
        bad=np.zeros(len(pairs),bool)
        if overlap.any():bad[overlap]=np.isin(triangle_pairs(triangles[overlap,0],triangles[overlap,1])['relation'],BAD)
        self.tests+=int(overlap.sum());self.proposals+=1
        after=frozenset(tuple(map(int,p))for p in pairs[bad]);covered=frozenset(tuple(map(int,p))for p in pairs)
        if not after<=self.remaining:self.rejected_new+=1;self.last_new_pairs=sorted(after-self.remaining);return None
        ids=ids.copy();values=values.copy();ids.flags.writeable=False;values.flags.writeable=False
        return CrossingProposal(self.owner,self.revision,ids,values,(self.remaining-covered)|after)

    def check_faces(self,ids,faces,active):
        ids=np.asarray(ids,dtype=np.int64);faces=np.asarray(faces,dtype=np.int64);active=np.asarray(active,bool)
        if ids.ndim!=1 or not 0<len(ids)<=10000 or len(np.unique(ids))!=len(ids)or ids.min()<0 or ids.max()>=len(self.faces)or faces.shape!=(len(ids),3)or active.shape!=(len(ids),)or faces.min()<0 or faces.max()>=len(self.positions):raise ValueError('A bounded existing face patch, corresponding triangles and active mask are required')
        order=np.argsort(ids);ids,faces,active=ids[order],faces[order],active[order]
        if np.any(active&~self.active[ids]):raise ValueError('Removed source faces cannot be reactivated')
        if np.linalg.norm(self.positions[faces[active]]-self.source_triangles[ids[active]],axis=2).max(initial=0)>self.maximum:raise ValueError('Retriangulation exceeds the immutable-BVH corner displacement bound')
        pairs=self.pairs_for_faces(ids)
        if self.tests+len(pairs)>self.maximum_pair_tests:raise CrossingBudgetExceeded('Global crossing guard pair-test budget exhausted')
        current=self.faces[pairs].copy();live=self.active[pairs].copy()
        lookup=np.searchsorted(ids,pairs);inside=lookup<len(ids);match=inside&(ids[np.minimum(lookup,len(ids)-1)]==pairs)
        current[match]=faces[lookup[match]];live[match]=active[lookup[match]]
        eligible=live.all(1);tri=self.positions[current[eligible]]
        overlap=box_distance_squared(tri[:,0].min(1),tri[:,0].max(1),tri[:,1].min(1),tri[:,1].max(1))<=0
        bad=np.zeros(len(pairs),bool)
        if overlap.any():bad[np.flatnonzero(eligible)[overlap]]=np.isin(triangle_pairs(tri[overlap,0],tri[overlap,1])['relation'],BAD)
        self.tests+=int(overlap.sum());self.proposals+=1
        after=frozenset(tuple(map(int,p))for p in pairs[bad]);covered=frozenset(tuple(map(int,p))for p in pairs)
        if not after<=self.remaining:self.rejected_new+=1;self.last_new_pairs=sorted(after-self.remaining);return None
        ids=ids.copy();faces=faces.copy();active=active.copy()
        for array in [ids,faces,active]:array.flags.writeable=False
        return FaceChangeProposal(self.owner,self.revision,ids,faces,active,(self.remaining-covered)|after)

    def commit(self,proposal):
        if not isinstance(proposal,(CrossingProposal,FaceChangeProposal))or proposal.owner is not self.owner or proposal.revision!=self.revision:raise ValueError('Only a current proposal from this crossing guard may be committed')
        if not proposal.remaining<=self.remaining:raise ValueError('A proposal cannot add a global defect')
        if isinstance(proposal,CrossingProposal):self.positions[proposal.vertices]=proposal.values
        else:
            self.faces.flags.writeable=True;self.faces[proposal.ids]=proposal.faces;self.faces.flags.writeable=False
            self.active[proposal.ids]=proposal.active;self.topology_changed=True
        self.remaining=proposal.remaining;self.revision+=1

    def evidence(self):
        return dict(initial_global_defect_pairs=self.initial,remaining_global_defect_pairs=sorted(self.remaining),
            remaining_global_defect_pair_count=len(self.remaining),accepted_moves=self.revision,
            actual_pair_tests=self.tests,tested_proposals=self.proposals,rejected_new_crossings=self.rejected_new,
            removed_faces=int((~self.active).sum()),
            maximum_pair_tests=self.maximum_pair_tests,maximum_displacement_from_indexed_source_m=self.maximum,
            complete_source_defect_list_required=True,new_global_defects_permitted=False)
