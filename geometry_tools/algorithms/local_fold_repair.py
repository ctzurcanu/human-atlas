"""Bounded, fixed-topology fold repair using complete source contact witnesses."""
import time
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.commands.filter_hand_components import encoded_attribute
from .local_crossing_guard import LocalCrossingGuard

def repair_folds(mesh,pairs,*,maximum_motion_m=.0005,wall_seconds=4.,locked=(),target_pairs=None,face_owners=None):
    start=time.monotonic();f=mesh.faces;guard=LocalCrossingGuard(mesh.vertices,f,pairs,maximum_displacement_m=maximum_motion_m,maximum_pair_tests=800000);fixed=set(map(int,locked));incident=[guard.incident[guard.starts[i]:guard.starts[i+1]]for i in range(len(mesh.vertices))];t=mesh.triangles;fn0=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);area0=np.linalg.norm(fn0,axis=1)
    targets=guard.remaining if target_pairs is None else frozenset(tuple(sorted(map(int,p)))for p in target_pairs)
    if not targets<=guard.remaining:raise ValueError('Every target must belong to the complete source defect inventory')
    owners=None if face_owners is None else np.asarray(face_owners)
    if owners is not None and owners.shape!=(len(f),):raise ValueError('One explicit assembly owner per source face required')
    companion_moves=0
    for _ in range(min(48,2*len(targets)+2)):
        if not targets.intersection(guard.remaining):break
        best=None
        def consider(ids,values,with_companion=False):
            nonlocal best
            if fixed.intersection(ids):return
            motion=float(np.linalg.norm(values-mesh.vertices[ids],axis=1).max())
            if motion>maximum_motion_m:return
            touched=np.unique(np.concatenate([incident[v]for v in ids]));tri=guard.positions[f[touched]].copy()
            for v,q in zip(ids,values):tri[f[touched]==v]=q
            fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
            if np.any(np.linalg.norm(fn,axis=1)<.15*area0[touched]):return
            proposal=guard.check_move(ids,values)
            if proposal is None:
                if owners is None or with_companion:return
                own=np.unique(owners[touched])
                if len(own)!=1:return
                other={j for pair in guard.last_new_pairs for j in pair if owners[j]!=own[0]}
                if not other:return
                extra=np.setdiff1d(np.unique(f[sorted(other)]),ids)
                if not len(extra)or len(extra)+len(ids)>32:return
                delta=np.mean(values-guard.positions[ids],axis=0)
                consider(np.r_[ids,extra],np.concatenate([values,guard.positions[extra]+delta]),True)
                return
            remaining=len(targets.intersection(proposal.remaining))
            if remaining>=len(targets.intersection(guard.remaining)):return
            rank=(remaining,len(proposal.remaining),motion)
            if best is None or rank<best[0]:best=(rank,proposal,with_companion)
        for a,b in sorted(targets.intersection(guard.remaining)):
            if time.monotonic()-start>wall_seconds:raise TimeoutError('Local fold repair native deadline')
            for own,other in [(a,b),(b,a)]:
                tri=guard.positions[f[other]];normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);normal/=np.linalg.norm(normal);d=(guard.positions[f[own]]-tri[0])@normal
                for k,v in enumerate(f[own]):
                    rest=np.delete(d,k)
                    if rest[0]*rest[1]<=0:continue
                    p=guard.positions[v]+(np.sign(rest[0])*2e-6-d[k])*normal;consider(np.array([v]),p[None])
                for amount in [2e-6-d.min(),-2e-6-d.max()]:consider(f[own],guard.positions[f[own]]+amount*normal)
        if best is None:break
        guard.commit(best[1]);companion_moves+=int(best[2])
    if targets.intersection(guard.remaining):raise ValueError('Unresolved bounded source folds: '+str(dict(guard.evidence(),remaining_target_pairs=len(targets.intersection(guard.remaining)))))
    q=guard.positions;changed=np.flatnonzero(np.any(q!=mesh.vertices,axis=1));touch=np.unique(f[np.any(np.isin(f,changed),axis=1)]);tri=q[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);old=np.zeros_like(q);new=np.zeros_like(q)
    for k in range(3):np.add.at(old,f[:,k],fn0);np.add.at(new,f[:,k],fn)
    old/=np.linalg.norm(old,axis=1)[:,None];new/=np.linalg.norm(new,axis=1)[:,None];axis=np.cross(old[touch],new[touch]);cos=np.einsum('ij,ij->i',old[touch],new[touch]);attrs={k:v.copy()for k,v in mesh.attributes.items()};words=encoded_attribute('normals',attrs['normals'],np.dtype('<i2'));n=words[touch].astype(float)*attrs['normalSigns'][touch,None];n/=np.linalg.norm(n,axis=1)[:,None]
    if np.any(1+cos<1e-8):raise ValueError('Antiparallel vertex normal star requires separate review')
    n=n+np.cross(axis,n)+np.cross(axis,np.cross(axis,n))/(1+cos[:,None]);n/=np.linalg.norm(n,axis=1)[:,None];words[touch]=np.rint(n*32767).astype('<i2');attrs['normals']=words;attrs['normalSigns'][touch]=1
    return Mesh(q,f,mesh.frame,attrs,mesh.source),dict(complete=True,changedVertices=changed.tolist(),maximumMotion_m=float(np.linalg.norm(q-mesh.vertices,axis=1).max(initial=0)),facesRetained=True,guard=guard.evidence(),targetPairCount=len(targets),remainingTargetPairs=0,allSourceDefectsRemoved=not guard.remaining,coupledOwnerMoves=companion_moves,elapsedSeconds=time.monotonic()-start)
