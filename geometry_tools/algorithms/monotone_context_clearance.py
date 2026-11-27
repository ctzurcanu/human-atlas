"""Directed, source-group-preserving triangle clearance against a fixed mesh.

Every retained step reduces actual target crossings, or explicitly reduces
sampled skin-side violation without increasing crossings. Self/other-context
contacts never increase. Closed-target controls stay outside its volume.
It is a bounded numerical repair, not anatomical or lumen acceptance.
"""
import time
from collections import Counter
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts, self_contact_summary
from geometry_tools.spatial import BVH
from geometry_tools.open3d_backend import library


def _contacts(mesh, tree):
    hit, proof=pair_contacts(mesh,tree.mesh,0.,second_bvh=tree,max_tests=250000,max_pairs=20000)
    if not proof['complete_within_pair_budget']:
        raise ValueError('Incomplete actual target/context inventory')
    return hit['pairs'][np.isin(hit['relation'],['crossing','coplanar-overlap'])]


def monotone_context_clearance(mesh, physical, sections, protected, target, target_index,
                                contexts=(), *, wall_seconds=12., maximum_trials=300,
                                maximum_motion_m=.001, maximum_accepted_steps=6,
                                connected_contact_groups=False,
                                target_side='closed-exterior',geodesic_spread_m=0.,
                                corotate_sections=False,skin_penetration_descent=False):
    if not 0 < wall_seconds <= 18 or not 0 < maximum_motion_m <= .003 or maximum_trials<1 or not 1<=maximum_accepted_steps<=32:
        raise ValueError('Finite native time/work/motion limits required')
    if target_side not in ['closed-exterior','interior-sheet']:
        raise ValueError('Explicit closed-solid exterior or oriented skin-sheet interior required')
    if not 0<=geodesic_spread_m<=.02:
        raise ValueError('Explicit source-geodesic influence width within [0,20] mm required')
    if corotate_sections and not geodesic_spread_m:
        raise ValueError('Section corotation requires an explicit source-geodesic field')
    if skin_penetration_descent and target_side!='interior-sheet':
        raise ValueError('Skin-side descent requires an explicit oriented interior sheet')
    direction = -1. if target_side == 'interior-sheet' else 1.
    started=time.monotonic();q=mesh.vertices.copy();base=q.copy();f=mesh.faces
    ids=np.asarray(physical,int);locked=set(map(int,protected));covered=set();groups=[]
    for section in sections:
        members=set(map(int,ids[section]))
        if covered.intersection(members):raise ValueError('Source section controls overlap')
        covered.update(members);groups.append(np.flatnonzero(np.isin(ids,list(members))))
    groups.extend(np.flatnonzero(ids==k)for k in np.unique(ids)if int(k)not in covered)
    owner=np.full(len(q),-1,int)
    for k,raw in enumerate(groups):owner[raw]=k
    frozen=[bool(locked.intersection(map(int,ids[raw])))for raw in groups]
    incident=[np.flatnonzero(np.isin(f,raw).any(1))for raw in groups]
    spread_cache={}
    if geodesic_spread_m:
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import dijkstra
        centre=np.array([base[raw].mean(0)for raw in groups])
        edge=np.unique(np.sort(owner[f[:,[0,1,1,2,2,0]].reshape(-1,2)],axis=1),axis=0)
        edge=edge[edge[:,0]!=edge[:,1]]
        length=np.maximum(np.linalg.norm(centre[edge[:,0]]-centre[edge[:,1]],axis=1),1e-12)
        source_graph=coo_matrix((np.r_[length,length],(np.r_[edge[:,0],edge[:,1]],np.r_[edge[:,1],edge[:,0]])),shape=(len(groups),len(groups))).tocsr()
    pairs=_contacts(mesh,target);first=len(pairs);self_audit,_=self_contact_summary(mesh,example_limit=10000)
    if not self_audit['complete_within_pair_budget'] or self_audit['defect_pairs']!=len(self_audit['examples']):
        raise ValueError('Incomplete source self inventory')
    self_first=self_audit['defect_pairs'];steps=[];trials=0;exhausted=False;timed_out=False;rejected=Counter()
    def skin_metric(vertices):
        from .feature_normals import indexed_feature_normals,inward_sheet_clearance
        tri=vertices[f]
        samples=np.r_[vertices[np.unique(f)],tri.mean(1),
                      (tri[:,0]+tri[:,1])*.5,(tri[:,1]+tri[:,2])*.5,(tri[:,2]+tri[:,0])*.5]
        hit,_=target_index.closest_points(samples)
        normal,rim=indexed_feature_normals(target_index,hit)
        _,gap,known=inward_sheet_clearance(samples,hit,normal,rim)
        violation=np.maximum(.00002-gap[known],0.)
        return dict(squaredViolation_m2=float(violation@violation),
                    maximumViolation_m=float(violation.max(initial=0.)),
                    knownSamples=int(known.sum()),unknownRimSamples=int((~known).sum()))
    metric=skin_metric(q)if target_side=='interior-sheet'else None
    first_metric=None if metric is None else dict(metric)
    def skin_pending():
        return metric is not None and (metric['maximumViolation_m']>1e-8 or
               metric['unknownRimSamples']>0 or metric['knownSamples']==0)
    while len(pairs) or skin_pending():
        if len(steps)>=maximum_accepted_steps:break
        if time.monotonic()-started>=wall_seconds or trials>=maximum_trials:
            timed_out=time.monotonic()-started>=wall_seconds;exhausted=trials>=maximum_trials;break
        current=Mesh(q,f,mesh.frame);tree=BVH(current);proposals=[]
        if connected_contact_groups:
            # Move a contacted source patch together. Moving only one triangle
            # cannot clear a patch whose adjacent controls are still inside the
            # fixed solid; its full incident-face occupancy check must reject it.
            parent={}
            def root(k):
                parent.setdefault(k,k)
                while parent[k]!=k:
                    parent[k]=parent[parent[k]];k=parent[k]
                return k
            for face in np.unique(pairs[:,0]):
                keys=sorted(set(map(int,owner[f[face]])))
                for k in keys[1:]:parent[root(k)]=root(keys[0])
            patches={}
            for k in parent:patches.setdefault(root(k),[]).append(k)
            for keys in patches.values():
                if len(keys)<2 or any(frozen[k]for k in keys):continue
                raw=np.unique(np.concatenate([groups[k]for k in keys]))
                changed=np.unique(np.concatenate([incident[k]for k in keys]))
                selected=pairs[np.isin(pairs[:,0],changed)]
                triangles=target.mesh.triangles[np.unique(selected[:,1])]
                normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
                lengths=np.linalg.norm(normals,axis=1)
                normals=direction*normals[lengths>1e-20]/lengths[lengths>1e-20,None]
                mean=normals.sum(0)
                directions=([mean/np.linalg.norm(mean)]if np.linalg.norm(mean)>1e-12 else [])
                directions.extend(np.unique(np.round(normals,6),axis=0)[:64])
                for n in directions:
                    n=n/np.linalg.norm(n)
                    needed=max(.00001,.00002-float(((q[raw]-triangles.mean((0,1)))@n).min()))
                    for distance in sorted(set([needed,.00005,.0001,.00025,.0005])):
                        if distance<=maximum_motion_m:
                            proposals.append((-len(selected),distance,tuple(sorted(keys)),n*distance))
        for va,fb in pairs:
            tri=target.mesh.triangles[fb];n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n*=direction/np.linalg.norm(n)
            signed=(q[f[va]]-tri[0])@n
            needed=max(.00001,.00002-float(signed.min()))
            active=[int(k)for k in np.unique(owner[f[va]])if not frozen[int(k)]]
            for keys in [(k,)for k in active]+([tuple(active)]if len(active)>1 else[]):
                if not keys:continue
                changed=np.unique(np.concatenate([incident[k]for k in keys]))
                affected=int(np.isin(pairs[:,0],changed).sum())
                for distance in sorted(set([needed,.00005,.0001,.00025,.0005])):
                    if distance<=maximum_motion_m:
                        proposals.append((-affected,distance,keys,n*distance))
        proposals.sort(key=lambda p:(p[0],-p[1]if skin_penetration_descent else p[1]));seen=set();accepted=False
        for _,distance,keys,d in proposals:
            sig=(keys,tuple(np.rint(d/1e-12).astype(int)))
            if sig in seen:continue
            seen.add(sig)
            if time.monotonic()-started>=wall_seconds or trials>=maximum_trials:
                timed_out=time.monotonic()-started>=wall_seconds;exhausted=trials>=maximum_trials;break
            trials+=1
            if geodesic_spread_m:
                if keys not in spread_cache:
                    distance=dijkstra(source_graph,indices=list(keys),min_only=True,directed=False)
                    t=np.clip(distance/geodesic_spread_m,0,1)
                    weight=1-3*t*t+2*t*t*t;weight[frozen]=0.
                    spread_cache[keys]=weight
                weight=spread_cache[keys];moving=np.flatnonzero(weight>0)
                raw=np.unique(np.concatenate([groups[k]for k in moving]));trial=q.copy()
                trial[raw]+=weight[owner[raw],None]*d
                if corotate_sections:
                    valid=True
                    for k in moving:
                        neighbours=source_graph.indices[source_graph.indptr[k]:source_graph.indptr[k+1]]
                        if not len(neighbours):continue
                        gradient=np.linalg.lstsq(centre[neighbours]-centre[k],weight[neighbours]-weight[k],rcond=1e-10)[0]
                        F=np.eye(3)+np.outer(d,gradient)
                        if np.linalg.det(F)<=.1:
                            valid=False;break
                        U,_,Vt=np.linalg.svd(F);R=U@Vt
                        local=groups[k];c=q[local].mean(0)
                        trial[local]=(q[local]-c)@R.T+c+weight[k]*d
                    if not valid:
                        rejected['nonpositive-local-carrier']+=1;continue
            else:
                moving=keys;raw=np.unique(np.concatenate([groups[k]for k in keys]));trial=q.copy();trial[raw]+=d
            if np.linalg.norm(trial[raw]-base[raw],axis=1).max()>maximum_motion_m+1e-15:
                rejected['motion-bound']+=1;continue
            changed=np.unique(np.concatenate([incident[k]for k in moving]));small=Mesh(trial,f[changed],mesh.frame)
            # Reject collapsed or reversed source faces before any expensive
            # contact replay. Group congruence alone does not protect the
            # connecting charts between two translated source sections.
            old_tri=base[f[changed]];new_tri=small.triangles
            old_n=np.cross(old_tri[:,1]-old_tri[:,0],old_tri[:,2]-old_tri[:,0])
            new_n=np.cross(new_tri[:,1]-new_tri[:,0],new_tri[:,2]-new_tri[:,0])
            area=np.linalg.norm(old_n,axis=1);new_area=np.linalg.norm(new_n,axis=1)
            if np.any(new_area<area*.5) or np.any(np.einsum('ij,ij->i',old_n,new_n)<=0):
                rejected['source-face-orientation-or-area']+=1;continue
            now=_contacts(small,target);old=int(np.isin(pairs[:,0],changed).sum());remaining=len(pairs)-old+len(now)
            if remaining>len(pairs) or (remaining==len(pairs)and not skin_penetration_descent):
                rejected['target-crossings-not-reduced']+=1;continue
            next_metric=None
            if skin_penetration_descent:
                next_metric=skin_metric(trial)
                if (next_metric['squaredViolation_m2']>=metric['squaredViolation_m2']-1e-16 or
                    next_metric['maximumViolation_m']>metric['maximumViolation_m']+1e-8 or
                    next_metric['unknownRimSamples']>metric['unknownRimSamples']):
                    rejected['skin-violation-not-reduced']+=1;continue
            sample=trial[np.unique(f[changed])]
            if target_side == 'interior-sheet'and not skin_penetration_descent:
                from .feature_normals import indexed_feature_normals,inward_sheet_clearance
                hit,_=target_index.closest_points(sample)
                normal,rim=indexed_feature_normals(target_index,hit)
                _,gap,_=inward_sheet_clearance(sample,hit,normal,rim)
                if np.any(gap < .00002-1e-8):
                    rejected['skin-interior-clearance']+=1;continue
            elif target_side=='closed-exterior':
                tensor=library().core.Tensor(((sample-target_index.origin)/target_index.scale).astype(np.float32))
                if (target_index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5).any():
                    rejected['closed-target-occupancy']+=1;continue
            hit,proof=pair_contacts(small,current,0.,second_bvh=tree,max_tests=250000,max_pairs=20000)
            if not proof['complete_within_pair_budget']:continue
            external=np.isin(hit['relation'],['crossing','coplanar-overlap'])&~np.isin(hit['pairs'][:,1],changed)
            local,_=self_contact_summary(small,example_limit=0)
            if not local['complete_within_pair_budget']:continue
            old_self=int(np.isin(np.array(self_audit['examples'],int),changed).any(1).sum()) if self_audit['defect_pairs'] else 0
            new_self=self_audit['defect_pairs']-old_self+int(external.sum())+local['defect_pairs']
            if new_self>self_audit['defect_pairs']:
                rejected['new-source-self-crossings']+=1;continue
            good=True;counts=[]
            for ci,context in enumerate(contexts):
                was=len(_contacts(Mesh(q,f[changed],mesh.frame),context));count=len(_contacts(small,context))
                counts.append(dict(context=ci,before=was,after=count))
                if count>was:good=False;break
            if not good:
                rejected['recipient-crossing-regression']+=1;continue
            full=Mesh(trial,f,mesh.frame);next_pairs=_contacts(full,target)
            next_self,_=self_contact_summary(full,example_limit=10000)
            if len(next_pairs)!=remaining or not next_self['complete_within_pair_budget'] or next_self['defect_pairs']!=new_self:
                raise ValueError('Local source-group contact deltas did not replay globally')
            step=dict(groups=list(map(int,moving)),seedGroups=list(keys),before=len(pairs),after=remaining,displacement_m=d.tolist(),
                      selfBefore=self_audit['defect_pairs'],selfAfter=new_self,contexts=counts)
            if geodesic_spread_m:step['groupMotionWeights']=weight[moving].tolist()
            if metric is not None:
                if next_metric is None:next_metric=skin_metric(trial)
                step['sampledSkinViolationBefore']=dict(metric)
                step['sampledSkinViolationAfter']=dict(next_metric)
                metric=next_metric
            steps.append(step)
            q,pairs,self_audit=trial,next_pairs,next_self;accepted=True;break
        if not accepted:break
    elapsed=time.monotonic()-started;timed_out=timed_out or elapsed>=wall_seconds
    return q,dict(complete=len(pairs)==0 and not skin_pending() and not timed_out and not exhausted,
                  initialTargetPairs=first,finalTargetPairs=len(pairs),initialSelfPairs=self_first,
                  finalSelfPairs=self_audit['defect_pairs'],acceptedSteps=steps,trials=trials,
                  timedOut=timed_out,trialBudgetExhausted=exhausted,motionBound_m=maximum_motion_m,
                  stepBudgetExhausted=len(steps)>=maximum_accepted_steps and len(pairs)>0,
                  maximumMotion_m=float(np.linalg.norm(q-base,axis=1).max()),
                  sourceFacesAndSectionsRetained=True,nativeWallSeconds=wall_seconds,
                  connectedContactGroups=connected_contact_groups,
                  targetSide=target_side,sourceFaceOrientationGuard=True,
                  sourceGeodesicSpread_m=geodesic_spread_m,
                  sectionCorotation=corotate_sections,
                  skinPenetrationDescent=skin_penetration_descent,
                  initialSampledSkinViolation=first_metric,finalSampledSkinViolation=metric,
                  rejectionCounts=dict(rejected),
                  openSkinContainmentEstablished=False,
                  elapsedSeconds=elapsed)
