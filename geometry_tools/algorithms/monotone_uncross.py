"""Micrometre source-core repair with exact local-to-global contact replay.

Only a strict reduction in complete intrinsic triangle crossings is retained.
Source vertices/faces, fixed controls and declared aliases remain present.
Anatomical source sides propose directions; actual triangles decide acceptance.
"""
import time

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts, self_contact_summary


def monotone_uncross(mesh, reference, physical_ids, protected, contexts=(), *,
                     maximum_motion_m=.00001, wall_seconds=18., maximum_trials=600,
                     preserve_context_counts=False, rigid_groups=(), rigid_motion_m=.001,
                     maximum_accepted_steps=None, face_owners=None, defer_full_replay=False,
                     coupled_groups=False):
    if mesh.frame != reference.frame or not np.array_equal(mesh.faces, reference.faces):
        raise ValueError('Source-side reference must retain the same frame and faces')
    if not 0 < maximum_motion_m <= .00002 or not 0 < wall_seconds <= 20:
        raise ValueError('Explicit micrometre motion and native wall bounds required')
    start = time.monotonic()
    labels=np.asarray(face_owners,int)if face_owners is not None else None
    if labels is not None and labels.shape!=(len(mesh.faces),):raise ValueError('One source owner per face required')
    q, base, faces = mesh.vertices.copy(), mesh.vertices, mesh.faces
    ids = np.asarray(physical_ids, int)
    locked = set(map(int, protected))
    aliases = {int(k): np.flatnonzero(ids == k) for k in np.unique(ids)}
    grouped=set();group_keys={}
    for group in rigid_groups:
        members=set(map(int,group))
        if grouped.intersection(members):raise ValueError('Rigid source sections overlap')
        grouped.update(members);key=min(members)
        if locked.intersection(members):locked.update(members)
        combined=np.concatenate([aliases[k]for k in sorted(members)])
        for k in members:group_keys[k]=key
        aliases[key]=combined
    for k in list(aliases):
        if k in group_keys and group_keys[k]!=k:del aliases[k]
    incident = {k: np.flatnonzero(np.isin(faces, raw).any(1)) for k, raw in aliases.items()}
    first, _ = self_contact_summary(mesh, example_limit=10000)
    if not first['complete_within_pair_budget'] or first['defect_pairs'] != len(first['examples']):
        raise ValueError('Repair requires the complete actual contact inventory')
    audit = first
    trials, steps, timed_out = 0, [], False;budget_exhausted=False
    while audit['defect_pairs']:
        if maximum_accepted_steps is not None and len(steps)>=maximum_accepted_steps:break
        if time.monotonic()-start >= wall_seconds or trials >= maximum_trials:
            timed_out = time.monotonic()-start >= wall_seconds
            budget_exhausted=trials>=maximum_trials
            break
        current = Mesh(q, faces, mesh.frame)
        tree = BVH(current)
        defect_pairs = np.array(sorted(tuple(sorted(p))for p in audit['examples']), int)
        proposals = [];direct=[]
        for pair in defect_pairs:
            for fa, fb in [pair, pair[::-1]]:
                tb, rb = q[faces[fb]], reference.vertices[faces[fb]]
                n = np.cross(tb[1]-tb[0], tb[2]-tb[0])
                nr = np.cross(rb[1]-rb[0], rb[2]-rb[0])
                if min(np.linalg.norm(n), np.linalg.norm(nr)) < 1e-20:
                    continue
                n /= np.linalg.norm(n); nr /= np.linalg.norm(nr)
                for raw in faces[fa]:
                    original_key=int(ids[raw]);key=group_keys.get(original_key,original_key)
                    if original_key in locked or np.any(ids[faces[fb]] == original_key):
                        continue
                    source_side = nr@(reference.vertices[raw]-rb[0])
                    if abs(source_side) > 1e-11:
                        direction = n*np.sign(source_side)
                        signed = direction@(q[raw]-tb[0])
                        distance = max(.00000005, .00000005-signed)
                        bound=rigid_motion_m if key in grouped else maximum_motion_m
                        if distance <= bound:
                            direct.append((distance,key,direction*distance))
                        other_originals=set(map(int,ids[faces[fb]]))
                        other_keys={group_keys.get(k,k)for k in other_originals}
                        if coupled_groups and key in grouped:
                            for other in sorted(other_keys):
                                if other!=key and other in grouped and not locked.intersection(ids[aliases[other]]):
                                    half=distance*.5
                                    if half<=rigid_motion_m:
                                        direct.append((half,(key,other),direction*half))
                    # Bounded coordinate directions cover coplanar slivers.
                    for length in [.0000001, .0000005, .000001, .000003]:
                        for direction in [n, -n]:
                            proposals.append((length, key, length*direction))
                    if key in grouped:
                        for length in [.00001,.00005,.0001,.0003]:
                            for direction in [n,-n]:proposals.append((length,key,length*direction))
        proposals.sort(key=lambda item: item[0])
        if rigid_groups:
            affected={k:int(np.isin(defect_pairs,incident[k]).any(1).sum())for k in aliases}
            for _,k,_ in direct:
                if isinstance(k,tuple)and k not in affected:
                    changed=np.unique(np.r_[incident[k[0]],incident[k[1]]])
                    affected[k]=int(np.isin(defect_pairs,changed).any(1).sum())
            rank=lambda item:(-affected[item[1]],item[0])
            proposals=sorted(direct,key=rank)+sorted(proposals,key=rank)
        else:proposals=sorted(proposals+direct,key=lambda item:item[0])
        seen = set()
        accepted = False
        for length, key, displacement in proposals:
            signature = (key, tuple(np.rint(displacement/1e-12).astype(int)))
            if signature in seen:
                continue
            seen.add(signature)
            if time.monotonic()-start >= wall_seconds or trials >= maximum_trials:
                timed_out = time.monotonic()-start >= wall_seconds
                budget_exhausted=trials>=maximum_trials
                break
            trials += 1
            paired=isinstance(key,tuple)
            raw=np.r_[aliases[key[0]],aliases[key[1]]]if paired else aliases[key]
            trial = q.copy()
            if paired:
                trial[aliases[key[0]]]+=displacement
                trial[aliases[key[1]]]-=displacement
            else:trial[raw] += displacement
            bound=rigid_motion_m if paired or key in grouped else maximum_motion_m
            if np.linalg.norm(trial[raw]-base[raw], axis=1).max() > bound+1e-15:
                continue
            changed=np.unique(np.r_[incident[key[0]],incident[key[1]]])if paired else incident[key]
            small = Mesh(trial, faces[changed], mesh.frame)
            area=np.linalg.norm(np.cross(trial[faces[changed,1]]-trial[faces[changed,0]],
                                         trial[faces[changed,2]]-trial[faces[changed,0]]),axis=1)
            original_area=np.linalg.norm(np.cross(reference.vertices[faces[changed,1]]-reference.vertices[faces[changed,0]],
                                                  reference.vertices[faces[changed,2]]-reference.vertices[faces[changed,0]]),axis=1)
            if np.any(area<=original_area*.05):continue
            contacts, proof = pair_contacts(small, current, 0., second_bvh=tree,
                                            max_tests=100000, max_pairs=10000)
            if not proof['complete_within_pair_budget']:
                continue
            wrong = np.isin(contacts['relation'], ['crossing', 'coplanar-overlap'])
            wrong &= ~np.isin(contacts['pairs'][:, 1], changed)
            local, _ = self_contact_summary(small, example_limit=10000 if labels is not None else 0, max_tests=100000)
            if not local['complete_within_pair_budget']:
                continue
            old = int(np.isin(defect_pairs, changed).any(1).sum())
            remaining = audit['defect_pairs']-old+int(wrong.sum())+local['defect_pairs']
            if remaining >= audit['defect_pairs']:
                continue
            if labels is not None:
                from collections import Counter
                before_pairs=defect_pairs[np.isin(defect_pairs,changed).any(1)]
                new_pairs=np.c_[changed[contacts['pairs'][wrong,0]],contacts['pairs'][wrong,1]]
                if local['examples']:new_pairs=np.r_[new_pairs,changed[np.array(local['examples'],int)]]
                hist=lambda pairs:Counter(tuple(sorted(labels[pair].tolist()))for pair in pairs)
                was_hist,now_hist=hist(before_pairs),hist(new_pairs)
                if any(n>was_hist[k]for k,n in now_hist.items()):continue
            body_ok = True;context_counts=[]
            for ci,context in enumerate(contexts):
                context_hit, context_proof = pair_contacts(small, context.mesh, 0., second_bvh=context,
                                                  max_tests=100000, max_pairs=10000)
                now=context_proof['crossing_face_pairs']+context_proof['coplanar_overlap_face_pairs']
                was=0
                if preserve_context_counts:
                    before_hit,before_proof=pair_contacts(Mesh(q,faces[changed],mesh.frame),context.mesh,0.,second_bvh=context,
                                                max_tests=100000,max_pairs=10000)
                    if not before_proof['complete_within_pair_budget']:body_ok=False;break
                    was=before_proof['crossing_face_pairs']+before_proof['coplanar_overlap_face_pairs']
                    if labels is not None:
                        def owner_counts(hit):
                            use=np.isin(hit['relation'],['crossing','coplanar-overlap'])
                            return Counter(labels[changed[hit['pairs'][use,0]]].tolist())
                        bh,nh=owner_counts(before_hit),owner_counts(context_hit)
                        if any(n>bh[k]for k,n in nh.items()):body_ok=False;break
                if not context_proof['complete_within_pair_budget'] or now>was:
                    body_ok = False
                    break
                context_counts.append(dict(context=ci,before=was,after=now))
            if not body_ok:
                continue
            if defer_full_replay:
                unchanged=defect_pairs[~np.isin(defect_pairs,changed).any(1)]
                fresh=np.c_[changed[contacts['pairs'][wrong,0]],contacts['pairs'][wrong,1]]
                # Retain explicit local pairs even when source-owner guards
                # are not requested. Their IDs must be complete, not sampled.
                if labels is None:
                    local,_=self_contact_summary(small,example_limit=10000,max_tests=100000)
                if local['defect_pairs']!=len(local['examples']):
                    raise ValueError('Complete local pair IDs are required for deferred replay')
                if local['examples']:fresh=np.r_[fresh,changed[np.array(local['examples'],int)]]
                maintained=np.unique(np.sort(np.r_[unchanged,fresh],axis=1),axis=0)
                if len(maintained)!=remaining:
                    raise ValueError('Local pair-set delta has duplicates or missing pairs')
                final=dict(audit,defect_pairs=remaining,examples=maintained.tolist())
            else:
                final, _ = self_contact_summary(Mesh(trial, faces, mesh.frame), example_limit=10000)
                if not final['complete_within_pair_budget'] or final['defect_pairs'] != remaining:
                    raise ValueError('Local contact delta did not replay in the complete source mesh')
            steps.append(dict(physicalControl=key, rawControls=raw.tolist(),
                              displacement_m=displacement.tolist(), before=audit['defect_pairs'], after=remaining,
                              changedFaces=changed.tolist(),contexts=context_counts,
                              rigidSection=paired or key in grouped,coupledSourceSections=paired))
            q, audit, accepted = trial, final, True
            break
        if not accepted:
            break
    if not np.array_equal(q[np.isin(ids, list(locked))], base[np.isin(ids, list(locked))]):
        raise ValueError('Protected source geometry changed')
    if defer_full_replay:
        final,_=self_contact_summary(Mesh(q,faces,mesh.frame),example_limit=10000)
        actual={tuple(sorted(p))for p in final['examples']}
        maintained={tuple(sorted(p))for p in audit['examples']}
        if not final['complete_within_pair_budget'] or final['defect_pairs']!=len(actual) or actual!=maintained:
            raise ValueError('Independent complete source replay disagrees with maintained local pairs')
        audit=final
    timed_out=timed_out or time.monotonic()-start>=wall_seconds
    return q, dict(complete=audit['defect_pairs'] == 0 and not timed_out and not budget_exhausted,
                   initialSelfPairs=first['defect_pairs'], finalSelfPairs=audit['defect_pairs'],
                   acceptedSteps=steps, trials=trials, timedOut=timed_out,
                   trialBudgetExhausted=budget_exhausted,
                   stepBudgetExhausted=maximum_accepted_steps is not None and len(steps)>=maximum_accepted_steps and audit['defect_pairs']>0,
                   maximumSourceMotion_m=float(np.linalg.norm(q-base, axis=1).max()),
                   nativeMotionBound_m=maximum_motion_m, nativeWallSeconds=wall_seconds,
                   rigidSectionMotionBound_m=rigid_motion_m if rigid_groups else None,
                   contextPolicy='nonincrease' if preserve_context_counts else 'zero',
                   deferredFullReplay=defer_full_replay,independentFullPairSetReplay=True,
                   coupledSourceSections=coupled_groups,
                   sourceFacesAndVerticesRetained=True, protectedControlsLiteral=True,
                   elapsedSeconds=time.monotonic()-start)
