"""Repair inverted source annuli by congruent section poses, monotonically.

Mouths stay fixed; no source coordinate, branch or face is discarded. Only
changes reducing complete intrinsic contacts without increasing any declared
context contact count are retained. This is not anatomical acceptance.
"""
import time

import numpy as np
from scipy.spatial.transform import Rotation

from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts, self_contact_summary
from .port_section_graph import shortest_rotation


def section_normal(vertices):
    x=vertices-vertices.mean(0)
    n=np.cross(x,np.roll(x,-1,axis=0)).sum(0)
    return n/np.linalg.norm(n)


def contact_count(mesh, tree):
    _, proof=pair_contacts(mesh,tree.mesh,0.,second_bvh=tree,max_tests=150000,max_pairs=20000)
    if not proof['complete_within_pair_budget']:
        raise ValueError('Incomplete context triangle test')
    return proof['crossing_face_pairs']+proof['coplanar_overlap_face_pairs']


def repair_sections(mesh,reference,physical,profiles,strips,mouths,contexts,*,wall_seconds=18.):
    start=time.monotonic(); q=mesh.vertices.copy(); f=mesh.faces; fixed=set(map(int,mouths))
    groups=[np.flatnonzero(np.isin(physical,physical[s]))for s in profiles]
    locked=[bool(fixed.intersection(map(int,physical[s])))for s in profiles]
    audit,_=self_contact_summary(mesh,example_limit=10000)
    if not audit['complete_within_pair_budget']:raise ValueError('Incomplete source self test')
    initial=audit['defect_pairs'];steps=[];trials=0;timed_out=False
    for iteration in range(12):
        if not audit['defect_pairs']:break
        if time.monotonic()-start>wall_seconds:timed_out=True;break
        current=Mesh(q,f,mesh.frame);tree=BVH(current);pairs=np.array(audit['examples'],int)
        proposals=[]
        for strip in strips:
            a,b=strip['rings']
            for i,j in [(a,b),(b,a)]:
                s,t=profiles[i],profiles[j]
                n0=section_normal(reference.vertices[s]);n=section_normal(q[s])
                d0=reference.vertices[t].mean(0)-reference.vertices[s].mean(0)
                d=q[t].mean(0)-q[s].mean(0)
                before=float(n0@d0);now=float(n@d)
                if abs(before)<.5*np.linalg.norm(d0) or before*now>=0:continue
                desired=np.sign(before)*d
                if not locked[i]:
                    full=Rotation.from_matrix(shortest_rotation(n,desired)).as_rotvec()
                    for fraction in [1.,.5,.25]:
                        R=Rotation.from_rotvec(full*fraction).as_matrix();c=q[s].mean(0)
                        target=(q[groups[i]]-c)@R.T+c
                        proposals.append((i,target,'section-frame'))
                needed=max(.00001,abs(before)*.25-np.sign(before)*now)
                movement=np.sign(before)*n*min(needed,.003)
                for node,delta in [(j,movement),(i,-movement)]:
                    if not locked[node]:proposals.append((node,q[groups[node]]+delta,'section-advance'))
        accepted=False
        for node,target,kind in proposals:
            if time.monotonic()-start>wall_seconds:timed_out=True;break
            trials+=1;raw=groups[node];changed=np.flatnonzero(np.isin(f,raw).any(1))
            trial=q.copy();trial[raw]=target;small=Mesh(trial,f[changed],mesh.frame)
            hit,proof=pair_contacts(small,current,0.,second_bvh=tree,max_tests=150000,max_pairs=20000)
            if not proof['complete_within_pair_budget']:continue
            wrong=np.isin(hit['relation'],['crossing','coplanar-overlap'])&~np.isin(hit['pairs'][:,1],changed)
            local,_=self_contact_summary(small,example_limit=0)
            old=int(np.isin(pairs,changed).any(1).sum())
            remaining=audit['defect_pairs']-old+int(wrong.sum())+local['defect_pairs']
            if remaining>=audit['defect_pairs']:continue
            previous=Mesh(q,f[changed],mesh.frame);context_counts=[];good=True
            for name,context in contexts:
                was,now=contact_count(previous,context),contact_count(small,context)
                context_counts.append(dict(name=name,before=was,after=now))
                if now>was:good=False;break
            if not good:continue
            final,_=self_contact_summary(Mesh(trial,f,mesh.frame),example_limit=10000)
            if not final['complete_within_pair_budget']or final['defect_pairs']!=remaining:
                raise ValueError('Section-local contact result failed full replay')
            steps.append(dict(section=node,kind=kind,before=audit['defect_pairs'],after=remaining,
                              maximumMotion_m=float(np.linalg.norm(target-q[raw],axis=1).max()),contexts=context_counts))
            q,audit,accepted=trial,final,True;break
        if not accepted:break
    return q,dict(complete=audit['defect_pairs']==0 and not timed_out,initialSelfPairs=initial,
                  finalSelfPairs=audit['defect_pairs'],steps=steps,trials=trials,timedOut=timed_out,
                  mouthsLiteral=True,sectionsCongruent=True,sourceFacesRetained=True,
                  elapsedSeconds=time.monotonic()-start)
