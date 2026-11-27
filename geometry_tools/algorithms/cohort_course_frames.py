"""Align existing rigid source sections with their actual posed course.

Every section centre stays fixed. Shared physical controls move once, fixed
baseline controls stay literal, and original branch cores solve a harmonic
boundary problem. No face, branch, or recipient tissue is replaced.
"""
import numpy as np
from .port_section_graph import shortest_rotation
from .profile_core_transport import transport_profile_cores


def align_cohort_course_frames(source,current,faces,sections,baseline):
    v=np.asarray(source,float);q=np.asarray(current,float).copy();f=np.asarray(faces,int);base=np.asarray(baseline,float)
    if v.shape!=q.shape or base.shape!=q.shape:raise ValueError('Corresponding original, posed and fixed baseline controls required')
    owner=np.full(len(v),-1,int)
    for k,s in enumerate(sections):
        if np.any(owner[s]>=0):raise ValueError('Disjoint complete source sections required')
        owner[s]=k
    count=len(sections)
    for i in np.flatnonzero(owner<0):owner[i]=count;count+=1
    source_centres=np.array([v[owner==k].mean(0)for k in range(count)])
    target_centres=np.array([q[owner==k].mean(0)for k in range(count)])
    edges=np.unique(np.sort(owner[f[:,[0,1,1,2,2,0]].reshape(-1,2)],axis=1),axis=0)
    edges=edges[edges[:,0]!=edges[:,1]];neighbours=[[]for _ in range(count)]
    for a,b in edges:neighbours[a].append(int(b));neighbours[b].append(int(a))
    fixed=np.all(q==base,axis=1);records=[]
    for k,s in enumerate(sections):
        near=neighbours[k]
        if len(near)!=2 or fixed[s].any():continue
        dx=source_centres[near[1]]-source_centres[near[0]]
        dy=target_centres[near[1]]-target_centres[near[0]]
        if min(np.linalg.norm(dx),np.linalg.norm(dy))<1e-8:continue
        a=v[s]-source_centres[k];b=q[s]-target_centres[k]
        u,_,vh=np.linalg.svd(a.T@b);sign=np.eye(3);sign[-1,-1]=np.linalg.det(vh.T@u.T);R=vh.T@sign@u.T
        old=R@dx;cosine=float(old@dy/np.linalg.norm(old)/np.linalg.norm(dy))
        C=shortest_rotation(old,dy)
        q[s]=b@C.T+target_centres[k]
        records.append(dict(section=k,tangentCosineBefore=cosine))
    boundary=np.unique(np.r_[np.concatenate(sections),np.flatnonzero(fixed)])
    result=transport_profile_cores(v,f,[],np.empty((0,3)),np.empty((0,3,3)),
        physical_ids=np.arange(len(v)),prescribed_vertices=boundary,prescribed_positions=q[boundary])
    q=result['vertices'];q[fixed]=base[fixed]
    centre_error=max(float(np.linalg.norm(q[s].mean(0)-target_centres[k]))for k,s in enumerate(sections))
    return q,dict(method='Proper section rotation follows the actual posed source-incidence course; harmonic original cores',
        alignedSections=len(records),reversedTangentsBefore=sum(x['tangentCosineBefore']<0 for x in records),
        minimumTangentCosineBefore=min((x['tangentCosineBefore']for x in records),default=1.),
        maximumCentreMotion_m=centre_error,fixedOriginalControls=int(fixed.sum()),
        sourceFacesRemoved=0,skeletonChanged=False,anatomicalAcceptance=False,
        actualTriangleAndContextAuditRequired=True)
