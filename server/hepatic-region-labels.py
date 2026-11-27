"""Remove point-pinched face selections on an already closed source surface.

Only labels change. No source point, normal or triangle is reconstructed.
The lexicographic (pinched vertex fans, disagreeing edges) energy decreases.
"""
import numpy as np


def regularize_labels(vertices,faces,labels,max_passes=40):
    _,inv=np.unique(vertices,axis=0,return_inverse=True);f=inv[faces];labels=np.asarray(labels).copy();original=labels.copy();n=len(f)
    edges=np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1)
    _,ei,count=np.unique(edges,axis=0,return_inverse=True,return_counts=True)
    if np.any(count!=2):raise ValueError('Selection carrier must already be a closed two-manifold')
    owners=np.argsort(ei,kind='stable').reshape(-1,2)%n
    around=[[]for _ in range(int(f.max())+1)]
    for i,t in enumerate(f):
        for v in t:around[v].append(i)
    neighbor=np.empty((n,3),int)
    for k in range(3):
        pairs=owners[ei[k*n:(k+1)*n]];neighbor[:,k]=np.where(pairs[:,0]==np.arange(n),pairs[:,1],pairs[:,0])
    fans=[]
    for v,incident in enumerate(around):
        if not incident:fans.append(np.array([],int));continue
        order=[incident[0]];previous=-1
        while True:
            i=order[-1];j=int(np.flatnonzero(f[i]==v)[0]);candidates=[int(neighbor[i,j]),int(neighbor[i,(j+2)%3])];nxt=next(x for x in candidates if x!=previous)
            if nxt==order[0]:break
            if nxt in order:raise ValueError('Source has an invalid closed vertex fan')
            order.append(nxt);previous=i
        if len(order)!=len(incident):raise ValueError('Source has disconnected vertex fans')
        fans.append(np.asarray(order,int))
    def score(v):
        c=labels[fans[v]]
        if len(c)==0:return 0
        starts=c[c!=np.roll(c,1)]
        return len(starts)-len(np.unique(starts))
    scores=np.array([score(v)for v in range(len(fans))]);before=int(scores.sum());counts=np.bincount(labels);moves=0
    for iteration in range(max_passes):
        bad=np.flatnonzero(scores)
        if not len(bad):break
        changed=False
        for v in bad:
            if not scores[v]:continue
            best=None
            for face in fans[v]:
                old=int(labels[face])
                if counts[old]<=1:continue
                touched=f[face];before_local=int(scores[touched].sum());prior_edges=int(np.count_nonzero(labels[neighbor[face]]!=old))
                for color in np.unique(labels[neighbor[face]]):
                    color=int(color)
                    if color==old:continue
                    labels[face]=color;after_local=[score(x)for x in touched];delta=sum(after_local)-before_local;edge_delta=int(np.count_nonzero(labels[neighbor[face]]!=color))-prior_edges;labels[face]=old
                    if delta<0 or delta==0 and edge_delta<0:
                        key=(delta,edge_delta,int(face),color)
                        if best is None or key<best[0]:best=(key,int(face),old,color,after_local)
            if best:
                _,face,old,color,after_local=best;labels[face]=color;counts[old]-=1;counts[color]+=1;scores[f[face]]=after_local;moves+=1;changed=True
        if not changed:
            # A boundary sometimes needs two adjacent triangles to move
            # together; either individual move would create another pinch.
            best=None;seen=set()
            for v in np.flatnonzero(scores):
                for first in fans[v]:
                    for second in neighbor[first]:
                        pair=tuple(sorted((int(first),int(second))))
                        if pair in seen:continue
                        seen.add(pair);faces2=np.array(pair);old=labels[faces2].copy();touched=np.unique(f[faces2]);before_local=int(scores[touched].sum());edges2={tuple(sorted((int(i),int(j))))for i in faces2 for j in neighbor[i]};before_edges=sum(labels[i]!=labels[j]for i,j in edges2)
                        options=[np.unique(labels[np.r_[neighbor[i],i]])for i in faces2]
                        for c1 in options[0]:
                            for c2 in options[1]:
                                new=np.array([c1,c2]);delta_counts=np.bincount(new,minlength=len(counts))-np.bincount(old,minlength=len(counts))
                                if np.any(counts+delta_counts<=0):continue
                                labels[faces2]=new;after_local=np.array([score(x)for x in touched]);delta=int(after_local.sum())-before_local;edge_delta=sum(labels[i]!=labels[j]for i,j in edges2)-before_edges;labels[faces2]=old
                                if delta<0 or delta==0 and edge_delta<0:
                                    key=(delta,edge_delta,*pair,int(c1),int(c2))
                                    if best is None or key<best[0]:best=(key,faces2,new,delta_counts,touched,after_local)
            if best:
                _,faces2,new,dc,touched,after_local=best;labels[faces2]=new;counts+=dc;scores[touched]=after_local;moves+=2
            else:
                # Try a complete existing vertex star at a multi-label
                # junction. This remains an index selection of the carrier.
                best=None
                for v in np.flatnonzero(scores):
                    for faces2 in [fans[v],np.unique(np.r_[fans[v],neighbor[fans[v]].ravel()])]:
                        old=labels[faces2].copy();touched=np.unique(f[faces2]);before_local=int(scores[touched].sum());edges2={tuple(sorted((int(i),int(j))))for i in faces2 for j in neighbor[i]};before_edges=sum(labels[i]!=labels[j]for i,j in edges2)
                        for color in np.unique(labels[np.r_[faces2,neighbor[faces2].ravel()]]):
                            new=np.full(len(faces2),color);dc=np.bincount(new,minlength=len(counts))-np.bincount(old,minlength=len(counts))
                            if np.any(counts+dc<=0):continue
                            labels[faces2]=new;after_local=np.array([score(x)for x in touched]);delta=int(after_local.sum())-before_local;edge_delta=sum(labels[i]!=labels[j]for i,j in edges2)-before_edges;labels[faces2]=old
                            if delta<0 or delta==0 and edge_delta<0:
                                key=(delta,edge_delta,len(faces2),int(color))
                                if best is None or key<best[0]:best=(key,faces2,new,dc,touched,after_local)
                if best:
                    _,faces2,new,dc,touched,after_local=best;labels[faces2]=new;counts+=dc;scores[touched]=after_local;moves+=len(faces2)
                else:break
    result=dict(source_triangles=n,initial_pinched_fans=before,remaining_pinched_fans=int(scores.sum()),changed_triangle_labels=int(np.count_nonzero(labels!=original)),moves=moves,passes=iteration+1,geometry_unchanged=True,all_regions_retained=bool(np.all(counts>0)))
    if result['remaining_pinched_fans']:raise ValueError('Selection boundary regularization did not converge: '+str(result))
    return labels,result
