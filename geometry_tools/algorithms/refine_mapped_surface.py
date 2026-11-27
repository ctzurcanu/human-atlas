"""Conforming local refinement before an explicit nonlinear coordinate map.

Original controls stay literal. New controls lie on the original linear
facets, then pass through the supplied map. Edge controls are shared and
protected open rims cannot be subdivided. This is tessellation, not a repair
certificate: callers must audit the actual mapped triangles and interfaces.
"""
from collections import Counter
import numpy as np


def refine_mapped_surface(source,mapped,faces,marked_edges,map_points,*,protect_open_rims=True,centroid_fan=True):
    source=np.asarray(source,float);mapped=np.asarray(mapped,float);faces=np.asarray(faces,np.int64)
    if source.shape!=mapped.shape or source.ndim!=2 or source.shape[1]!=3:
        raise ValueError('Matching source/mapped 3D controls are required')
    counts=Counter(tuple(sorted([int(a),int(b)]))for t in faces for a,b in zip(t,np.roll(t,-1)))
    marked={tuple(sorted(map(int,e)))for e in marked_edges}
    if any(e not in counts for e in marked):raise ValueError('Every marked edge must belong to the actual surface')
    if protect_open_rims and any(counts[e]==1 for e in marked):raise ValueError('A protected open rim cannot be subdivided')
    extra=[];edge_ids={}
    for a,b in sorted(marked):edge_ids[(a,b)]=len(source)+len(extra);extra.append((source[a]+source[b])*.5)
    result=[];parents=[]
    for i,tri in enumerate(faces):
        if not centroid_fan:
            a,b,c=map(int,tri);ab=edge_ids.get(tuple(sorted((a,b))));bc=edge_ids.get(tuple(sorted((b,c))));ca=edge_ids.get(tuple(sorted((c,a))))
            count=sum(q is not None for q in [ab,bc,ca])
            if count==0:rows=[[a,b,c]]
            elif count==3:rows=[[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]
            else:
                # Rotate a marked edge to AB. With two marks arrange the
                # second at BC. This produces two/three source subfacets,
                # rather than adding an unnecessary centroid fan.
                while ab is None or count==2 and bc is None:
                    a,b,c=b,c,a;ab,bc,ca=bc,ca,ab
                rows=[[a,ab,c],[ab,b,c]]if count==1 else[[a,ab,c],[ab,bc,c],[ab,b,bc]]
            result.extend(rows);parents.extend([i]*len(rows));continue
        boundary=[];changed=False
        for a,b in zip(tri,np.roll(tri,-1)):
            boundary.append(int(a));mid=edge_ids.get(tuple(sorted([int(a),int(b)])))
            if mid is not None:boundary.append(mid);changed=True
        if not changed:result.append(tri.tolist());parents.append(i);continue
        center=len(source)+len(extra);extra.append(source[tri].mean(0))
        for a,b in zip(boundary,np.roll(boundary,-1)):result.append([center,int(a),int(b)]);parents.append(i)
    new_source=np.asarray(extra,float).reshape(-1,3);new_mapped=np.asarray(map_points(new_source),float)
    if new_mapped.shape!=new_source.shape or not np.isfinite(new_mapped).all():raise ValueError('The map returned invalid new controls')
    return dict(source_positions=np.r_[source,new_source],positions=np.r_[mapped,new_mapped],
                faces=np.asarray(result,np.int64),parent_faces=np.asarray(parents,np.int64),
                new_controls=len(extra),shared_edge_controls=len(edge_ids),original_controls_literal=True,
                source_facets_preserved=True,anatomical_acceptance=False)
