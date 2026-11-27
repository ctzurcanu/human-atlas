"""Propagate exact geometric edge subdivisions through reciprocal surfaces.

This never bridges distinct curves: both boundaries must replay the same
oriented polygon after inserting collinear controls within the stated bound.
"""
from collections import Counter,defaultdict
import numpy as np

def boundary(faces):
    edges=Counter((int(t[k]),int(t[(k+1)%3]))for t in faces for k in range(3))
    return Counter({e:n-edges[e[::-1]]for e,n in edges.items()if n>edges[e[::-1]]})

def reconcile_boundaries(positions,actual,expected,plans,tolerance=2e-12):
    v=np.asarray(positions,float);a=boundary(actual);b=boundary(expected);maximum=0.
    def insert_new(source,target):
        nonlocal maximum
        edges=np.array(list(target),int);start=v[edges[:,0]];delta=v[edges[:,1]]-start;length2=np.einsum('ij,ij->i',delta,delta)
        new=set(np.asarray(list(source),int).ravel())-set(edges.ravel())
        for control in new:
            t=np.clip(np.einsum('ij,ij->i',v[control]-start,delta)/length2,0,1);distance=np.linalg.norm(v[control]-start-t[:,None]*delta,axis=1);i=int(np.argmin(distance));maximum=max(maximum,float(distance[i]))
            if distance[i]>tolerance or not 0<t[i]<1:
                raise ValueError('Boundary differs geometrically; collinear refinement is not a repair; '
                    f'control={control}, closest_edge={edges[i].tolist()}, '
                    f'distance_m={distance[i]:.17g}, fraction={t[i]:.17g}, '
                    f'control_m={v[control].tolist()}, edge_m={v[edges[i]].tolist()}')
            plans[tuple(sorted(map(int,edges[i])))].add(int(control))
    insert_new(a,b);insert_new(b,a)
    def expanded(source):
        result=Counter()
        for (x,y),count in source.items():
            points=sorted(plans.get(tuple(sorted((x,y))),()),key=lambda i:float((v[i]-v[x])@(v[y]-v[x])))
            chain=[x]+points+[y]
            for i,j in zip(chain,chain[1:]):result[(i,j)]+=count
        return result
    if expanded(a)!=expanded(b):raise ValueError('Refined directed boundaries still differ')
    return maximum

def refine_surfaces(bank,surfaces,plans):
    output={};ancestry={};changed=0
    for sid,faces in surfaces.items():
        triangles=[];parents=[]
        for ordinal,face in enumerate(faces):
            polygon=[];edge_membership=[]
            for edge,(x,y) in enumerate(zip(face,np.roll(face,-1))):
                x,y=int(x),int(y);polygon.append(x);edge_membership.append({edge,(edge-1)%3})
                points=plans.get(tuple(sorted((x,y))),())
                points=sorted(points,key=lambda i:float((np.array(bank.rows[i])-bank.rows[x])@(np.array(bank.rows[y])-bank.rows[x])))
                polygon.extend(points);edge_membership.extend([{edge}for _ in points])
            if len(polygon)==3:triangles.append(face);parents.append(ordinal);continue
            # No new centroid: for a very thin parent it could be snapped to a
            # boundary control. Clip convex ears while retaining every edge
            # subdivision. Edge membership detects a collinear remainder by
            # construction, rather than by a geometric welding tolerance.
            original=np.array([bank.rows[int(i)]for i in face]);normal=np.cross(original[1]-original[0],original[2]-original[0])
            if not np.linalg.norm(normal)>0:raise ValueError('Cannot refine a zero-area source triangle')
            normal/=np.linalg.norm(normal);active=list(range(len(polygon)));changed+=1
            while len(active)>3:
                candidates=[]
                for slot,middle in enumerate(active):
                    previous,following=active[slot-1],active[(slot+1)%len(active)]
                    remaining=active[:slot]+active[slot+1:]
                    if set.intersection(*(edge_membership[i]for i in remaining)):continue
                    a,b,c=(np.asarray(bank.rows[polygon[i]])for i in(previous,middle,following))
                    area=float(np.cross(b-a,c-a)@normal)
                    if area>0:candidates.append((area,slot,(polygon[previous],polygon[middle],polygon[following])))
                if not candidates:raise ValueError('Subdivision polygon has no positive noncollapsing ear')
                _,slot,triangle=max(candidates,key=lambda row:row[0]);triangles.append(triangle);parents.append(ordinal);active.pop(slot)
            triangle=tuple(polygon[i]for i in active);xyz=np.asarray([bank.rows[i]for i in triangle])
            if not np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0])@normal>0:raise ValueError('Subdivision has a collapsed final triangle')
            triangles.append(triangle);parents.append(ordinal)
        output[sid]=np.asarray(triangles,int);ancestry[sid]=np.asarray(parents,int)
    return output,ancestry,changed

def bisect_surfaces(surfaces,midpoints):
    """Conforming midpoint refinement without adding face-centroid controls.

    Every geometric edge is split in all incident surfaces. Surface orientation
    and parent-face ancestry are retained for subsequent nonlinear embedding.
    """
    output={};ancestry={};changed=0
    for sid,faces in surfaces.items():
        triangles=[];parents=[]
        for ordinal,face in enumerate(faces):
            marked=[tuple(sorted((int(face[i]),int(face[(i+1)%3]))))in midpoints for i in range(3)];count=sum(marked)
            if not count:triangles.append(face);parents.append(ordinal);continue
            changed+=1
            shift=marked.index(True)if count==1 else(marked.index(False)+1)%3 if count==2 else 0
            a,b,c=map(int,np.roll(face,-shift));ab=midpoints.get(tuple(sorted((a,b))));bc=midpoints.get(tuple(sorted((b,c))));ca=midpoints.get(tuple(sorted((c,a))))
            new=[(a,ab,c),(ab,b,c)]if count==1 else[(b,bc,ab),(a,ab,c),(ab,bc,c)]if count==2 else[(a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)]
            triangles.extend(new);parents.extend([ordinal]*len(new))
        output[sid]=np.asarray(triangles,int);ancestry[sid]=np.asarray(parents,int)
    return output,ancestry,changed
