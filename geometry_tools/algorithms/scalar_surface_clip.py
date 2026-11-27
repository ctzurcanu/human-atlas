"""Split an existing surface by an explicit piecewise-linear scalar field.

Original controls remain a literal prefix. A new control belongs to one
original edge and is reused by both sides and every incident face. Neither
side receives a cap: an anatomical opening is not a separate solid.
"""
import numpy as np


def clip_scalar_surface(positions,faces,values,level):
    v=np.asarray(positions,float);f=np.asarray(faces);s=np.asarray(values,float)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all():raise ValueError('Finite positions required')
    if f.ndim!=2 or f.shape[1]!=3 or not np.issubdtype(f.dtype,np.integer)or not len(f)or f.min()<0 or f.max()>=len(v):
        raise ValueError('Valid nonempty integer source triangles required')
    if s.shape!=(len(v),)or not np.isfinite(s).all()or not np.isfinite(level):raise ValueError('One finite scalar per source control required')
    d=s-float(level)
    if np.any(np.all(d[f]==0,axis=1)):raise ValueError('A whole source triangle lies on the isovalue; its ownership is unresolved')
    rows=v.tolist();cache={};edge_parents=[];fractions=[];positive=[];negative=[];positive_parents=[];negative_parents=[]
    def crossing(a,b):
        if d[a]==0:return int(a)
        if d[b]==0:return int(b)
        a,b=sorted((int(a),int(b)));key=(a,b)
        if key not in cache:
            mu=float(d[a]/(d[a]-d[b]))
            if not 0<mu<1:
                raise ValueError('An edge crossing is not representable in Float64; declare a plane-classification error bound before clipping')
            cache[key]=len(rows);rows.append((v[a]+mu*(v[b]-v[a])).tolist());edge_parents.append([a,b]);fractions.append(mu)
        return cache[key]
    def polygon(t,side):
        out=[]
        for a,b in zip(t,np.roll(t,-1)):
            if side*d[a]>=0:out.append(int(a))
            if d[a]*d[b]<0:out.append(crossing(a,b))
        return out
    def triangles(poly):
        if len(poly)<3:return []
        if len(poly)==3:return [poly]
        assert len(poly)==4
        a,b,c,d=poly;options=[[[a,b,c],[a,c,d]],[[a,b,d],[b,c,d]]]
        def quality(fs):
            # Only these six corners participate in this diagonal test.
            # Rebuilding the complete growing bank per crossing facet made
            # otherwise linear section queries quadratic in source size.
            p=np.asarray([[rows[int(i)]for i in face]for face in fs]);n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);length=np.linalg.norm(np.roll(p,-1,axis=1)-p,axis=2).max(1)
            return float(np.min(np.linalg.norm(n,axis=1)/length**2))
        return max(options,key=quality)
    for parent,t in enumerate(f):
        for side,out,owners in [(1,positive,positive_parents),(-1,negative,negative_parents)]:
            child=triangles(polygon(t,side));out.extend(child);owners.extend([parent]*len(child))
    return dict(positions=np.asarray(rows,float),positive_faces=np.asarray(positive,int).reshape(-1,3),
        negative_faces=np.asarray(negative,int).reshape(-1,3),positive_parents=np.asarray(positive_parents,int),
        negative_parents=np.asarray(negative_parents,int),new_edge_parents=np.asarray(edge_parents,int).reshape(-1,2),
        fractions=np.asarray(fractions,float),original_controls_literal=True,caps_added=0)


def interpolate_clip_controls(source,new_edge_parents,fractions):
    """Carry the exact logical edge subdivision to corresponding interfaces."""
    v=np.asarray(source,float);edges=np.asarray(new_edge_parents,int);mu=np.asarray(fractions,float)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all()or edges.shape!=(len(mu),2):raise ValueError('One finite corresponding position bank required')
    if len(edges)and(edges.min()<0 or edges.max()>=len(v)or np.any(mu<=0)or np.any(mu>=1)):raise ValueError('Invalid original edge interpolation')
    return np.r_[v,v[edges[:,0]]+mu[:,None]*(v[edges[:,1]]-v[edges[:,0]])]
