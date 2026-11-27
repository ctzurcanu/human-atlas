"""Exact affine-weight identities for overlapping source-triangle charts.

This overlays caller-supplied triangles; it invents no tissue or boundary.
Fractions retain the literal binary input values without integer quantization.
Every shared point is evaluated once from immutable original controls.
"""
from fractions import Fraction
import numpy as np


def weight_key(weights):
    return tuple((int(i),w.numerator,w.denominator)for i,w in sorted(weights.items())if w)


def weights_from_key(key):
    return {int(i):Fraction(int(n),int(d))for i,n,d in key}


def source_weights(point,source_ids,positions):
    """Project onto the explicitly declared simplex, using exact arithmetic.

    The caller must check the actual 3D replay error. This is not a nearest
    simplex search or permission to project unrelated geometry.
    """
    ids=tuple(sorted(map(int,source_ids)));p=np.asarray(positions,float);q=np.asarray(point,float)
    if not 1<=len(ids)<=3:raise ValueError('One declared source vertex, edge or triangle required')
    if len(ids)==1:return {ids[0]:Fraction(1)}
    a,b=p[list(ids[:2])];axis=int(np.argmax(abs(b-a)))
    def number(x):return Fraction.from_float(float(x))
    if len(ids)==2:
        t=(number(q[axis])-number(a[axis]))/(number(b[axis])-number(a[axis]))
        return {ids[0]:1-t,ids[1]:t}
    c=p[ids[2]];normal=np.cross(b-a,c-a);drop=int(np.argmax(abs(normal)));axes=[i for i in range(3)if i!=drop]
    u,v=axes;aa=number(b[u])-number(a[u]);ab=number(c[u])-number(a[u]);ba=number(b[v])-number(a[v]);bb=number(c[v])-number(a[v]);det=aa*bb-ab*ba
    if not det:raise ValueError('Declared source triangle is degenerate')
    x=number(q[u])-number(a[u]);y=number(q[v])-number(a[v]);s=(x*bb-ab*y)/det;t=(aa*y-x*ba)/det
    return {ids[0]:1-s-t,ids[1]:s,ids[2]:t}


class RationalSurfaceBank:
    def __init__(self,positions):
        self.source=np.asarray(positions,float).copy()
        if self.source.ndim!=2 or self.source.shape[1]!=3 or not np.isfinite(self.source).all():raise ValueError('Finite source controls required')
        self.rows=self.source.tolist();self.weights={};self.lookup={};self.aliases={};self.maximum_replay_error_m=0.

    def evaluate(self,weights):
        ids=sorted(weights);origin=self.source[ids[0]]
        return origin+sum((float(weights[i])*(self.source[i]-origin)for i in ids[1:]),np.zeros(3))

    @staticmethod
    def validate(weights):
        if not weights or any(not isinstance(w,Fraction)for w in weights.values())or sum(weights.values(),Fraction(0))!=1 or min(weights.values())<0:
            raise ValueError('Nonnegative exact source weights summing to one required')

    def register(self,index,weights,*,maximum_replay_error_m):
        weights={int(i):w for i,w in weights.items()if w};self.validate(weights)
        if not 0<=index<len(self.source)or not np.isfinite(maximum_replay_error_m)or maximum_replay_error_m<0:raise ValueError('Explicit valid source index and arithmetic bound required')
        q=self.evaluate(weights);error=float(np.linalg.norm(q-self.source[index]));self.maximum_replay_error_m=max(self.maximum_replay_error_m,error)
        if error>maximum_replay_error_m:raise ValueError('Declared simplex replay exceeds its arithmetic bound')
        key=weight_key(weights)
        if key in self.lookup:
            alias=self.lookup[key];self.aliases[index]=alias;return alias
        self.lookup[key]=index;self.weights[index]=weights;self.rows[index]=q.tolist();return index

    def add(self,weights):
        weights={int(i):w for i,w in weights.items()if w};self.validate(weights);key=weight_key(weights)
        if key not in self.lookup:
            index=len(self.rows);self.lookup[key]=index;self.weights[index]=weights;self.rows.append(self.evaluate(weights).tolist())
        return self.lookup[key]

    def mix(self,a,b,t):
        if t==0:return a
        if t==1:return b
        if not 0<t<1:raise ValueError('Intersection leaves its actual segment')
        wa,wb=self.weights[a],self.weights[b];ids=set(wa)|set(wb)
        return self.add({i:(1-t)*wa.get(i,Fraction(0))+t*wb.get(i,Fraction(0))for i in ids})

    def xy(self,index,source_face):
        a,b,c=map(int,source_face);w=self.weights[index]
        if set(w)-{a,b,c}:raise ValueError('Point is not in the declared source facet')
        return w.get(b,Fraction(0)),w.get(c,Fraction(0))

    def signed_area2(self,polygon,source_face):
        xy=[self.xy(i,source_face)for i in polygon]
        return sum((a[0]*b[1]-a[1]*b[0]for a,b in zip(xy,xy[1:]+xy[:1])),Fraction(0))

    def triangle_intersection(self,first,second,source_face):
        """Return one convex polygon, positive in source barycentric axes."""
        polygon=list(map(int,first));clip=list(map(int,second))
        if self.signed_area2(polygon,source_face)<0:polygon.reverse()
        if self.signed_area2(clip,source_face)<0:clip.reverse()
        if not self.signed_area2(polygon,source_face)or not self.signed_area2(clip,source_face):raise ValueError('Zero-area input chart')
        for a,b in zip(clip,clip[1:]+clip[:1]):
            x,y=self.xy(a,source_face),self.xy(b,source_face)
            def value(i):
                q=self.xy(i,source_face);return (y[0]-x[0])*(q[1]-x[1])-(y[1]-x[1])*(q[0]-x[0])
            result=[]
            for c,d in zip(polygon,polygon[1:]+polygon[:1]):
                dc,dd=value(c),value(d)
                if dc>=0:result.append(c)
                if dc*dd<0:result.append(self.mix(c,d,dc/(dc-dd)))
            polygon=[]
            for i in result:
                if not polygon or polygon[-1]!=i:polygon.append(i)
            if len(polygon)>1 and polygon[0]==polygon[-1]:polygon.pop()
            if len(polygon)<3:return []
        area=self.signed_area2(polygon,source_face)
        if area==0:return []
        if area<0:raise ValueError('Reversed intersection polygon')
        return polygon

    def triangulate(self,polygon,source_face):
        if len(polygon)==3:return [list(polygon)]
        weights={}
        for index in polygon:
            for i,w in self.weights[index].items():weights[i]=weights.get(i,Fraction(0))+w/len(polygon)
        center=self.add(weights);faces=[[center,a,b]for a,b in zip(polygon,polygon[1:]+polygon[:1])]
        if any(self.signed_area2(f,source_face)<=0 for f in faces):raise ValueError('Collapsed overlay triangulation')
        return faces

    def segment_fraction(self,index,first,second):
        """Exact affine collinearity, not a distance-based edge weld."""
        q,a,b=self.weights[index],self.weights[first],self.weights[second];ids=set(q)|set(a)|set(b);t=None
        for i in ids:
            delta=b.get(i,Fraction(0))-a.get(i,Fraction(0))
            if delta:t=(q.get(i,Fraction(0))-a.get(i,Fraction(0)))/delta;break
        if t is None:return None
        if all(q.get(i,Fraction(0))==a.get(i,Fraction(0))+t*(b.get(i,Fraction(0))-a.get(i,Fraction(0)))for i in ids):return t
        return None
