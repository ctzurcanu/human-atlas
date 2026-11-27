"""Batched closest points between actual triangles, including edge/face hits.

Float64 arithmetic in local coordinates. No nearest-vertex approximation or
clinical tolerance is used to turn a gap into a join. Degenerate inputs retain
their segment/point distance and are separately identified.
"""
import numpy as np
from fractions import Fraction


def exact_intersection_segment(first,second):
    """Construct the noncoplanar intersection of two binary-valued triangles.

    Fraction arithmetic establishes the endpoints before the single Float64
    output conversion. The reported coordinate error bounds that conversion;
    it is not an anatomical tolerance or an allowance for a geometric gap.
    Coplanar triangles require a separate polygon operation and are explicit.
    """
    inputs=[np.asarray(t,float)for t in (first,second)]
    if any(t.shape!=(3,3)or not np.isfinite(t).all()for t in inputs):
        raise ValueError('Two finite 3 by 3 triangles required')
    a,b=[[[Fraction(float(x))for x in p]for p in t]for t in inputs]
    sub=lambda p,q:[x-y for x,y in zip(p,q)]
    dot=lambda p,q:sum(x*y for x,y in zip(p,q))
    cross=lambda p,q:[p[1]*q[2]-p[2]*q[1],p[2]*q[0]-p[0]*q[2],p[0]*q[1]-p[1]*q[0]]
    na=cross(sub(a[1],a[0]),sub(a[2],a[0]));nb=cross(sub(b[1],b[0]),sub(b[2],b[0]))
    empty=lambda relation:dict(relation=relation,endpoints=np.empty((0,3)),coordinate_error_upper_m=np.empty((0,3)))
    if not any(na)or not any(nb):return empty('degenerate')
    da=[dot(sub(p,b[0]),nb)for p in a];db=[dot(sub(p,a[0]),na)for p in b]
    line=cross(na,nb)
    if not any(line):return empty('coplanar'if not any(da)else'disjoint')
    axis=max(range(3),key=lambda k:abs(line[k]))
    def section(t,d):
        points=[t[i]for i in range(3)if d[i]==0]
        for i,j in [(0,1),(1,2),(2,0)]:
            if d[i]*d[j]<0:
                mu=d[i]/(d[i]-d[j]);points.append([x+mu*(y-x)for x,y in zip(t[i],t[j])])
        points=list(dict.fromkeys(tuple(p)for p in points))
        return sorted(points,key=lambda p:p[axis])
    sa,sb=section(a,da),section(b,db)
    if not sa or not sb:return empty('disjoint')
    lo=max(sa[0],sb[0],key=lambda p:p[axis]);hi=min(sa[-1],sb[-1],key=lambda p:p[axis])
    if lo[axis]>hi[axis]:return empty('disjoint')
    exact=[lo,hi]
    # Each exact endpoint lies on both exact source planes. This assertion
    # concerns the rational construction, before global-coordinate rounding.
    assert all(dot(sub(p,a[0]),na)==0 and dot(sub(p,b[0]),nb)==0 for p in exact)
    points=np.asarray([[float(x)for x in p]for p in exact])
    error=np.asarray([[np.nextafter(float(abs(Fraction(float(q))-x)),np.inf)if Fraction(float(q))!=x else 0.
                       for q,x in zip(rounded,p)]for rounded,p in zip(points,exact)])
    proper=lo[axis]<hi[axis]and min(da)<0<max(da)and min(db)<0<max(db)
    return dict(relation='crossing'if proper else'touching',endpoints=points,coordinate_error_upper_m=error)


def exact_transverse_crossing(first,second):
    """Exact rational predicate for a suspected noncoplanar interior crossing.

    Float inputs are interpreted as their exact binary values. Both planes
    must cut triangle interiors, and their exact intersection intervals must
    overlap with positive length. Shared edges alone cannot pass this test.
    This slow path is used only to confirm suspected defects, not all pairs.
    """
    a=[[Fraction(float(x))for x in p]for p in first];b=[[Fraction(float(x))for x in p]for p in second]
    sub=lambda p,q:[x-y for x,y in zip(p,q)]
    dot=lambda p,q:sum(x*y for x,y in zip(p,q))
    cross=lambda p,q:[p[1]*q[2]-p[2]*q[1],p[2]*q[0]-p[0]*q[2],p[0]*q[1]-p[1]*q[0]]
    na=cross(sub(a[1],a[0]),sub(a[2],a[0]));nb=cross(sub(b[1],b[0]),sub(b[2],b[0]));line=cross(na,nb)
    if not any(line):return False  # Coplanar area is checked separately.
    da=[dot(sub(p,b[0]),nb)for p in a];db=[dot(sub(p,a[0]),na)for p in b]
    if not(min(da)<0<max(da)and min(db)<0<max(db)):return False
    axis=max(range(3),key=lambda k:abs(line[k]))
    def interval(t,d):
        values=[t[i][axis]for i in range(3)if d[i]==0]
        for i,j in [(0,1),(1,2),(2,0)]:
            if d[i]*d[j]<0:values.append(t[i][axis]+(t[j][axis]-t[i][axis])*d[i]/(d[i]-d[j]))
        return min(values),max(values)
    ia,ib=interval(a,da),interval(b,db)
    return max(ia[0],ib[0])<min(ia[1],ib[1])


def barycentric(points, triangles):
    p,t=np.asarray(points,float),np.asarray(triangles,float)
    e0,e1=t[:,1]-t[:,0],t[:,2]-t[:,0];n=np.cross(e0,e1);d=np.einsum('ij,ij->i',n,n)
    a=np.divide(np.einsum('ij,ij->i',np.cross(p-t[:,0],e1),n),d,out=np.zeros(len(t)),where=d>0)
    b=np.divide(np.einsum('ij,ij->i',np.cross(e0,p-t[:,0]),n),d,out=np.zeros(len(t)),where=d>0)
    return np.c_[1-a-b,a,b],d>0


def point_triangle(points, triangles):
    p,t=np.asarray(points,float),np.asarray(triangles,float);n=len(p)
    best=np.full(n,np.inf);out=np.zeros_like(p);weights=np.zeros((n,3))
    for i,j in [(0,1),(1,2),(2,0)]:
        e=t[:,j]-t[:,i];den=np.einsum('ij,ij->i',e,e)
        mu=np.clip(np.divide(np.einsum('ij,ij->i',p-t[:,i],e),den,out=np.zeros(n),where=den>0),0,1)
        q=t[:,i]+mu[:,None]*e;d=np.einsum('ij,ij->i',q-p,q-p);take=d<best
        out[take]=q[take];best[take]=d[take];weights[take]=0;weights[take,i]=1-mu[take];weights[take,j]=mu[take]
    w,valid=barycentric(p,t);q=np.einsum('ni,nij->nj',w,t);d=np.einsum('ij,ij->i',q-p,q-p)
    take=valid&np.all(w>=0,axis=1)&(d<best);out[take]=q[take];weights[take]=w[take];best[take]=d[take]
    return best,out,weights


def _coplanar_polygon(a,b,normal):
    axis=int(np.argmax(abs(normal)));axes=[i for i in range(3)if i!=axis]
    origin=a[0];aa=(a-origin)[:,axes];bb=(b-origin)[:,axes];poly=aa.tolist()
    cross=lambda x,y:x[0]*y[1]-x[1]*y[0]
    sign=np.sign(cross(bb[1]-bb[0],bb[2]-bb[0]))
    if sign==0:return 0.,None
    for x,y in zip(bb,np.roll(bb,-1,axis=0)):
        result=[]
        if not poly:return 0.,None
        for u,v in zip(poly,poly[1:]+poly[:1]):
            du=sign*cross(y-x,np.asarray(u)-x);dv=sign*cross(y-x,np.asarray(v)-x)
            if du>=0:result.append(u)
            if(du>=0)!=(dv>=0):
                mu=du/(du-dv);result.append((np.asarray(u)+mu*(np.asarray(v)-u)).tolist())
        poly=result
    if len(poly)<3:return 0.,None
    q=np.asarray(poly);area=abs(np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1)))/2
    area/=max(abs(normal[axis]),np.finfo(float).tiny)
    mid2=q.mean(0);mid=np.zeros(3);mid[axes]=mid2;mid[axis]=-np.dot(normal[axes],mid2)/normal[axis]
    return float(area),mid+origin


def triangle_pairs(first,second):
    """One witness per pair. Relations: near, touching, crossing, coplanar-overlap.

    Euclidean distance is minimized over six vertex/face and nine interior
    edge/edge candidates, with edge/face intersections included explicitly.
    """
    a,b=np.asarray(first,float),np.asarray(second,float)
    if a.shape!=b.shape or a.ndim!=3 or a.shape[1:]!=(3,3):raise ValueError('Matching triangle batches required')
    count=len(a);origin=a[:,0].copy();a=a-origin[:,None];b=b-origin[:,None]
    best=np.full(count,np.inf);pa=np.zeros((count,3));pb=pa.copy();wa=pa.copy();wb=pa.copy()
    def keep(d,x,y,u,v):
        take=d<best;best[take]=d[take];pa[take]=x[take];pb[take]=y[take];wa[take]=u[take];wb[take]=v[take]
    eye=np.eye(3)
    for i in range(3):
        d,q,w=point_triangle(a[:,i],b);keep(d,a[:,i],q,np.broadcast_to(eye[i],(count,3)),w)
        d,q,w=point_triangle(b[:,i],a);keep(d,q,b[:,i],w,np.broadcast_to(eye[i],(count,3)))
    for i,j in [(0,1),(1,2),(2,0)]:
      u=a[:,j]-a[:,i];A=np.einsum('ij,ij->i',u,u)
      for k,l in [(0,1),(1,2),(2,0)]:
        v=b[:,l]-b[:,k];r=a[:,i]-b[:,k];cross=np.cross(u,v);den=np.einsum('ij,ij->i',cross,cross)
        s=np.divide(-np.einsum('ij,ij->i',np.cross(r,v),cross),den,out=np.zeros(count),where=den>0);t=np.divide(-np.einsum('ij,ij->i',np.cross(r,u),cross),den,out=np.zeros(count),where=den>0)
        valid=(den>0)&(s>=0)&(s<=1)&(t>=0)&(t<=1)
        x=a[:,i]+s[:,None]*u;y=b[:,k]+t[:,None]*v;d=np.einsum('ij,ij->i',x-y,x-y);d[~valid]=np.inf
        w1=np.zeros((count,3));w2=w1.copy();w1[:,i]=1-s;w1[:,j]=s;w2[:,k]=1-t;w2[:,l]=t;keep(d,x,y,w1,w2)
    na=np.cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);nb=np.cross(b[:,1]-b[:,0],b[:,2]-b[:,0]);la=np.linalg.norm(na,axis=1);lb=np.linalg.norm(nb,axis=1)
    valid=(la>0)&(lb>0);un=np.divide(na,la[:,None],out=np.zeros_like(na),where=la[:,None]>0);vn=np.divide(nb,lb[:,None],out=np.zeros_like(nb),where=lb[:,None]>0)
    hit_sum=np.zeros((count,3));hit_count=np.zeros(count,int)
    for source,target,n in [(a,b,nb),(b,a,na)]:
      for i,j in [(0,1),(1,2),(2,0)]:
        d0=np.einsum('ij,ij->i',source[:,i]-target[:,0],n);d1=np.einsum('ij,ij->i',source[:,j]-target[:,0],n);den=d0-d1
        mu=np.divide(d0,den,out=np.zeros(count),where=den!=0);point=source[:,i]+mu[:,None]*(source[:,j]-source[:,i]);w,good=barycentric(point,target)
        take=good&(np.linalg.norm(n,axis=1)>0)&(den!=0)&(mu>=0)&(mu<=1)&np.all(w>=0,axis=1)
        hit_sum[take]+=point[take];hit_count[take]+=1
    intersect=hit_count>0
    mid=np.divide(hit_sum,hit_count[:,None],out=np.zeros_like(hit_sum),where=hit_count[:,None]>0)
    previous={k:v[intersect].copy()for k,v in [('best',best),('pa',pa),('pb',pb),('wa',wa),('wb',wb)]}
    best[intersect]=0;pa[intersect]=pb[intersect]=mid[intersect];u,_=barycentric(mid,a);v,_=barycentric(mid,b);wa[intersect]=u[intersect];wb[intersect]=v[intersect]
    scale=np.maximum(np.max(abs(a),axis=(1,2)),np.max(abs(b),axis=(1,2)))
    roundoff=128*np.finfo(float).eps*np.maximum(scale,np.finfo(float).tiny)
    relation=np.full(count,'near',dtype='<U32');touch=np.sqrt(best)<=roundoff;relation[touch]='touching'
    parallel=np.linalg.norm(np.cross(un,vn),axis=1)<=128*np.finfo(float).eps
    proper=intersect&~parallel&np.all(wa>128*np.finfo(float).eps,axis=1)&np.all(wb>128*np.finfo(float).eps,axis=1)
    rejected=np.array([i for i in np.flatnonzero(proper)if not exact_transverse_crossing(a[i],b[i])],int)
    if len(rejected):
        at=np.searchsorted(np.flatnonzero(intersect),rejected)
        for name,array in [('best',best),('pa',pa),('pb',pb),('wa',wa),('wb',wb)]:array[rejected]=previous[name][at]
        proper[rejected]=False;touch[rejected]=np.sqrt(best[rejected])<=roundoff[rejected];relation[rejected]=np.where(touch[rejected],'touching','near')
    relation[proper]='crossing';area=np.zeros(count)
    coplanar=touch&parallel&valid&(np.max(abs(np.einsum('nki,ni->nk',b-a[:,0,None],un)),axis=1)<=roundoff)
    for i in np.flatnonzero(coplanar):
        ar,point=_coplanar_polygon(a[i],b[i],un[i]);bound=512*np.finfo(float).eps*scale[i]*scale[i]
        if ar>bound:
            area[i]=ar;relation[i]='coplanar-overlap';pa[i]=pb[i]=point;wa[i]=barycentric(point[None],a[i:i+1])[0][0];wb[i]=barycentric(point[None],b[i:i+1])[0][0];best[i]=0
    return dict(distance_m=np.sqrt(np.maximum(best,0)),first=pa+origin,second=pb+origin,
                first_barycentric=wa,second_barycentric=wb,relation=relation,coplanar_overlap_area_m2=area,
                normal_cosine=np.einsum('ij,ij->i',un,vn),degenerate=~valid,numerical_roundoff_m=roundoff)
