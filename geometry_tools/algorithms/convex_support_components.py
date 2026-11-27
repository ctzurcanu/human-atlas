"""Source-indexed convex constraint domains for disjoint concave outlines.

These domains are mathematical support tests, never new anatomy or caps.
No new boundary coordinates are introduced by the constrained triangulation.
"""
import time
import numpy as np
import shapely


def separating_components(first, second, *, clearance_m, time_limit_seconds=2.):
    start=time.monotonic();a=np.asarray(first,float);b=np.asarray(second,float)
    if (a.ndim!=2 or b.ndim!=2 or a.shape[1]!=2 or b.shape[1]!=2
            or not np.isfinite(a).all() or not np.isfinite(b).all()
            or min(len(a),len(b))<2 or clearance_m<=0 or not 0<time_limit_seconds<=5):
        raise ValueError('Finite source-indexed outlines and positive clearance required')

    def shape(p):return shapely.linestrings(p)if len(p)==2 else shapely.polygons(p)

    def support(x,y):
        line=shapely.get_coordinates(shapely.shortest_line(shapely.convex_hull(shape(y)),shapely.convex_hull(shape(x))))
        if len(line)!=2:return None
        n=line[1]-line[0];length=np.linalg.norm(n)
        if length<1e-14:return None
        n/=length;amin=float((x@n).min());bmax=float((y@n).max())
        return (n,amin,bmax)if amin-bmax>=clearance_m-5e-12 else None

    whole=support(a,b)
    if whole is not None:return [(np.arange(len(a)),np.arange(len(b)),*whole)]

    def parts(p):
        s=shape(p)
        if not shapely.is_valid(s):raise ValueError('Invalid source outline; never repair by triangulation')
        if len(p)==2 or shapely.equals(s,shapely.convex_hull(s)):return [np.arange(len(p))]
        index={tuple(q):i for i,q in enumerate(p)};triangles=shapely.get_parts(shapely.constrained_delaunay_triangles(s));result=[]
        for tri in triangles:
            coords=shapely.get_coordinates(tri)[:-1]
            if len(coords)!=3 or any(tuple(q)not in index for q in coords):
                raise ValueError('Constrained support triangulation introduced non-source controls')
            result.append(np.array([index[tuple(q)]for q in coords]))
        if not result or not np.isclose(sum(shapely.area(t)for t in triangles),shapely.area(s),rtol=1e-10,atol=1e-18):
            raise ValueError('Incomplete source support partition')
        return result

    result=[]
    for i in parts(a):
        for j in parts(b):
            if time.monotonic()-start>time_limit_seconds:raise TimeoutError('Native convex support decomposition deadline')
            sep=support(a[i],b[j])
            if sep is None:raise ValueError('Source domains do not establish requested support clearance')
            result.append((i,j,*sep))
            if len(result)>4096:raise ValueError('Explicit support-pair budget exceeded')
    return result
