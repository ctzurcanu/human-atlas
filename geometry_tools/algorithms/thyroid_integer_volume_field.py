"""Staged thyroid-only deterministic microscopic volume candidate.

PCG3D comes from Jarzynski/Olano, JCGT 9(3), 2020, section 6.1.
Only the unsigned hash and fixed-point site construction are bit-defined.
World transforms, histology and a live renderer are separate acceptance gates.
"""
from collections import defaultdict
import heapq
import numpy as np
from scipy.spatial import cKDTree

PITCH_NM=353554  # Rounded upward from 500 um / sqrt(2).
LARGE_CONTAINER_NM=250000
PADDING_NM=2000
MINIMUM_RADIUS_NM=25000
SIZE_OFFSET=np.array([1913,373,1117],dtype=np.uint32)

GLSL_SITE=r'''
// Candidate only: PCG3D, Jarzynski/Olano JCGT 2020, section 6.1.
// The CPU oracle uses these same unsigned operations modulo 2^32.
uvec3 thyroidIntegerHash(uvec3 v){
    v=v*1664525u+1013904223u;
    v.x+=v.y*v.z;v.y+=v.z*v.x;v.z+=v.x*v.y;
    v^=v>>16u;
    v.x+=v.y*v.z;v.y+=v.z*v.x;v.z+=v.x*v.y;
    return v;
}
ivec4 thyroidIntegerSite(ivec3 id){
    uvec3 random=thyroidIntegerHash(uvec3(id));
    uvec3 size=thyroidIntegerHash(uvec3(id)+uvec3(1913u,373u,1117u));
    bool large=((id.x+id.y+id.z)&1)==0;
    uint container=large?250000u:101554u;
    uint span=container-2000u-25000u+1u;
    uint radius=25000u+max(max(size.x%span,size.y%span),size.z%span);
    // A conservative cube inside the remaining spherical container.
    uint jitter=(container-radius-2000u)*577u/1000u;
    uint width=2u*jitter+1u;
    ivec3 offset=ivec3(random%uvec3(width))-ivec3(int(jitter));
    return ivec4(id*353554+offset,int(radius));
}
'''

GLSL_LOOKUP_FIELD=r'''
uniform highp usampler3D thyroidOwnerPages;
uniform highp usampler2D thyroidOwnerBricks;
uniform ivec3 thyroidBrickLow;
uniform ivec3 thyroidPageDimensions;
uint thyroidSiteGroup(ivec3 id){
    ivec3 brick=ivec3(floor(vec3(id)/8.0));
    ivec3 page=brick-thyroidBrickLow;
    if(any(lessThan(page,ivec3(0)))||any(greaterThanEqual(page,thyroidPageDimensions)))return 0u;
    uint address=texelFetch(thyroidOwnerPages,page,0).r;
    if(address==0u)return 0u;
    ivec3 local=id-brick*8;
    int word=local.x*64+local.y*8+local.z;
    return texelFetch(thyroidOwnerBricks,ivec2(word,int(address)-1),0).r;
}
// q is the declared organ-local rotated physical frame, in metres.
// group0 is interfollicular stroma. The discrete group lookup identifies every
// fully enclosed follicle, without allocating an individual mesh/draw call.
uvec2 thyroidVolumeOwner(vec3 q){
    ivec3 base=ivec3(floor(q/0.000353554));
    for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++)for(int z=-1;z<=1;z++){
        ivec3 id=base+ivec3(x,y,z);uint group=thyroidSiteGroup(id);
        if(group==0u)continue;
        ivec4 site=thyroidIntegerSite(id);
        vec3 delta=q-vec3(site.xyz)*1e-9;float radius=float(site.w)*1e-9;
        if(dot(delta,delta)<radius*radius){
            float coreRadius=radius*.94;
            return uvec2(group,dot(delta,delta)<coreRadius*coreRadius?1u:2u);
        }
    }
    return uvec2(0u,0u);
}
'''


def pcg3d(ids):
    """Literal sequential unsigned operations; signed input wraps modulo2^32."""
    v=np.asarray(ids,dtype=np.int64).astype(np.uint32).copy()
    v=v*np.uint32(1664525)+np.uint32(1013904223)
    for phase in range(2):
        v[:,0]+=v[:,1]*v[:,2]
        v[:,1]+=v[:,2]*v[:,0]
        v[:,2]+=v[:,0]*v[:,1]
        if phase==0:v^=v>>np.uint32(16)
    return v


def integer_sites(ids):
    ids=np.asarray(ids,dtype=np.int64)
    if np.max(abs(ids))*PITCH_NM>np.iinfo(np.int32).max-250000:
        raise ValueError('Local thyroid site coordinate would overflow int32')
    random=pcg3d(ids);size=pcg3d(ids+SIZE_OFFSET.astype(np.int64))
    large=(ids.sum(1)&1)==0
    container=np.where(large,LARGE_CONTAINER_NM,PITCH_NM-LARGE_CONTAINER_NM-PADDING_NM)
    span=(container-PADDING_NM-MINIMUM_RADIUS_NM+1).astype(np.uint32)
    radius=MINIMUM_RADIUS_NM+np.max(size%span[:,None],axis=1).astype(np.int64)
    jitter=(container-radius-PADDING_NM)*577//1000
    offset=(random%(2*jitter+1).astype(np.uint32)[:,None]).astype(np.int64)-jitter[:,None]
    centers=ids*PITCH_NM+offset
    assert np.all(np.linalg.norm(offset,axis=1)+radius+PADDING_NM<=container+1e-8)
    return np.column_stack([centers,radius]).astype('<i4'),large


def surface_band_and_interior_sites(triangles_lattice,margin=.75):
    """Conservative lattice enumeration; no full dense volume tessellation.

Ray intervals include interior grid nodes. Expanded source-triangle boxes
include every possible jittered-inside site whose nominal node is outside.
Final center/sphere containment is checked separately against finite triangles.
"""
    tri=np.asarray(triangles_lattice,dtype=float)
    low=np.floor(tri.min((0,1))-margin).astype(int)
    high=np.ceil(tri.max((0,1))+margin).astype(int)
    shape=high-low+1
    if np.prod(shape)>20000000:raise ValueError('Unexpected thyroid lattice bounds')
    mask=np.zeros(shape,dtype=bool);events=defaultdict(list)
    for t in tri:
        band_low=np.maximum(0,np.floor(t.min(0)-margin).astype(int)-low)
        band_high=np.minimum(shape-1,np.ceil(t.max(0)+margin).astype(int)-low)
        mask[tuple(slice(a,b+1) for a,b in zip(band_low,band_high))]=True
        a,b,c=t;u=b[1:]-a[1:];v=c[1:]-a[1:]
        determinant=u[0]*v[1]-u[1]*v[0]
        if abs(determinant)<1e-12:continue
        yl,yh=int(np.ceil(t[:,1].min())),int(np.floor(t[:,1].max()))
        zl,zh=int(np.ceil(t[:,2].min())),int(np.floor(t[:,2].max()))
        for y in range(yl,yh+1):
            for z in range(zl,zh+1):
                d=np.array([y,z])-a[1:]
                s=(d[0]*v[1]-d[1]*v[0])/determinant
                w=(u[0]*d[1]-u[1]*d[0])/determinant
                if min(s,w,1-s-w)<-1e-12:continue
                x=a[0]+s*(b[0]-a[0])+w*(c[0]-a[0])
                events[(y,z)].append(float(x))
    odd=0
    for (y,z),hits in events.items():
        hits.sort();unique=[]
        for x in hits:
            if not unique or abs(x-unique[-1])>1e-8:unique.append(x)
        if len(unique)%2:
            # The conservative triangle band already includes tangencies;
            # do not silently synthesize an interior interval for an odd ray.
            odd+=1;continue
        for a,b in zip(unique[::2],unique[1::2]):
            first=max(0,int(np.ceil(a))-low[0]);last=min(shape[0]-1,int(np.floor(b))-low[0])
            if last>=first:mask[first:last+1,y-low[1],z-low[2]]=True
    return np.argwhere(mask).astype(np.int64)+low,{'boundsLow':low.tolist(),'boundsHigh':high.tolist(),
                                                 'conservativeCandidates':int(mask.sum()),'gridPointsInBounds':int(mask.size),
                                                 'ambiguousRegularGridRays':odd}


def ray_contains(points,triangles,bucket_width=.0007):
    """Oriented finite X-ray occupancy with explicit edge ambiguities."""
    points=np.asarray(points,dtype=float);tri=np.asarray(triangles,dtype=float)
    origin=tri.mean((0,1));p=points-origin;tri=tri-origin
    buckets=defaultdict(list)
    for i,key in enumerate(np.floor(p[:,1:]/bucket_width).astype(int)):
        buckets[tuple(key)].append(i)
    buckets={key:np.array(value,dtype=int) for key,value in buckets.items()}
    count=np.zeros(len(points),dtype=int);ambiguous=np.zeros(len(points),dtype=bool)
    for a,b,c in tri:
        u=b-a;v=c-a;det=u[1]*v[2]-u[2]*v[1]
        if abs(det)<1e-20:continue
        lower=np.floor(np.minimum.reduce([a[1:],b[1:],c[1:]])/bucket_width).astype(int)
        upper=np.floor(np.maximum.reduce([a[1:],b[1:],c[1:]])/bucket_width).astype(int)
        pieces=[buckets[(y,z)] for y in range(lower[0],upper[0]+1) for z in range(lower[1],upper[1]+1) if (y,z) in buckets]
        if not pieces:continue
        ids=np.concatenate(pieces);d=p[ids]-a
        s=(d[:,1]*v[2]-d[:,2]*v[1])/det;w=(u[1]*d[:,2]-u[2]*d[:,1])/det
        margin=np.minimum.reduce([s,w,1-s-w])
        x=a[0]+s*u[0]+w*v[0];hit=(margin>=-1e-11)&(x>p[ids,0]+1e-12)
        ambiguous[ids[(abs(margin)<1e-9)&hit]]=True
        count[ids[hit]]+=1 if det>0 else -1
    if np.any((count!=0)&(count!=1)&~ambiguous):raise ValueError('Inconsistent oriented ray count')
    return count==1,ambiguous


def near_boundary_distance(points,triangles,maximum):
    """Exact finite point/triangle distance where within maximum; else inf."""
    points=np.asarray(points,dtype=float);tree=cKDTree(points);result=np.full(len(points),np.inf)
    for tri in triangles:
        a,b,c=np.asarray(tri,dtype=float);center=(a+b+c)/3
        candidates=tree.query_ball_point(center,np.linalg.norm(tri-center,axis=1).max()+maximum)
        if not candidates:continue
        ids=np.asarray(candidates);p=points[ids];d=p-a;u=b-a;v=c-a
        uu,uv,vv=u@u,u@v,v@v;det=uu*vv-uv*uv
        if det<=0:raise ValueError('Degenerate finite triangle')
        du,dv=d@u,d@v;s=(du*vv-dv*uv)/det;w=(dv*uu-du*uv)/det
        n=np.cross(u,v);squared=np.where((s>=0)&(w>=0)&(s+w<=1),(d@n)**2/(n@n),np.inf)
        for start,end in [(a,b),(b,c),(c,a)]:
            edge=end-start;t=np.clip((p-start)@edge/(edge@edge),0,1)
            delta=p-start-t[:,None]*edge;squared=np.minimum(squared,np.einsum('ij,ij->i',delta,delta))
        result[ids]=np.minimum(result[ids],np.sqrt(squared))
    return result


class PowerOwnership:
    """Full 3D variable-radius ownership, not only exterior/cap coloring."""
    def __init__(self,centers,radii,groups):
        self.centers=np.asarray(centers,dtype=float);self.radii=np.asarray(radii,dtype=float)
        self.groups=np.asarray(groups);self.tree=cKDTree(self.centers);self.maximum_radius=float(self.radii.max())

    def query(self,points):
        points=np.asarray(points,dtype=float);distance,index=self.tree.query(points)
        power=distance**2-self.radii[index]**2
        bounds=np.sqrt(np.maximum(0,power+self.maximum_radius**2))+1e-14
        candidates=self.tree.query_ball_point(points,bounds)
        for i,ids in enumerate(candidates):
            ids=np.asarray(ids);d=self.centers[ids]-points[i];value=np.einsum('ij,ij->i',d,d)-self.radii[ids]**2
            rank=np.lexsort((ids,value));index[i]=ids[rank[0]];power[i]=value[rank[0]]
        return index,self.groups[index],power

    def segment_clearance(self,segments,outer_radius=0.,margin=0.):
        """Complete swept straight-segment/sphere exclusion, not sampled rays.

A query ball encloses the complete segment Minkowski radius. Exact clamped
projection then checks every potentially contacting spherical follicle. Curve
chords require an additional caller-supplied deviation bound. This gate alone
does not establish that a vascular wall remains inside the organ or connects
to another vessel with correct shared annuli.
"""
        segments=np.asarray(segments,dtype=float)
        if segments.ndim!=3 or segments.shape[1:]!=(2,3):raise ValueError('Expected Nx2x3 segments')
        radius=np.broadcast_to(np.asarray(outer_radius,dtype=float),len(segments))
        if np.any(radius<0) or margin<0:raise ValueError('Radii and margin must be nonnegative')
        rows=[]
        for (a,b),r in zip(segments,radius):
            delta=b-a;length2=delta@delta
            candidates=self.tree.query_ball_point((a+b)/2,np.sqrt(length2)/2+self.maximum_radius+r+margin)
            if not candidates:
                rows.append({'potentiallyContactingFollicles':0,'pass':True,'minimumDetectedGapMeters':None});continue
            ids=np.asarray(candidates);centers=self.centers[ids]
            t=np.clip((centers-a)@delta/length2,0,1) if length2>0 else np.zeros(len(ids))
            points=a+t[:,None]*delta
            gap=np.linalg.norm(centers-points,axis=1)-self.radii[ids]-r
            k=int(np.argmin(gap));violating=ids[gap<=margin]
            rows.append({'potentiallyContactingFollicles':len(ids),'pass':not len(violating),
                         'minimumDetectedGapMeters':float(gap[k]),'nearestSiteIndex':int(ids[k]),
                         'segmentFraction':float(t[k]),'violatingSiteIndices':violating.tolist()})
        return rows


def power_neighbor_graph(centers,radii,neighbors=10):
    """Conservative true-power-face graph, not arbitrary proximity links.

At a pair's bisector along the center segment, no third site may have lower
power. The disjoint spheres make both centers strictly interior to their
power cells, so an accepted pair has a shared-face connection. Intersection
with the concave organ can still change connectivity and needs extraction QA.
"""
    centers=np.asarray(centers);radii=np.asarray(radii)
    tree=cKDTree(centers);_,knn=tree.query(centers,k=min(neighbors+1,len(centers)))
    a=np.repeat(np.arange(len(centers)),knn.shape[1]-1);b=knn[:,1:].reshape(-1)
    pairs=np.unique(np.sort(np.c_[a,b],axis=1),axis=0);a,b=pairs.T
    delta=centers[b]-centers[a];d2=np.einsum('ij,ij->i',delta,delta)
    t=(d2+radii[a]**2-radii[b]**2)/(2*d2)
    if np.any((t<=0)|(t>=1)):raise ValueError('Sites violate disjoint interior ownership')
    x=centers[a]+t[:,None]*delta
    pair_power=np.einsum('ij,ij->i',x-centers[a],x-centers[a])-radii[a]**2
    oracle=PowerOwnership(centers,radii,np.arange(len(centers)))
    selected=[]
    for first in range(0,len(pairs),4096):
        owner,_,actual=oracle.query(x[first:first+4096])
        take=actual>=pair_power[first:first+4096]-1e-16
        selected.extend(pairs[first:first+4096][take].tolist())
    graph=[[] for _ in centers]
    for a,b in selected:graph[a].append(b);graph[b].append(a)
    return graph,{'proximityPairsTested':len(pairs),'sharedPowerFaceEdges':len(selected),
                  'unclippedPowerAdjacencyVerified':True,'concaveOrganClippedConnectivityVerified':False}


def connected_count_groups(graph,centers,minimum=20,maximum=40):
    """Count-constrained connected site groups on actual power-face edges.

Small end remainders merge through a real graph edge when the cap permits.
An infeasible graph component remains pending; it is never silently discarded
or renamed a lobule. The graph is not an anatomical reconstruction.
"""
    centers=np.asarray(centers);n=len(graph);parent=np.full(n,-2,dtype=int);order=[];roots=[]
    # A depth-first tree keeps dense 3D graph branches long rather than using
    # one center fan or rectilinear fixed-size boxes.
    for root in np.lexsort(centers.T[::-1]):
        if parent[root]!=-2:continue
        roots.append(int(root));parent[root]=-1;order.append(int(root))
        def neighbors(v):return iter(sorted(graph[v],key=lambda j:float(np.linalg.norm(centers[j]-centers[v]))))
        stack=[(int(root),neighbors(root))]
        while stack:
            v,iterator=stack[-1]
            try:j=next(iterator)
            except StopIteration:stack.pop();continue
            if parent[j]==-2:
                parent[j]=v;order.append(j);stack.append((j,neighbors(j)))
    children=[[] for _ in graph]
    for i,p in enumerate(parent):
        if p>=0:children[p].append(i)
    owners=np.zeros(n,dtype=np.uint16);pending={};groups=[];oversize=[]
    for v in reversed(order):
        members=[v]
        for child in children[v]:members.extend(pending.pop(child,[]))
        if len(members)<minimum:pending[v]=members;continue
        if len(members)>maximum:
            oversize.append({'treeRoot':v,'sites':len(members)});pending[v]=members;continue
        group=len(groups)+1
        if group>65535:raise ValueError('u16 groups overflow')
        owners[members]=group;groups.append(np.asarray(members,dtype=int))
    # A remaining small component is attached by a literal graph edge; avoid
    # filling an unsupported isolated branch with an arbitrary remote owner.
    unresolved=[]
    for root,members in pending.items():
        adjacent={int(owners[j]) for i in members for j in graph[i] if owners[j]}
        candidates=[g for g in adjacent if len(groups[g-1])+len(members)<=maximum]
        if candidates:
            group=min(candidates,key=lambda g:float(np.linalg.norm(centers[groups[g-1]].mean(0)-centers[members].mean(0))))
            groups[group-1]=np.r_[groups[group-1],members];owners[members]=group
        else:unresolved.append({'treeRoot':int(root),'sites':len(members),'adjacentGroups':sorted(adjacent)})
    # The still-connected tree backbone can be much larger than one end
    # remainder. Absorb it through proven power edges, keeping each group
    # connected and capped. This uses no arbitrary spatial-plane assignment.
    before_flood=int((owners==0).sum());sizes=np.array([len(g) for g in groups],dtype=int)
    seeds=np.array([centers[g].mean(0) for g in groups]);heap=[];appended=[[] for _ in groups]
    def push(point,group):
        if sizes[group-1]<maximum:
            distance=float(np.linalg.norm(centers[point]-seeds[group-1]))
            heapq.heappush(heap,(distance,point,group))
    for i in np.flatnonzero(owners==0):
        for group in {int(owners[j]) for j in graph[i] if owners[j]}:push(int(i),group)
    while heap:
        distance,point,group=heapq.heappop(heap)
        if owners[point] or sizes[group-1]>=maximum:continue
        owners[point]=group;sizes[group-1]+=1;appended[group-1].append(point)
        for j in graph[point]:
            if not owners[j]:push(j,group)
    groups=[np.r_[g,added].astype(int) for g,added in zip(groups,appended)]
    # Repair the few capacity-saturated pockets with connected subtree cuts of
    # a small induced graph. A subtree and its complement are both connected
    # by literal power-face edges, so repair cannot introduce remote owners.
    repairs=[]
    def partition_pool(pool):
        remaining=set(map(int,pool));result=[]
        while len(remaining)>maximum:
            chosen=None
            for attempt in range(min(8,len(remaining))):
                root=sorted(remaining,key=lambda i:tuple(centers[i]))[attempt]
                par={root:-1};sequence=[root]
                def nearby(i):return iter(sorted((j for j in graph[i] if j in remaining),
                                                 key=lambda j:float(np.linalg.norm(centers[j]-centers[i]))))
                stack=[(root,nearby(root))]
                while stack:
                    i,it=stack[-1]
                    try:j=next(it)
                    except StopIteration:stack.pop();continue
                    if j not in par:par[j]=i;sequence.append(j);stack.append((j,nearby(j)))
                if len(sequence)!=len(remaining):return None
                sub={i:[i] for i in sequence};options=[]
                for i in reversed(sequence):
                    if i!=root:
                        if minimum<=len(sub[i])<=maximum and len(remaining)-len(sub[i])>=minimum:
                            q=centers[sub[i]];diameter=float(np.linalg.norm(q.max(0)-q.min(0)))
                            options.append((diameter,abs(len(sub[i])-30),i))
                        sub[par[i]].extend(sub[i])
                if options:
                    _,_,i=min(options);chosen=sub[i];break
            if chosen is None:return None
            result.append(np.asarray(chosen,dtype=int));remaining.difference_update(chosen)
        if len(remaining)<minimum:return None
        result.append(np.asarray(sorted(remaining),dtype=int));return result
    todo=set(np.flatnonzero(owners==0));pockets=[]
    while todo:
        start=todo.pop();pocket=[start];queue=[start]
        while queue:
            for j in graph[queue.pop()]:
                if j in todo:todo.remove(j);pocket.append(j);queue.append(j)
        pockets.append(pocket)
    for pocket in sorted(pockets,key=len,reverse=True):
        adjacent={int(owners[j]) for i in pocket for j in graph[i] if owners[j]}
        used=[];pool=list(pocket);parts=partition_pool(pool) if len(pool)>=minimum else None
        if parts is None:
            ordered=sorted(adjacent,key=lambda g:float(np.linalg.norm(centers[groups[g-1]].mean(0)-centers[pocket].mean(0))))
            for g in ordered[:8]:
                used.append(g);pool.extend(groups[g-1].tolist());parts=partition_pool(pool)
                if parts is not None:break
        if parts is None:continue
        if len(parts)<len(used):raise ValueError('Unexpected group-capacity reduction')
        for i,part in enumerate(parts):
            if i<len(used):group=used[i];groups[group-1]=part
            else:group=len(groups)+1;groups.append(part)
            if group>65535:raise ValueError('u16 groups overflow')
            owners[part]=group
        repairs.append({'unassignedSites':len(pocket),'reusedGroups':len(used),'resultGroups':len(parts)})
    left=int((owners==0).sum())
    return owners,groups,{'components':len(roots),'groupCounts':[len(g) for g in groups],
                          'unassignedBeforePowerEdgeFlood':before_flood,'unassignedSites':left,
                          'unresolvedRemainders':[] if not left else unresolved,
                          'oversizeTreeEventCount':len(oversize),'maximumTreeRemainder':max([r['sites'] for r in oversize],default=0),
                          'connectedPocketRepairs':repairs,
                          'connectedPowerSiteGroups':True,'allGroups20To40':not left and all(minimum<=len(g)<=maximum for g in groups)}
