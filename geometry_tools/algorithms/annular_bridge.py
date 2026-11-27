"""Join declared finite-wall mouths without moving controls or plugging lumens.

Unequal genuine polygon cycles are connected by a common arc parameter;
interior stations follow a rounded Hermite course with transported frames.
This helper supplies no anatomical mouth identification or acceptance.
"""
from collections import defaultdict
import numpy as np
from .tube_continuity import boundary_loops
from .tube_sections import rotation_minimizing_frames

def mouth_loops(vertices,faces):
    v=np.asarray(vertices,float);f=np.asarray(faces,int);loops=boundary_loops(f)
    if len(loops)!=2:raise ValueError('A tissue mouth needs exactly two boundary cycles')
    tri=v[f];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]).sum(0);normal/=np.linalg.norm(normal)
    center=v[np.unique(f)].mean(0)
    def area(ring):
        p=v[ring]-center;return abs(np.cross(p,np.roll(p,-1,axis=0)).sum(0)@normal)
    outer=max(loops,key=area);inner=next(ring for ring in loops if ring is not outer)
    # Boundary loop of an outward annular chart has opposite inner winding.
    if np.dot(np.cross(v[outer]-center,np.roll(v[outer],-1,axis=0)-center).sum(0),normal)<0:outer=outer[::-1]
    if np.dot(np.cross(v[inner]-center,np.roll(v[inner],-1,axis=0)-center).sum(0),normal)<0:inner=inner[::-1]
    return np.asarray(outer),np.asarray(inner),normal

def find_annular_patch(vertices,faces,center,normal,*,maximum_distance=.012):
    """Find an actual two-cycle tissue end near an explicitly named source.

    A nearby sidewall or a solid end disk is rejected by the two-cycle test.
    The caller still owns its anatomical identity and recipient assignment.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    v,inv=np.unique(np.asarray(vertices,float),axis=0,return_inverse=True);f=inv[np.asarray(faces,int)];tri=v[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);fn/=np.linalg.norm(fn,axis=1)[:,None]
    normal=np.asarray(normal,float);normal/=np.linalg.norm(normal);center=np.asarray(center,float)
    selected=np.flatnonzero((abs(fn@normal)>.75)&(np.linalg.norm(tri.mean(1)-center,axis=1)<maximum_distance));sf=f[selected];edges=defaultdict(list)
    for i,t in enumerate(sf):
        for x,y in zip(t,np.roll(t,-1)):edges[tuple(sorted((int(x),int(y))))].append(i)
    pairs=[p for p in edges.values()if len(p)==2]
    if not len(selected)or not pairs:raise ValueError('No annular end chart at the declared source')
    graph=coo_matrix((np.ones(len(pairs)*2),(np.r_[[p[0]for p in pairs],[p[1]for p in pairs]],np.r_[[p[1]for p in pairs],[p[0]for p in pairs]])),shape=(len(sf),len(sf))).tocsr();_,labels=connected_components(graph,directed=False);choices=[]
    for label in set(labels):
        take=np.flatnonzero(labels==label)
        if len(take)<6:continue
        try:outer,inner,n=mouth_loops(v,sf[take])
        except ValueError:continue
        cc=v[np.unique(sf[take])].mean(0);choices.append((np.linalg.norm(cc-center),selected[take]))
    if not choices:raise ValueError('No two-cycle annular end; a sidewall or solid disk cannot be a conduit mouth')
    return min(choices,key=lambda x:x[0])[1]

def trim_named_inlet(vertices,faces,plane_y,seed):
    """Remove only the declared inlet-side connected component of a cut.

    Other caudally directed branches are retained and their temporary cut
    faces cancel. The returned inlet annulus is a real wall section, not a
    cap across fluid. The source and output still require geometry audits.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from .shared_closed_cell_clip import SharedCutBank,split_closed_cell
    v,inv=np.unique(np.asarray(vertices,float),axis=0,return_inverse=True);f=inv[np.asarray(faces,int)]
    bank=SharedCutBank(v);cut=split_closed_cell(bank,f,[0,1,0],plane_y,'declared-inlet',allow_holes=True);v=np.asarray(bank.positions);nf=cut.negative_faces
    e=np.r_[nf[:,[0,1]],nf[:,[1,2]],nf[:,[2,0]]];graph=coo_matrix((np.ones(len(e)*2),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(v),len(v))).tocsr();_,labels=connected_components(graph,directed=False)
    used=np.unique(nf);point=used[np.argmin(np.linalg.norm(v[used]-np.asarray(seed),axis=1))];discard=labels[point];drop=labels[nf[:,0]]==discard
    removed={tuple(sorted(map(int,t)))for t in nf[drop]};retained=np.r_[cut.positive_faces,nf[~drop]];counts=defaultdict(int)
    for t in retained:counts[tuple(sorted(map(int,t)))]+=1
    final=np.array([t for t in retained if counts[tuple(sorted(map(int,t)))]==1]);keys={tuple(sorted(map(int,t)))for t in cut.interface_faces if tuple(sorted(map(int,t)))in removed};cap=np.array([i for i,t in enumerate(final)if tuple(sorted(map(int,t)))in keys])
    mouth_loops(v,final[cap]) # explicitly reject a solid disk or branched rim
    return v,final,cap,dict(discarded_inlet_triangles=int(drop.sum()),other_caudal_components_retained=True,plane_y_m=plane_y)

def arc(points):
    length=np.linalg.norm(np.roll(points,-1,axis=0)-points,axis=1)
    if np.any(length<=0):raise ValueError('Collapsed genuine mouth edge')
    return np.r_[0,np.cumsum(length)]/length.sum()

def sample(points,t):
    u=arc(points);i=np.minimum(np.searchsorted(u,t,side='right')-1,len(points)-1);w=(t-u[i])/(u[i+1]-u[i])
    return points[i]*(1-w[:,None])+points[(i+1)%len(points)]*w[:,None]

def zipper(a,b,va,vb,*,angles_a=None,angles_b=None):
    if angles_a is not None and angles_b is not None:
        aa=np.mod(angles_a,1);bb=np.mod(angles_b,1)
        ia=np.argsort(aa);ib=np.argsort(bb);a=np.asarray(a)[ia];b=np.asarray(b)[ib];aa=aa[ia];bb=bb[ib]
        i=j=0;x=a[-1];y=b[-1];f=[]
        while i<len(a) or j<len(b):
            ta=aa[i]if i<len(a)else np.inf;tb=bb[j]if j<len(b)else np.inf
            if abs(ta-tb)<1e-12:
                f.extend([[x,a[i],b[j]],[x,b[j],y]]);x=a[i];y=b[j];i+=1;j+=1
            elif ta<tb:f.append([x,a[i],y]);x=a[i];i+=1
            else:f.append([x,b[j],y]);y=b[j];j+=1
        return np.asarray(f,int)
    aa=arc(va)if angles_a is None else np.r_[angles_a,angles_a[0]+1]
    bb=arc(vb)if angles_b is None else np.r_[angles_b,angles_b[0]+1]
    i=j=0;f=[];na,nb=len(a),len(b)
    while i<na or j<nb:
        x,y=a[i%na],b[j%nb];ta=aa[i+1]if i<na else np.inf;tb=bb[j+1]if j<nb else np.inf
        if abs(ta-tb)<1e-12:
            f.extend([[x,a[(i+1)%na],b[(j+1)%nb]],[x,b[(j+1)%nb],y]]);i+=1;j+=1
        elif ta<tb:f.append([x,a[(i+1)%na],y]);i+=1
        else:f.append([x,b[(j+1)%nb],y]);j+=1
    return np.asarray(f,int)

def radial_profile(points,angles):
    """Common angular rays preserve nesting of an irregular thin wall."""
    p=np.asarray(points,float);e=np.roll(p,-1,axis=0)-p;result=[]
    def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]
    for angle in angles:
        d=np.array([np.cos(angle),np.sin(angle)]);denom=cross(d,e[:,:2]);t=np.divide(cross(p[:,:2],d),denom,out=np.full(len(p),np.nan),where=abs(denom)>1e-18);radius=np.divide(cross(p[:,:2],e[:,:2]),denom,out=np.full(len(p),np.nan),where=abs(denom)>1e-18);hit=np.flatnonzero((t>=-1e-10)&(t<=1+1e-10)&(radius>0))
        if not len(hit)or np.ptp(radius[hit])>max(1e-10,float(radius[hit].max())*1e-8):raise ValueError('A genuine mouth is not star-shaped about the selected fluid centre')
        j=int(hit[0]);result.append(p[j]+np.clip(t[j],0,1)*e[j])
    return np.asarray(result)

def fluid_kernel_center(v,outer,inner,normal):
    """Choose a point in the visibility kernel of both genuine rim polygons."""
    from scipy.optimize import linprog
    c=v[inner].mean(0);u=v[outer[0]]-c;u-=normal*np.dot(u,normal);u/=np.linalg.norm(u);w=np.cross(normal,u);basis=np.column_stack([u,w]);rows=[];limits=[];scale=float(np.max(np.linalg.norm((v[outer]-c)@basis,axis=1)))
    for ring in [outer,inner]:
        p=((v[ring]-c)@basis)/scale;e=np.roll(p,-1,axis=0)-p
        nxt=np.roll(p,-1,axis=0)
        if np.sum(p[:,0]*nxt[:,1]-p[:,1]*nxt[:,0])<0:p=p[::-1];e=np.roll(p,-1,axis=0)-p
        for x,d in zip(p,e):
            d=d/np.linalg.norm(d);rows.append([d[1],-d[0],1.]);limits.append(d[1]*x[0]-d[0]*x[1])
    q=linprog([0,0,-1],A_ub=rows,b_ub=limits,bounds=[(None,None),(None,None),(0,None)],method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not q.success or q.x[2]<1e-10:raise ValueError('No common visibility kernel for the genuine fluid rim')
    return c+basis@(q.x[:2]*scale)

def attach_bridge(parent_v,parent_f,parent_cap,child_v,child_f,child_cap,*,stations=12,controls=24,straight_shear=False):
    """Return the parent plus a finite collar, meeting the unchanged child.

    The returned final annulus is literally the child's reversed tissue cap.
    Hide that reciprocal chart only when the adjoining child is displayed.
    """
    pv,pf=np.asarray(parent_v,float),np.asarray(parent_f,int);cv,cf=np.asarray(child_v,float),np.asarray(child_f,int)
    pc,cc=np.asarray(parent_cap,int),np.asarray(child_cap,int)
    po,pi,n0=mouth_loops(pv,pf[pc]);co,ci,nc=mouth_loops(cv,cf[cc]);n1=-nc
    c0=fluid_kernel_center(pv,po,pi,n0);c1=fluid_kernel_center(cv,co,ci,nc);distance=np.linalg.norm(c1-c0)
    if distance<1e-7:raise ValueError('Coincident mouths need a literal join, not a collar')
    co,ci=co[::-1],ci[::-1]
    # Keep the nearest cyclic correspondence of the two real inner/outer rims.
    def align_inner(v,outer,inner):
        costs=[np.mean(np.linalg.norm(v[outer]-sample(v[np.roll(inner,k)],arc(v[outer])[:-1]),axis=1))for k in range(len(inner))]
        return np.roll(inner,int(np.argmin(costs)))
    pi=align_inner(pv,po,pi);ci=align_inner(cv,co,ci)
    t=np.linspace(0,1,stations);length=min(distance*.5,.015);b0=c0;b1=c0+n0*length;b2=c1+nc*length;b3=c1
    path=(1-t[:,None])**3*b0+3*(1-t[:,None])**2*t[:,None]*b1+3*(1-t[:,None])*t[:,None]**2*b2+t[:,None]**3*b3
    tangents=3*(1-t[:,None])**2*(b1-b0)+6*(1-t[:,None])*t[:,None]*(b2-b1)+3*t[:,None]**2*(b3-b2);tangents/=np.linalg.norm(tangents,axis=1)[:,None]
    frames=rotation_minimizing_frames(path,tangents=tangents,initial_reference=pv[po[0]]-c0)
    if np.dot(n0,n1)>.99999:
        # Parallel anatomical sections admit a shear loft. Tilting these broad,
        # thin rings to a short diagonal centreline can make a valid inner wall
        # cross its outer wall on the concave bend. Preserve their common
        # section planes and strictly ordered axial stations instead.
        u=pv[po[0]]-c0;u-=n0*np.dot(u,n0);u/=np.linalg.norm(u)
        frames[:]=np.column_stack([u,np.cross(n0,u),n0])
    if straight_shear:
        u=pv[po[0]]-c0;u-=n1*np.dot(u,n1);u/=np.linalg.norm(u)
        frames[:]=np.column_stack([u,np.cross(n1,u),n1]);path=(1-t[:,None])*c0+t[:,None]*c1
    def local(v,ids,c,frame):return(v[ids]-c)@frame
    pprofile=local(pv,po,c0,frames[0]);cprofile=local(cv,co,c1,frames[-1])
    candidates=[]
    for k in range(len(co)):
        q=sample(np.roll(cprofile,k,axis=0),np.arange(controls)/controls)
        candidates.append(np.sum((sample(pprofile,np.arange(controls)/controls)[:,:2]-q[:,:2])**2))
    k=int(np.argmin(candidates));co=np.roll(co,k);ci=np.roll(ci,k)
    vertices=np.r_[pv,cv].tolist();offset=len(pv);layers=[];ring_angles=[]
    for inner,pids,cids in [(False,po,co),(True,pi,ci)]:
        rings=[pids.copy()];angles=np.arange(controls)*2*np.pi/controls;lp=local(pv,pids,c0,frames[0]);lc=local(cv,cids,c1,frames[-1]);p0=radial_profile(lp,angles);p1=radial_profile(lc,angles)
        def parameters(profile):
            q=np.unwrap(np.arctan2(profile[:,1],profile[:,0]))/(2*np.pi)
            if np.any(np.diff(q)<=0):raise ValueError('A mouth cycle reverses angular order: '+str(dict(minimum=float(np.diff(q).min()),axis_dot=float(n0@n1),profile=profile.tolist())))
            return q
        aa=[parameters(lp)]
        for i,s in enumerate(t[1:-1],1):
            profile=(1-s)*p0+s*p1
            if not straight_shear:profile[:,2]=(1-s)**3*p0[:,2]+s**3*p1[:,2]
            points=path[i]+profile@frames[i].T;ring=np.arange(len(vertices),len(vertices)+controls);vertices.extend(points.tolist());rings.append(ring);aa.append(np.arange(controls)/controls)
        rings.append(cids+offset);aa.append(parameters(lc));layers.append(rings);ring_angles.append(aa)
    v=np.asarray(vertices,float);outer=[];inner=[]
    for rings,aa,target in zip(layers,ring_angles,[outer,inner]):
        for i,(a,b)in enumerate(zip(rings[:-1],rings[1:])):target.extend(zipper(a,b,v[a],v[b],angles_a=aa[i],angles_b=aa[i+1]).tolist())
    outer=np.asarray(outer,int);inner=np.asarray(inner,int)[:,[0,2,1]];end=cf[cc][:,[0,2,1]]+offset
    keep=np.ones(len(pf),bool);keep[pc]=False;rest=pf[keep];faces=np.r_[rest,outer,inner,end]
    return dict(vertices=v,faces=faces,parent_face_ids=np.flatnonzero(keep),outer_faces=outer,inner_faces=inner,interface_faces=np.arange(len(rest)+len(outer)+len(inner),len(faces)),child_interface_faces=cc,child_vertex_offset=offset,
                parent_source_controls_literal=bool(np.array_equal(v[:len(pv)],pv)),child_source_controls_literal=bool(np.array_equal(v[offset:offset+len(cv)],cv)),path=path.tolist(),source_cycle_controls=[len(po),len(co)],lumen_disks_added=0)
