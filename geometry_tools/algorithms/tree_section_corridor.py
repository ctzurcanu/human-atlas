"""Finite configuration search around fixed surfaces on a source-incidence tree.

One translation per original section/core group. No radii, faces or branches
are removed. Vertex clearance is proposal evidence, not full triangle clearance.
"""
import time
import numpy as np
from geometry_tools.open3d_backend import library
from .feature_normals import indexed_feature_normals,inward_sheet_clearance


def source_span_hits(points,index):
    """Native segment screen on every literal edge of posed source triangles.

    Float64 triangle replay remains mandatory. This is a proposal screen,
    including segment interiors, rather than a point-only clearance claim.
    """
    p=np.asarray(points,float)
    start=p[..., [0,1,2],:];end=p[..., [1,2,0],:]
    rays=np.concatenate([(start-index.origin)/index.scale,(end-start)/index.scale],axis=-1)
    hit=index.scene.cast_rays(library().core.Tensor(rays.reshape(-1,6).astype(np.float32)),nthreads=2)['t_hit'].numpy()
    hit=hit.reshape(p.shape[:-2]+(3,))
    return np.any((hit>1e-6)&(hit<1-1e-6),axis=-1)


def source_closure_tree(faces,labels):
    """Keep source-incidence cycles as rigid closures; retain every edge."""
    f=np.asarray(faces,int)
    _, owner = np.unique(labels, return_inverse=True)
    n = int(owner.max())+1
    edges = np.unique(np.sort(owner[f][:, [0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
    edges = edges[edges[:,0] != edges[:,1]]
    original_nodes=n
    if len(edges) != n-1:
        # Keep every cycle as one rigid closure. Dropping a source graph edge
        # would leave a finite junction unconstrained; no tissue is removed.
        neighbours=[[]for _ in range(n)]
        for a,b in edges:neighbours[a].append(int(b));neighbours[b].append(int(a))
        parent={0:-1};visit=[0]
        for a in visit:
            for b in neighbours[a]:
                if b not in parent:parent[b]=a;visit.append(b)
        if len(visit)!=n:raise ValueError('Disconnected source control graph')
        union=np.arange(n)
        def root(k):
            while union[k]!=k:union[k]=union[union[k]];k=int(union[k])
            return k
        for a,b in edges:
            if parent[int(a)]==b or parent[int(b)]==a:continue
            path=[];k=int(a)
            while k>=0:path.append(k);k=parent[k]
            tail=[];k=int(b)
            while k not in path:tail.append(k);k=parent[k]
            cycle=path[:path.index(k)+1]+tail
            r=root(cycle[0])
            for k in cycle:union[root(k)]=r
        compact=np.array([root(k)for k in range(n)])
        _,owner=np.unique(compact[owner],return_inverse=True);n=int(owner.max())+1
        edges=np.unique(np.sort(owner[f][:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
        edges=edges[edges[:,0]!=edges[:,1]]
        if len(edges)!=n-1:raise ValueError('Cyclic source closure unresolved')
    return owner,edges,original_nodes


def tree_section_corridor(vertices, faces, labels, bone, skin, *, wall_seconds=5.,recipient_surfaces=(),span_context=None,
                          guide_vertices=None,grid_spacing_m=.001,local_radius_m=None,protected_raw=()):
    started = time.monotonic()
    v, f = np.asarray(vertices, float), np.asarray(faces, int)
    if not 0 < wall_seconds <= 5 or not np.isfinite(v).all():
        raise ValueError('Finite source geometry and at most five native seconds required')
    owner,edges,original_nodes=source_closure_tree(f,labels);n=int(owner.max())+1
    protected_raw=np.asarray(protected_raw,int)
    if protected_raw.size and (protected_raw.min()<0 or protected_raw.max()>=len(v)):
        raise ValueError('Protected raw points must belong to this source component')
    protected=set(owner[protected_raw].tolist())
    neighbours = [[] for _ in range(n)]
    for a,b in edges:
        neighbours[a].append(int(b));neighbours[b].append(int(a))
    order = [0]; parent = {0:-1}
    for a in order:
        for b in neighbours[a]:
            if b not in parent:parent[b]=a;order.append(b)
    if len(order) != n:
        return v.copy(), dict(complete=False, status='DisconnectedSourceGraph')
    guide_translation=np.zeros_like(v)
    if guide_vertices is None:
        if grid_spacing_m!=.001 or local_radius_m is not None:raise ValueError('Fine local grids require an explicit saved source route')
        grid = np.stack(np.meshgrid(*[np.arange(-4,5)*.001]*3,indexing='ij'),axis=-1).reshape(-1,3)
        grid = grid[np.linalg.norm(grid,axis=1) <= .005+1e-15]
    else:
        guide=np.asarray(guide_vertices,float)
        if guide.shape!=v.shape or not np.isfinite(guide).all()or not 0<grid_spacing_m<=.001 or local_radius_m is None or not grid_spacing_m<=local_radius_m<=.002:
            raise ValueError('Corresponding source guide and finite local-grid spacing/radius required')
        guide_translation=guide-v
        for k in range(n):
            shifts=guide_translation[owner==k]
            if np.max(abs(shifts-shifts[0]))>1e-10:raise ValueError('Guide must preserve each complete source closure')
        steps=int(np.floor(local_radius_m/grid_spacing_m))
        if (2*steps+1)**3>4096:raise ValueError('Local grid exceeds its finite state budget')
        grid=np.stack(np.meshgrid(*[np.arange(-steps,steps+1)*grid_spacing_m]*3,indexing='ij'),axis=-1).reshape(-1,3)
        grid=grid[np.linalg.norm(grid,axis=1)<=local_radius_m+1e-15]
    # Batched actual rigid-group vertices; never a bounding-box corridor.
    p = (v[:,None]+guide_translation[:,None]+grid[None]).reshape(-1,3)
    h,_ = bone.closest_points(p)
    tensor = library().core.Tensor(((p-bone.origin)/bone.scale).astype(np.float32))
    inside = bone.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
    bg = h['distance_m']*np.where(inside,-1.,1.)
    h,_ = skin.closest_points(p)
    outward,rim=indexed_feature_normals(skin,h)
    _,sg,_=inward_sheet_clearance(p,h,outward,rim)
    good = (np.minimum(bg,sg)>=.00021).reshape(len(v),len(grid))
    recipient_queries=[]
    for index,closed,label in recipient_surfaces:
        # The AABB is only a conservative query cull; every potentially near
        # point is tested on the actual complete fixed source triangles.
        low=index.mesh.vertices.min(0)-.00021;high=index.mesh.vertices.max(0)+.00021
        chosen=np.flatnonzero(np.all(p>=low,axis=1)&np.all(p<=high,axis=1))
        if len(chosen):
            hit,_=index.closest_points(p[chosen]);clear=hit['distance_m']>=.00021
            if closed:
                tensor=library().core.Tensor(((p[chosen]-index.origin)/index.scale).astype(np.float32))
                clear&=index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()<.5
            good.ravel()[chosen]&=clear
        recipient_queries.append(dict(partId=label,closedSolid=bool(closed),sampleQueries=len(chosen)))
        if time.monotonic()-started>=wall_seconds:
            return v.copy(),dict(complete=False,status='NativeDeadline',elapsedSeconds=time.monotonic()-started)
    states=[];centres=[]
    for k in range(n):
        raw = np.flatnonzero(owner==k);centres.append(v[raw].mean(0))
        possible = np.flatnonzero(good[raw].all(0))
        shifts=grid+guide_translation[raw[0]]
        possible=possible[np.linalg.norm(shifts[possible],axis=1)<=.005+1e-15]
        if k in protected:possible=possible[np.all(shifts[possible]==0.,axis=1)]
        if span_context is not None and len(possible):
            internal=np.flatnonzero(np.all(owner[f]==k,axis=1))
            if len(internal):
                posed=v[f[internal]][None]+shifts[possible,None,None,:]
                possible=possible[~source_span_hits(posed,span_context).any(1)]
        if not len(possible):
            return v.copy(),dict(complete=False,status='NoSampledRigidGroupCorridor',node=k,
                                 sourceControls=raw.tolist(),gridStates=len(grid),
                                 elapsedSeconds=time.monotonic()-started)
        possible = possible[np.argsort(np.sum(shifts[possible]**2,axis=1))[:48]]
        states.append(shifts[possible])
    centres=np.array(centres);transitions={}
    for a,b in edges:
            if parent.get(int(b))!=a:a,b=b,a
            delta=states[b][None]-states[a][:,None]
            transition=5*np.sum(delta**2,axis=2)
            axis=centres[b]-centres[a];length=np.linalg.norm(axis)
            if length>1e-10:
                moved=axis+delta
                transition[(moved@axis <= .05*length*length)|
                           (np.linalg.norm(moved,axis=2)>2*length+.0005)]=np.inf
            # A tree edge carries all its original source triangles. Reject
            # state transitions which invert or flatten even one face.
            membership=owner[f];affected=np.flatnonzero(np.any(membership==a,axis=1)&np.any(membership==b,axis=1))
            if len(affected):
                original=v[f[affected]];normal=np.cross(original[:,1]-original[:,0],original[:,2]-original[:,0])
                area=np.linalg.norm(normal,axis=1);normal/=area[:,None]
                for face,points,face_normal,source_area in zip(affected,original,normal,area):
                    posed=points[None,None]+np.where((owner[f[face]]==b)[None,None,:,None],states[b][None,:,None,:],states[a][:,None,None,:])
                    signed=np.einsum('ijk,k->ij',np.cross(posed[:,:,1]-posed[:,:,0],posed[:,:,2]-posed[:,:,0]),face_normal)
                    transition[signed<.22*source_area]=np.inf
            if span_context is not None and len(affected):
                points=v[f[affected]]
                posed=points[None,None]+np.where((owner[f[affected]]==b)[None,None,:,:,None],states[b][None,:,None,None,:],states[a][:,None,None,None,:])
                transition[source_span_hits(posed,span_context).any(2)]=np.inf
            transitions[int(a),int(b)]=transition
            if time.monotonic()-started >= wall_seconds:
                return v.copy(),dict(complete=False,status='NativeDeadline',elapsedSeconds=time.monotonic()-started)
    bans=[set()for _ in range(n)];replays=[];q=v.copy()
    for repair in range(6):
        costs=[np.sum(s*s,axis=1)for s in states];choices={}
        for k,forbidden in enumerate(bans):
            if forbidden:costs[k][list(forbidden)]=np.inf
        for a in reversed(order):
            for b in neighbours[a]:
                if parent.get(b)!=a:continue
                proposal=transitions[a,b]+costs[b][None]
                pick=np.argmin(proposal,axis=1);costs[a]+=proposal[np.arange(len(pick)),pick];choices[a,b]=pick
        if not np.isfinite(costs[0]).any():
            return v.copy(),dict(complete=False,status='NoContinuousSampledRoute',sourceSelfReplays=replays,elapsedSeconds=time.monotonic()-started)
        picked={0:int(np.argmin(costs[0]))}
        for a in order[1:]:picked[a]=int(choices[parent[a],a][picked[parent[a]]])
        translation=np.array([states[k][picked[k]]for k in range(n)]);q=v+translation[owner]
        from geometry_tools.contact import self_contact_summary
        from geometry_tools.model import Mesh
        audit,_=self_contact_summary(Mesh(q,f,'diagnostic:source-route'),example_limit=1000)
        if not audit['complete_within_pair_budget']:raise ValueError('Incomplete source-route self replay')
        replays.append(audit['defect_pairs'])
        if not audit['defect_pairs']:break
        # Finite conservative state elimination; never drop a source edge,
        # triangle or branch. Failure establishes no optimality claim.
        for first,second in audit['examples']:
            groups=np.unique(owner[f[[first,second]]])
            movable=[int(k)for k in groups if np.linalg.norm(translation[k])>1e-12]
            if not movable:return v.copy(),dict(complete=False,status='PreexistingUnmovedSourceCrossing',sourceSelfReplays=replays,elapsedSeconds=time.monotonic()-started)
            k=max(movable,key=lambda j:np.linalg.norm(translation[j]));bans[k].add(picked[k])
        if time.monotonic()-started>=wall_seconds:
            return v.copy(),dict(complete=False,status='NativeDeadline',sourceSelfReplays=replays,elapsedSeconds=time.monotonic()-started)
    else:
        return v.copy(),dict(complete=False,status='SourceFoldUnresolved',sourceSelfReplays=replays,elapsedSeconds=time.monotonic()-started)
    if time.monotonic()-started>=wall_seconds:
        return v.copy(),dict(complete=False,status='NativeDeadline',elapsedSeconds=time.monotonic()-started)
    return q,dict(complete=True,status='SampledRouteProposed',sourceGroups=n,
                  originalSourceGroups=original_nodes,cyclesRetainedAsRigidClosures=True,
                  maximumGroupTranslation_m=float(np.linalg.norm(translation,axis=1).max()),
                  maximumStatesPerNode=48,gridStates=len(grid),allSourceFacesRetained=True,
                  triangleClearancePending=True,openSkinContainmentNotEstablished=True,
                  skinOpenRimStateWitnesses=int(rim.sum()),
                  protectedSourceGroups=len(protected),protectedSourceRawControls=len(protected_raw),
                  fixedRecipientSurfaces=recipient_queries,openRecipientSheetsAreUnsigned=True,
                  originalFaceProjectionGuard=.22,
                  sourceSelfReplays=replays,sourceSpanScreened=span_context is not None,
                  savedRouteGuideUsed=guide_vertices is not None,gridSpacing_m=grid_spacing_m,localRadius_m=local_radius_m,
                  elapsedSeconds=time.monotonic()-started)
