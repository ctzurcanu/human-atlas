"""Bounded layered route search around a curved guide, without source surgery.

Dynamic programming keeps the full transverse lattice at each station. It does
not contract source-incidence cycles into one rigid body. Clearance samples
and optional segment tests are proposal screens; emitted triangles still need
independent full contact/containment, source-calibre and anatomical review.
"""
import time
import numpy as np
from .tube_sections import arclength,rotation_minimizing_frames


def route_station_correspondence(guide,route,*,grid_m,half_width_m):
    """Recover a literal saved lattice route's source-station correspondence.

    Normalizing the route's total length loses this correspondence when a
    detour lengthens one part of the route. Last visits are ordered because
    lattice edges advance by at most one guide station. No geometry is moved.
    """
    from scipy.spatial import cKDTree
    g=np.asarray(guide,float);q=np.asarray(route,float)
    if g.ndim!=2 or g.shape[1]!=3 or q.ndim!=2 or q.shape[1]!=3 or len(q)<3:
        raise ValueError('Finite saved guide and route required')
    n=int(np.floor(half_width_m/grid_m));width=2*n+1
    if len(g)*width*width>120000:raise ValueError('Saved route exceeds finite reconstruction budget')
    offsets=(np.indices((width,width)).reshape(2,-1).T-n)*grid_m
    frame=rotation_minimizing_frames(g)
    lattice=g[:,None]+np.einsum('sij,kj->ski',frame[:,:,:2],offsets)
    residual,ids=cKDTree(lattice.reshape(-1,3)).query(q);station=ids//(width*width)
    if (residual>1e-12).any()or station[0]!=0 or station[-1]!=len(g)-1 or np.any(abs(np.diff(station))>1):
        raise ValueError('Saved points are not a literal adjacent-station lattice route')
    s=arclength(q);last=[]
    for k in range(len(g)):
        found=np.flatnonzero(station==k)
        if not len(found):raise ValueError('Saved route skips a source guide station')
        last.append(s[found[-1]])
    mapped=np.array(last);mapped[0]=0.
    if np.any(np.diff(mapped)<=0):raise ValueError('Saved route has no strictly ordered station correspondence')
    return mapped,dict(literalLatticeMaximumResidual_m=float(residual.max()),
        sourceStations=len(g),routePoints=len(q),backwardStationSteps=int((np.diff(station)<0).sum()),
        sourceStationLastVisitsOrdered=True,sourceGeometryChanged=False)


def layered_tube_corridor(guide,radii,clearance,*,grid_m=.0005,half_width_m=.004,
                          segment_clear=None,fixed_prefix=2,fixed_distal=True,wall_seconds=16.,
                          maximum_points=120000,transverse_neighbour_cells=1,allow_lateral_detours=False,
                          configuration_clearance=None):
    started=time.monotonic();g=np.asarray(guide,float);r=np.asarray(radii,float)
    if (g.ndim!=2 or g.shape[1]!=3 or len(g)<3 or r.shape!=(len(g),) or
        not np.isfinite(g).all() or not np.isfinite(r).all() or (r<=0).any() or
        not 0<grid_m<=.001 or not grid_m<=half_width_m<=.01 or
        not 0<wall_seconds<=20 or not 0<=fixed_prefix<len(g)):
        raise ValueError('Finite guide, source radii, explicit lattice and native deadline required')
    if not isinstance(transverse_neighbour_cells,int)or not 1<=transverse_neighbour_cells<=4:
        raise ValueError('Explicit transverse transition stencil within one to four cells required')
    def deadline():
        if time.monotonic()-started>wall_seconds:raise TimeoutError('Native layered tube-route deadline')
    def available(points,stations):
        if configuration_clearance is None:return np.asarray(clearance(points),float)
        # The callback checks actual immutable section outlines. Add the
        # unmodified radius only to retain the existing surplus-cost units;
        # the feasibility comparison cancels it exactly. No inscribed-radius
        # certificate is implied by this oriented configuration proposal.
        station=np.broadcast_to(np.asarray(stations,float),(len(points),))
        radius=np.interp(station,np.arange(len(r)),r)
        margin=np.asarray(configuration_clearance(points,station),float)
        if margin.shape!=(len(points),):raise ValueError('One configuration margin per state required')
        return margin+radius
    frames=rotation_minimizing_frames(g);n=int(np.floor(half_width_m/grid_m));width=2*n+1
    ij=np.indices((width,width)).reshape(2,-1).T-n;offsets=ij*grid_m;zero=int(np.flatnonzero(np.all(ij==0,axis=1))[0])
    q=g[:,None]+np.einsum('sij,kj->ski',frames[:,:,:2],offsets)
    if q.shape[0]*q.shape[1]>maximum_points:raise ValueError('Layered route exceeds its explicit point budget')
    gap=available(q.reshape(-1,3),np.repeat(np.arange(len(g)),len(ij))).reshape(q.shape[:2]);deadline()
    good=np.isfinite(gap)&(gap>=r[:,None])
    if not allow_lateral_detours:
        for k in range(fixed_prefix):good[k,np.arange(len(ij))!=zero]=False
        if fixed_distal:good[-1,np.arange(len(ij))!=zero]=False
    missing=np.flatnonzero(~good.any(1))
    if len(missing):
        return None,dict(complete=False,status='NoSampledStationCorridor',stations=missing.tolist(),
                         maximumGap_m=np.max(np.where(np.isfinite(gap),gap,-np.inf),axis=1)[missing].tolist(),
                         requiredRadius_m=r[missing].tolist(),elapsedSeconds=time.monotonic()-started)
    if allow_lateral_detours:
        # A forward-only station DAG can reject a real corridor whose route
        # first travels sideways around a tissue sheet. Use the SAME bounded
        # lattice and actual context tests, with reversible lateral edges.
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import dijkstra
        ids=np.full(good.shape,-1,int);ids[good]=np.arange(int(good.sum()))
        first=zero
        last=zero
        first_station=max(0,fixed_prefix-1)
        if fixed_prefix and not np.all(good[:fixed_prefix,zero]):
            raise ValueError('An explicit protected source-prefix station is blocked')
        if first_station:
            prefix_a,prefix_b=g[:first_station],g[1:first_station+1]
            mid=available((prefix_a+prefix_b)*.5,np.arange(first_station)+.5)
            safe=np.isfinite(mid)&(mid>=(r[:first_station]+r[1:first_station+1])*.5)
            if segment_clear is not None:safe&=np.asarray(segment_clear(prefix_a,prefix_b),bool)
            if not safe.all():
                return None,dict(complete=False,status='BlockedProtectedSourcePrefix',
                    elapsedSeconds=time.monotonic()-started)
        start_ids=np.array([ids[first_station,first]])if fixed_prefix else ids[0,good[0]]
        end_ids=np.array([ids[-1,last]])if fixed_distal else ids[-1,good[-1]]
        if np.any(start_ids<0)or np.any(end_ids<0):raise ValueError('Explicit graph mouth is blocked')
        cube=ids.reshape(len(g),width,width);edges_a=[];edges_b=[];edge_cost=[];queries=0
        live=q[good];surplus=gap[good]-np.broadcast_to(r[:,None],good.shape)[good]
        station_grid=np.broadcast_to(np.arange(len(g))[:,None],good.shape)[good]
        # Thirteen undirected half-neighbours form a full 26-neighbour graph.
        for ds in [0,1]:
            for di in [-1,0,1]:
                for dj in [-1,0,1]:
                    if ds==0 and (di<0 or di==0 and dj<=0):continue
                    shift=(ds,di,dj);aa=[];bb=[]
                    for size,d in zip(cube.shape,shift):
                        aa.append(slice(max(0,-d),min(size,size-d)))
                        bb.append(slice(max(0,d),min(size,size+d)))
                    a,b=cube[tuple(aa)].ravel(),cube[tuple(bb)].ravel();ok=(a>=0)&(b>=0);a,b=a[ok],b[ok]
                    keep=np.ones(len(a),bool)
                    for low in range(0,len(a),16000):
                        x,y=live[a[low:low+16000]],live[b[low:low+16000]]
                        station=(station_grid[a[low:low+16000]]+station_grid[b[low:low+16000]])*.5
                        mg=available((x+y)*.5,station)
                        radius=np.interp(station,np.arange(len(r)),r)
                        valid=np.isfinite(mg)&(mg>=radius)
                        if segment_clear is not None:valid&=np.asarray(segment_clear(x,y),bool)
                        keep[low:low+16000]=valid;queries+=len(x);deadline()
                    a,b=a[keep],b[keep]
                    cost=np.linalg.norm(live[a]-live[b],axis=1)*(1+(.0005/np.maximum((surplus[a]+surplus[b])*.5+.0001,.0001))**2)
                    edges_a.append(a);edges_b.append(b);edge_cost.append(cost)
        a,b,c=np.concatenate(edges_a),np.concatenate(edges_b),np.concatenate(edge_cost)
        graph=coo_matrix((np.r_[c,c],(np.r_[a,b],np.r_[b,a])),shape=(len(live),len(live))).tocsr()
        if len(start_ids)==1:
            distance,parent=dijkstra(graph,directed=False,indices=int(start_ids[0]),return_predecessors=True)
            sources=np.full(len(live),int(start_ids[0]),int)
        else:
            # Free mouths are a set of eligible attachment states. Choosing
            # their nearest state before routing can falsely trap the route
            # in the donor mouth's disconnected component.
            distance,parent,sources=dijkstra(graph,directed=False,indices=start_ids,
                return_predecessors=True,min_only=True)
        end_cost=distance[end_ids]+np.linalg.norm(live[end_ids]-g[-1],axis=1)
        end_id=int(end_ids[int(np.argmin(end_cost))]);deadline()
        if not np.isfinite(distance[end_id]):
            from scipy.sparse.csgraph import connected_components
            from scipy.spatial import cKDTree
            _,labels=connected_components(graph,directed=False)
            source_ids=np.flatnonzero(np.isfinite(distance));target_ids=np.flatnonzero(labels==labels[end_id])
            near,other=cKDTree(live[target_ids]).query(live[source_ids])
            k=int(np.argmin(near));pair=live[[source_ids[k],target_ids[int(other[k])]]]
            return None,dict(complete=False,status='DisconnectedSampledClearanceGraph',gridPoints=int(good.size),
                gridSpacing_m=grid_m,halfWidth_m=half_width_m,
                actualSectionOutlineProposal=configuration_clearance is not None,
                freeStates=len(live),screenedSegments=queries,
                sourceRegionStates=len(source_ids),targetRegionStates=len(target_ids),
                sourceRegionBounds_m=np.array([live[source_ids].min(0),live[source_ids].max(0)]).tolist(),
                targetRegionBounds_m=np.array([live[target_ids].min(0),live[target_ids].max(0)]).tolist(),
                closestRegionPoints_m=pair.tolist(),closestRegionDistance_m=float(near[k]),
                elapsedSeconds=time.monotonic()-started)
        start_id=int(sources[end_id]);selected=[end_id]
        while selected[-1]!=start_id:
            k=int(parent[selected[-1]])
            if k<0 or k in selected:raise ValueError('Incomplete graph route predecessor chain')
            selected.append(k)
        selected=np.array(selected[::-1]);route=np.r_[g[:first_station],live[selected]]
        return route,dict(complete=True,status='FiniteSampledGraphRoute',stations=len(g),routePoints=len(route),
            gridPoints=int(good.size),freeStates=len(live),gridSpacing_m=grid_m,halfWidth_m=half_width_m,
            screenedSegments=queries,minimumNodeRadiusSurplus_m=float(surplus[selected].min()),
            sourceRadiiReduced=False,fullTriangleAuditRequired=True,continuousInscribedRadiusProven=False,
            fixedPrefixStations=fixed_prefix,fixedDistal=fixed_distal,
            proposedEndpointDisplacements_m=np.linalg.norm(route[[0,-1]]-g[[0,-1]],axis=1).tolist(),
            eligibleStartStates=len(start_ids),eligibleEndStates=len(end_ids),
            actualSectionOutlineProposal=configuration_clearance is not None,
            lateralDetoursAllowed=True,anatomicalAcceptance=False,nativeWallSeconds=wall_seconds,
            elapsedSeconds=time.monotonic()-started)
    # Each state can arrive from the same or one neighbouring transverse cell.
    # World-space edge lengths and actual segment screens account for guide
    # curvature and rotation of its transported cross-section frames.
    reach=range(-transverse_neighbour_cells,transverse_neighbour_cells+1)
    steps=np.array([[a,b]for a in reach for b in reach],int)
    predecessor_ij=ij[:,None]+steps[None]
    valid=np.all(abs(predecessor_ij)<=n,axis=2)
    predecessor=np.clip((predecessor_ij[:,:,0]+n)*width+predecessor_ij[:,:,1]+n,0,len(ij)-1)
    costs=np.where(good[0],np.linalg.norm(offsets,axis=1),np.inf)
    history=np.full((len(g),len(ij)),-1,int)
    segment_queries=0
    for station in range(1,len(g)):
        previous_costs=costs.copy()
        a=q[station-1,predecessor];b=np.broadcast_to(q[station,:,None],a.shape)
        live=valid&good[station,:,None]&good[station-1,predecessor]&np.isfinite(costs[predecessor])
        candidate_segments=int(live.sum());midpoint_rejects=segment_rejects=0
        if live.any():
            # Include interior points; node-only screening can tunnel across
            # an arbitrarily thin fixed tissue sheet between lattice planes.
            mid=(a[live]+b[live])*.5;mg=available(mid,np.full(len(mid),station-.5));deadline()
            mid_radius=(r[station-1]+r[station])*.5
            ok=np.isfinite(mg)&(mg>=mid_radius)
            midpoint_rejects=int((~ok).sum())
            if segment_clear is not None:
                segment_ok=np.asarray(segment_clear(a[live],b[live]),bool)
                segment_rejects=int((ok&~segment_ok).sum());ok&=segment_ok
            segment_queries+=len(mid);live[live]=ok
        length=np.linalg.norm(b-a,axis=2)
        proximity=(.0005/np.maximum(gap[station]-r[station]+.0001,.0001))**2
        preference=.08*np.sum(offsets**2,axis=1)/(half_width_m**2)
        energy=costs[predecessor]+length*(1+proximity[:,None]+preference[:,None])
        energy[~live]=np.inf;choice=energy.argmin(1);costs=energy[np.arange(len(ij)),choice]
        history[station]=predecessor[np.arange(len(ij)),choice]
        if not np.isfinite(costs).any():
            picked=int(np.argmin(previous_costs));partial=np.empty(station,int);partial[-1]=picked
            for k in range(station-1,0,-1):partial[k-1]=history[k,partial[k]]
            return None,dict(complete=False,status='NoContinuousSampledCorridor',station=station,
                             nodeStates=int(good[station].sum()),segmentQueries=segment_queries,
                             candidateSegments=candidate_segments,midpointRejected=midpoint_rejects,segmentRejected=segment_rejects,
                             partialRoute_m=q[np.arange(station),partial].tolist(),
                             nextFreeStates_m=q[station,good[station]].tolist(),
                             elapsedSeconds=time.monotonic()-started)
        deadline()
    selected=np.empty(len(g),int);selected[-1]=zero if fixed_distal else int(np.argmin(costs+np.linalg.norm(offsets,axis=1)))
    if not np.isfinite(costs[selected[-1]]):raise ValueError('Distal route state is unreachable')
    for k in range(len(g)-1,0,-1):selected[k-1]=history[k,selected[k]]
    route=q[np.arange(len(g)),selected]
    if not np.array_equal(route[:fixed_prefix],g[:fixed_prefix]) or fixed_distal and not np.array_equal(route[-1],g[-1]):
        raise ValueError('An explicit fixed route station moved')
    return route,dict(complete=True,status='FiniteSampledRoute',stations=len(g),statesPerStation=len(ij),
        gridPoints=int(np.prod(gap.shape)),gridSpacing_m=grid_m,halfWidth_m=half_width_m,
        transverseNeighbourCells=transverse_neighbour_cells,
        segmentQueries=segment_queries,segmentInteriorsScreened=True,
        minimumNodeRadiusSurplus_m=float(np.min(gap[np.arange(len(g)),selected]-r)),
        sourceRadiiReduced=False,fullTriangleAuditRequired=True,continuousInscribedRadiusProven=False,
        fixedPrefixStations=fixed_prefix,fixedDistal=fixed_distal,
        proposedEndpointDisplacements_m=np.linalg.norm(route[[0,-1]]-g[[0,-1]],axis=1).tolist(),
        actualSectionOutlineProposal=configuration_clearance is not None,
        anatomicalAcceptance=False,nativeWallSeconds=wall_seconds,elapsedSeconds=time.monotonic()-started)
