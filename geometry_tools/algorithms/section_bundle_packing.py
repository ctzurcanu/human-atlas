"""Bounded proper-pose packing of separate vessel sections in fixed tissue.

Caller-supplied pose candidates must already satisfy bone/cover constraints.
This is section evidence only: longitudinal mouths and 3D geometry still need
explicit transport, source replay and whole-surface contact verification.
"""
import numpy as np
import shapely


def pack_section_polygons(polygons, poses, obstacles, *, reference_pose,
                          tissue_margin_m=.00025, mutual_margin_m=.00010,
                          maximum_displacement_m=.008, beam_width=64,
                          allowed_masks=None, obstacle_groups=(),return_configurations=False,
                          diversity_m=0.,preserve_translation_homotopy=False,
                          extra_pose_costs=None,pair_centre_limits=(),centre_connector_barriers=(),
                          primary_diversity_m=0.,maximum_per_primary_cell=2,
                          spread_primary_candidates=False):
    poses = np.asarray(poses, float)
    reference_pose = np.asarray(reference_pose, float)
    if poses.ndim != 2 or poses.shape[1] != 3 or not np.isfinite(poses).all():
        raise ValueError('Finite X/Z/roll proper poses required')
    if reference_pose.shape not in [(3,),(len(polygons),3)] or not np.isfinite(reference_pose).all():
        raise ValueError('Finite reference X/Z/roll required')
    if not (1 <= beam_width <= 256) or min(tissue_margin_m, mutual_margin_m) < 0 or maximum_displacement_m <= 0:
        raise ValueError('Bounded beam and physical clearances required')
    if not np.isfinite(diversity_m)or diversity_m<0:raise ValueError('Finite nonnegative layout diversity scale required')
    if not np.isfinite(primary_diversity_m)or primary_diversity_m<0 or not isinstance(maximum_per_primary_cell,int)or maximum_per_primary_cell<1:raise ValueError('Finite nonnegative primary diversity scale and positive per-cell budget required')
    if extra_pose_costs is not None:
        extra_pose_costs=np.asarray(extra_pose_costs,float)
        if extra_pose_costs.shape!=(len(polygons),len(poses))or not np.isfinite(extra_pose_costs).all()or np.any(extra_pose_costs<0):raise ValueError('Finite nonnegative cost for every actual profile pose required')
    limits={}
    for first,second,distance in pair_centre_limits:
        if not(0<=first<second<len(polygons))or not np.isfinite(distance)or distance<=0:raise ValueError('Explicit profile pair indices and positive finite centre limit required')
        limits[(int(first),int(second))]=float(distance)
    connector_barriers={}
    for first,second,barrier in centre_connector_barriers:
        if not(0<=first<second<len(polygons))or not shapely.is_valid(barrier)or shapely.is_empty(barrier):raise ValueError('Explicit member indices and a valid nonempty connector barrier required')
        connector_barriers[(int(first),int(second))]=barrier

    def shape(points):
        p = np.asarray(points, float)
        if p.ndim != 2 or p.shape[1] != 2 or len(p) < 3 or not np.isfinite(p).all():
            raise ValueError('Literal finite closed source contours required')
        geometry = shapely.polygons(p)
        if not shapely.is_valid(geometry) or shapely.area(geometry) <= 0:
            raise ValueError('Invalid source polygon; no repair or closure inferred')
        return p, geometry

    forbidden = shapely.union_all([shape(p)[1] for p in obstacles])
    extra=[(shapely.union_all([shape(p)[1]for p in group]),margin)for group,margin in obstacle_groups]
    members = []
    c, s = np.cos(poses[:, 2]), np.sin(poses[:, 2])
    for member_index,points in enumerate(polygons):
        ref=reference_pose if reference_pose.ndim==1 else reference_pose[member_index]
        cr, sr = np.cos(ref[2]), np.sin(ref[2])
        p, original = shape(points)
        placed = np.stack([c[:, None]*p[None, :, 0]+s[:, None]*p[None, :, 1],
                           -s[:, None]*p[None, :, 0]+c[:, None]*p[None, :, 1]], axis=2) + poses[:, None, :2]
        reference = p @ np.array([[cr, sr], [-sr, cr]]).T + ref[:2]
        movement = np.linalg.norm(placed-reference, axis=2)
        geometries = shapely.polygons(placed)
        valid = (movement.max(axis=1) <= maximum_displacement_m)
        if not shapely.is_empty(forbidden):
            valid &= shapely.disjoint(geometries, forbidden)
            valid &= shapely.distance(geometries, forbidden) >= tissue_margin_m
        if allowed_masks is not None:
            mask=np.asarray(allowed_masks[member_index],bool)
            if mask.shape!=(len(poses),):raise ValueError('One explicit candidate mask per actual profile required')
            valid &= mask
        for geometry,margin in extra:
            if not shapely.is_empty(geometry):valid &= shapely.disjoint(geometries,geometry)&(shapely.distance(geometries,geometry)>=margin)
        indices = np.flatnonzero(valid)
        cost = np.mean(movement**2, axis=1)
        if extra_pose_costs is not None:cost+=extra_pose_costs[member_index]
        indices = indices[np.argsort(cost[indices], kind='stable')]
        if not len(indices):
            return dict(feasible=False, reason='No individual proper pose', member=len(members), complete=False)
        members.append(dict(source=p, source_area=float(shapely.area(original)), geometry=geometries,
                            placed=placed, indices=indices, cost=cost))
    if not members:
        raise ValueError('At least one actual section required')
    beam = [(0., ())]
    for index, member in enumerate(members):
        expanded = []
        indices = member['indices']
        for score, chosen in beam:
            valid = np.ones(len(indices), bool)
            for prior_member_index,(prior, prior_index) in enumerate(zip(members[:index], chosen)):
                first = prior['geometry'][prior_index]
                valid &= shapely.disjoint(member['geometry'][indices], first)
                valid &= shapely.distance(member['geometry'][indices], first) >= mutual_margin_m
                if(prior_member_index,index)in limits:
                    valid&=np.linalg.norm(member['placed'][indices].mean(1)-prior['placed'][prior_index].mean(0),axis=1)<=limits[(prior_member_index,index)]
                if(prior_member_index,index)in connector_barriers:
                    ends=member['placed'][indices].mean(1);begin=np.broadcast_to(prior['placed'][prior_index].mean(0),ends.shape);connectors=shapely.linestrings(np.stack([begin,ends],axis=1));valid&=shapely.disjoint(connectors,connector_barriers[(prior_member_index,index)])
                if preserve_translation_homotopy:
                    if np.any(abs(poses[indices,2])>1e-12)or abs(poses[prior_index,2])>1e-12 or np.any(abs(reference_pose[...,2])>1e-12):raise ValueError('Translation homotopy requires explicit zero-roll poses')
                    forbidden_pair=shapely.convex_hull(shapely.multipoints((prior['source'][:,None]-member['source'][None]).reshape(-1,2)))
                    reference_first=reference_pose if reference_pose.ndim==1 else reference_pose[prior_member_index]
                    reference_second=reference_pose if reference_pose.ndim==1 else reference_pose[index]
                    beginning=reference_second[:2]-reference_first[:2];endings=poses[indices,:2]-poses[prior_index,:2];segments=shapely.linestrings(np.stack([np.broadcast_to(beginning,endings.shape),endings],axis=1));valid &= shapely.disjoint(segments,shapely.buffer(forbidden_pair,mutual_margin_m))
            possible = indices[valid]
            if index==0 and spread_primary_candidates and len(possible)>beam_width:
                # Reserve actual spatial alternatives before the first beam
                # truncation. Later diversity cannot recover discarded poses.
                centres=member['placed'][possible].mean(1);selected=[0]
                distance=np.sum((centres-centres[0])**2,axis=1);distance[0]=-1
                for _ in range(1,beam_width):
                    pick=int(np.argmax(distance));selected.append(pick)
                    distance=np.minimum(distance,np.sum((centres-centres[pick])**2,axis=1));distance[selected]=-1
                possible=possible[selected]
            else:possible=possible[:beam_width]
            expanded.extend((score+float(member['cost'][j]), chosen+(int(j),)) for j in possible)
        if not expanded:
            return dict(feasible=False, reason='Bounded beam found no mutually clear packing', member=index, complete=False,
                        individualFeasiblePoseCounts=[len(m['indices']) for m in members])
        ordered=sorted(expanded,key=lambda row:(row[0],row[1]))
        if diversity_m or primary_diversity_m:
            seen=set();beam=[];primary_counts={}
            for entry in ordered:
                key=None if not diversity_m else tuple(np.rint(np.concatenate([m['placed'][i].mean(0)for m,i in zip(members,entry[1])])/diversity_m).astype(int))
                primary=None if not primary_diversity_m else tuple(np.rint(members[0]['placed'][entry[1][0]].mean(0)/primary_diversity_m).astype(int))
                if key is not None and key in seen:continue
                if primary is not None and primary_counts.get(primary,0)>=maximum_per_primary_cell:continue
                if key is not None:seen.add(key)
                if primary is not None:primary_counts[primary]=primary_counts.get(primary,0)+1
                beam.append(entry)
                if len(beam)==beam_width:break
        else:beam=ordered[:beam_width]
    score, chosen = beam[0]
    selected = [m['placed'][i] for m, i in zip(members, chosen)]
    geometries = [shapely.polygons(p) for p in selected]
    result=dict(feasible=True, complete=True, score_m2=score, selectedPoseIndices=list(chosen),
                selectedPoses=poses[list(chosen)].tolist(), polygons=[p.tolist() for p in selected],
                individualFeasiblePoseCounts=[len(m['indices']) for m in members],
                minimumTissueClearance_m=None if shapely.is_empty(forbidden)else float(min(shapely.distance(p, forbidden) for p in geometries)),
                minimumMutualClearance_m=None if len(geometries)<2 else float(min(shapely.distance(a,b) for j,a in enumerate(geometries) for b in geometries[j+1:])),
                maximumAreaRelativeError=float(max(abs(shapely.area(g)/m['source_area']-1) for g,m in zip(geometries,members))),
                beamWidth=beam_width, allTransformsProperUnitScale=True,translationHomotopyChecked=preserve_translation_homotopy,threeDimensionalAcceptance=False,
                explicitPairCentreLimits=[list(x)for x in pair_centre_limits],extraPoseCostsApplied=extra_pose_costs is not None,centreConnectorBarrierPairs=[list(k)for k in connector_barriers],centreConnectorsAreHeuristicsNotConduits=bool(connector_barriers),primaryDiversity_m=primary_diversity_m,maximumPerPrimaryCell=maximum_per_primary_cell,primaryCandidatesSpatiallySpread=bool(spread_primary_candidates))
    if return_configurations:
        result['configurations']=[dict(score_m2=s,poseIndices=list(indices),poses=poses[list(indices)].tolist(),centresXZ_m=[m['placed'][i].mean(0).tolist()for m,i in zip(members,indices)])for s,indices in beam]
    return result
