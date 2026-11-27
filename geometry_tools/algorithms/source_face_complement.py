"""Complement contact fragments on their actual original source triangles.

An occupied body can share the same exterior plane as its host. Subtracting
the host then loses that zero-thickness common exterior. Recovering free
boundary face-by-face preserves those real source patches without inventing
a cap. Coordinates and parent charts are explicit; all callers must run full
topology, self, reciprocal-interface and occupied-shape verification afterwards.
"""
from collections import defaultdict
import numpy as np
import manifold3d as md
from scipy.spatial import cKDTree
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank


def complement_source_faces(positions, source_positions, source_faces,
                            contact_faces, contact_parents, *, tolerance=2e-12,
                            planar_scale=1e6, backend='clipper'):
    v=np.asarray(positions,float);sv=np.asarray(source_positions,float);sf=np.asarray(source_faces,int)
    cf=np.asarray(contact_faces,int);parents=np.asarray(contact_parents,int)
    if len(cf)!=len(parents)or np.any(parents<0)or np.any(parents>=len(sf)):
        raise ValueError('Every contact fragment needs its actual source face')
    bank=CommonPositionBank(v,tolerance);source_map=np.array([bank.insert(p)for p in sv],int)
    groups=defaultdict(list)
    for i,parent in enumerate(parents):groups[int(parent)].append(i)
    free=[];free_parent=[];maximum_plane=0.;maximum_snap=0.;covered=0;zero_fragments=0;collinear_controls=[];vanished_spurs=[]
    # Clipper's internal fixed-point rounding is made sub-picometre in metres
    # by performing each local 2D operation in micrometres. No source point is
    # rounded onto a clinical or anatomical grid.
    scale=float(planar_scale)
    if not np.isfinite(scale) or scale<=0:raise ValueError('Positive finite planar scale required')
    if backend not in ['clipper','geos','triangle']:raise ValueError('Unknown planar complement backend')
    for parent,tri_ids in enumerate(sf):
        if parent not in groups:
            free.append(source_map[tri_ids]);free_parent.append(parent);continue
        p=sv[tri_ids];e=p[1]-p[0];normal=np.cross(e,p[2]-p[0]);length=np.linalg.norm(normal)
        if length==0:raise ValueError('A source boundary triangle has zero area')
        n=normal/length;u=e/np.linalg.norm(e);w=np.cross(n,u);basis=np.stack([u,w],axis=1)
        original=(p-p[0])@basis*scale
        fragments=cf[np.asarray(groups[parent],int)];points=v[fragments]
        residual=float(np.abs((points-p[0])@n).max());maximum_plane=max(maximum_plane,residual)
        arithmetic=128*np.finfo(float).eps*max(1.,float(np.abs(p).max()))
        if residual>tolerance+arithmetic:raise ValueError('A contact fragment left its source parent plane')
        xy=(points-p[0])@basis*scale
        a=xy[:,1]-xy[:,0];b=xy[:,2]-xy[:,0];signed=a[:,0]*b[:,1]-a[:,1]*b[:,0]
        xy[signed<0]=xy[signed<0][:,[0,2,1]]
        if backend=='triangle':
            import triangle as tr
            from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature
            source_ids=np.unique(np.r_[source_map[tri_ids],fragments.ravel()])
            source_xy=(np.asarray([bank.rows[int(i)]for i in source_ids])-p[0])@basis*scale
            lookup={int(i):k for k,i in enumerate(source_ids)}
            source_edges=list(zip(source_map[tri_ids],np.roll(source_map[tri_ids],-1)))
            edges=source_edges+list(boundary_signature(fragments))
            segments=np.array(sorted({tuple(sorted((lookup[int(a)],lookup[int(b)]))) for a,b in edges}),int)
            triangulation=tr.triangulate(dict(vertices=source_xy,segments=segments),'pQ')
            generated=np.asarray(triangulation.get('triangles',[]),int).reshape(-1,3)
            coordinates=triangulation['vertices'];centers=coordinates[generated].mean(1)
            kept=np.zeros(len(generated),bool)
            for cut in xy:
                aa,bb=cut[1]-cut[0],cut[2]-cut[0];det=aa[0]*bb[1]-aa[1]*bb[0]
                if det==0:continue
                delta=centers-cut[0]
                u=(delta[:,0]*bb[1]-delta[:,1]*bb[0])/det
                w=(aa[0]*delta[:,1]-aa[1]*delta[:,0])/det
                kept|=(u>=0)&(w>=0)&(u+w<=1)
            tree=cKDTree(source_xy);distance,index=tree.query(coordinates)
            mapped=[]
            for coord,d,i in zip(coordinates,distance,index):
                if d<=tolerance*scale:
                    mapped.append(int(source_ids[int(i)]));maximum_snap=max(maximum_snap,float(d/scale))
                else:mapped.append(bank.insert(p[0]+(coord@basis.T)/scale))
            for tri in np.asarray(mapped,int)[generated[~kept]]:
                pp=np.asarray([bank.rows[int(i)]for i in tri]);nn=np.cross(pp[1]-pp[0],pp[2]-pp[0])
                if np.all(nn==0):zero_fragments+=1;continue
                if np.dot(nn,normal)<0:tri=tri[[0,2,1]]
                free.append(tri);free_parent.append(parent)
            continue
        if backend=='geos':
            # Floating precision overlay plus constrained triangulation avoids
            # the automatic near-collinear simplification of integer Clipper.
            # No buffer/repair or unconstrained triangle filtering is used.
            from shapely import Polygon,union_all,difference,constrained_delaunay_triangles
            whole=Polygon(original);mask=union_all([Polygon(tri) for tri in xy],grid_size=0.)
            region=difference(whole,mask,grid_size=0.)
            if not region.is_valid:raise ValueError('Invalid floating source-face complement')
            if region.is_empty:covered+=1;continue
            source_ids=np.unique(np.r_[source_map[tri_ids],fragments.ravel()])
            source_xy=(np.asarray([bank.rows[int(i)]for i in source_ids])-p[0])@basis*scale
            tree=cKDTree(source_xy)
            for polygon in constrained_delaunay_triangles(region).geoms:
                flat=np.asarray(polygon.exterior.coords)[:-1]
                if len(flat)!=3:raise ValueError('Constrained triangulator emitted a non-triangle')
                distance,index=tree.query(flat);ids=[]
                for coord,d,i in zip(flat,distance,index):
                    if d<=tolerance*scale:
                        ids.append(int(source_ids[int(i)]));maximum_snap=max(maximum_snap,float(d/scale))
                    else:ids.append(bank.insert(p[0]+(coord@basis.T)/scale))
                tri=np.asarray(ids,int);pp=np.asarray([bank.rows[int(i)]for i in tri]);nn=np.cross(pp[1]-pp[0],pp[2]-pp[0])
                if np.all(nn==0):zero_fragments+=1;continue
                if np.dot(nn,normal)<0:tri=tri[[0,2,1]]
                free.append(tri);free_parent.append(parent)
            continue
        whole=md.CrossSection([original]);mask=md.CrossSection(list(xy));region=whole-mask
        polygons=region.to_polygons()
        if not polygons:covered+=1;continue
        source_ids=np.unique(np.r_[source_map[tri_ids],fragments.ravel()]);source_xy=(np.asarray([bank.rows[int(i)]for i in source_ids])-p[0])@basis*scale
        tree=cKDTree(source_xy);flat=np.concatenate(polygons);distance,index=tree.query(flat)
        ids=[]
        for q,d,i in zip(flat,distance,index):
            if d<=tolerance*scale:
                ids.append(int(source_ids[int(i)]));maximum_snap=max(maximum_snap,float(d/scale))
            else:ids.append(bank.insert(p[0]+(q@basis.T)/scale))
        # Adjacent Clipper controls can map to the same existing source point.
        # Triangulate the actual canonical contours, rather than triangulating
        # a tiny polygon first and only then identifying its corners.
        clean=[];clean_ids=[];offset=0
        for polygon in polygons:
            ring=ids[offset:offset+len(polygon)];offset+=len(polygon)
            ring=[i for k,i in enumerate(ring)if k==0 or i!=ring[k-1]]
            if len(ring)>1 and ring[0]==ring[-1]:ring.pop()
            # Clipper can leave a backtracking, zero-width spur on a source
            # edge. Remove only its numerically collinear contour visit, never
            # a source position. Shared-edge conformity can reinsert any
            # necessary boundary control on the surviving face afterwards.
            changed=True
            while changed and len(ring)>3:
                changed=False
                for k,control in enumerate(ring):
                    a=np.asarray(bank.rows[ring[k-1]]);b=np.asarray(bank.rows[control]);c=np.asarray(bank.rows[ring[(k+1)%len(ring)]])
                    length=np.linalg.norm(c-a)
                    residual=np.linalg.norm(np.cross(b-a,c-a))/length if length else np.linalg.norm(b-a)
                    if residual<=tolerance:
                        collinear_controls.append(dict(parentFace=parent,control=int(control),lineResidualMeters=float(residual)))
                        ring.pop(k);changed=True;break
            if len(ring)<3:continue
            clean.append((np.asarray([bank.rows[int(i)]for i in ring])-p[0])@basis*scale)
            clean_ids.extend(ring)
        if not clean:continue
        local=np.asarray(md.triangulate(clean,epsilon=1e-8,allow_convex=False),int)
        mapped=np.asarray(clean_ids,int)[local]
        for tri in mapped:
            pp=np.asarray([bank.rows[int(i)]for i in tri])
            nn=np.cross(pp[1]-pp[0],pp[2]-pp[0])
            if np.all(nn==0):zero_fragments+=1;continue
            if np.dot(nn,normal)<=0:
                longest=max(np.linalg.norm(pp[1]-pp[0]),np.linalg.norm(pp[2]-pp[1]),np.linalg.norm(pp[0]-pp[2]))
                altitude=np.linalg.norm(nn)/longest
                if altitude<=tolerance:
                    vanished_spurs.append(dict(parentFace=parent,controls=tri.tolist(),altitudeMeters=float(altitude),areaMetersSquared=float(np.linalg.norm(nn)*.5)))
                    continue
                raise ValueError(f'A free source fragment reversed on actual source face {parent}, controls {tri.tolist()}')
            free.append(tri);free_parent.append(parent)
    return bank.positions,np.asarray(free,int).reshape(-1,3),np.asarray(free_parent,int),dict(
        sourceFaces=len(sf),sourceFacesWithContacts=len(groups),sourceFacesCompletelyCovered=covered,
        freeFragments=len(free),maximumContactParentPlaneResidualMeters=maximum_plane,
        maximumExistingControlLookupResidualMeters=maximum_snap,
        maximumCommonBankCanonicalizationMeters=bank.maximum_snap,
        originalPositionPrefixLiteral=np.array_equal(bank.positions[:len(v)],v),
        constructionToleranceMeters=tolerance,localPlanarScale=scale,exactZeroAreaFragmentsRemoved=zero_fragments,
        planarBackend=backend,
        numericalCollinearContourVisitsRemoved=collinear_controls,
        numericalBacktrackingSpurFragmentsRemoved=vanished_spurs,
        commonExteriorSourcePatchesRetained=True,anatomicalCapsAdded=False)


def conform_surface_boundary_edges(positions,surfaces,cells,*,tolerance=2e-12,maximum_rounds=8,record_face_parents=False,allow_collinear_cycles=False):
    """Split hanging boundary edges at existing incident boundary controls.

    One canonical surface changes, so its reciprocal owner changes identically.
    This inserts no coordinate and never fills a positive-area boundary loop.
    """
    v=np.asarray(positions,float);out={sid:np.asarray(f,int).copy()for sid,f in surfaces.items()};edits=[];cycle_proofs=[]
    ancestry={sid:np.arange(len(f),dtype=int)for sid,f in out.items()}if record_face_parents else None
    for iteration in range(maximum_rounds):
        pending={}
        for cell,uses in cells.items():
            fs=[];owners=[]
            for use in uses:
                f=out[use['surface']];fs.append(f if use['side']==1 else f[:,[0,2,1]])
                owners.extend((use['surface'],i)for i in range(len(f)))
            faces=np.concatenate(fs);edges=np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]
            keys,inverse,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True);boundary=keys[counts==1]
            if not len(boundary):continue
            controls=np.unique(boundary);tree=cKDTree(v[controls]);edge_set={tuple(q)for q in boundary}
            adjacency=defaultdict(list)
            for a,b in edges[counts[inverse]==1]:adjacency[int(a)].append(int(b))
            cycle_controls={}
            if allow_collinear_cycles:
                # Opposite chains can have interleaved controls (A-P-B and
                # A-Q-B). Neither chain contains the other's intermediate
                # endpoint, so a complete reverse path per edge cannot exist.
                # A closed degree-two component that lies on one line proves
                # the common boundary without inferring a nearby attachment.
                neighbors=defaultdict(set)
                for a,b in boundary:neighbors[int(a)].add(int(b));neighbors[int(b)].add(int(a))
                todo=set(neighbors)
                while todo:
                    component=set();queue=[min(todo)]
                    for node in queue:
                        if node in component:continue
                        component.add(node);queue.extend(neighbors[node]-component)
                    todo-=component
                    if len(component)<4 or len(component)>256 or any(len(neighbors[i])!=2 for i in component):continue
                    ids=np.array(sorted(component),int);q=v[ids];axis=int(np.argmax(np.ptp(q,axis=0)))
                    a=q[np.argmin(q[:,axis])];b=q[np.argmax(q[:,axis])];direction=b-a;length=float(np.linalg.norm(direction))
                    if length==0:continue
                    fraction=(q-a)@direction/(length*length);distance=np.linalg.norm(q-(a+fraction[:,None]*direction),axis=1)
                    residual=float(distance.max())
                    if residual>tolerance:continue
                    for control in ids:cycle_controls[int(control)]=ids
                    cycle_proofs.append(dict(cell=cell,round=iteration,controls=ids.tolist(),maximum_line_residual_m=residual,
                        closed_degree_two_boundary=True,positive_area_caps_added=0))
            for face_id,t in enumerate(faces):
                for a,b in zip(t,np.roll(t,-1)):
                    a,b=int(a),int(b)
                    if tuple(sorted((a,b)))not in edge_set:continue
                    d=v[b]-v[a];length=float(np.linalg.norm(d))
                    if length==0:continue
                    candidates=controls[tree.query_ball_point((v[a]+v[b])*.5,length*.5+tolerance)]
                    q=(v[candidates]-v[a])@d/(length*length)
                    distance=np.linalg.norm(v[candidates]-(v[a]+q[:,None]*d),axis=1)
                    # Proximity alone cannot establish the same boundary.
                    # Require an actual oppositely directed, monotone boundary
                    # chain joining the edge's endpoints before subdividing.
                    allowed={int(i)for i,t,distance in zip(candidates,q,distance)
                             if 0.<=t<=1. and distance<=tolerance}
                    allowed.update([a,b]);fractions=dict(zip(map(int,candidates),q));fractions[a],fractions[b]=0.,1.
                    paths={b:[b]};queue=[b]
                    for node in queue:
                        if a in paths:break
                        for peer in adjacency[node]:
                            if peer in allowed and peer not in paths and fractions[peer]<fractions[node]:
                                paths[peer]=paths[node]+[peer];queue.append(peer)
                    ids=paths.get(a,[])[1:-1]
                    if not len(ids)and a in cycle_controls and b in cycle_controls:
                        ids=[int(i)for i,s,distance in zip(candidates,q,distance)
                             if int(i)in cycle_controls[a]and 0.<s<1. and distance<=tolerance]
                    if not len(ids):continue
                    sid,ordinal=owners[face_id];pending.setdefault((sid,ordinal),set()).update(map(int,ids))
        if not pending:
            proof=dict(rounds=iteration,edgeRefinements=edits,coordinatesMoved=False,newControls=0,positiveAreaCapsAdded=0,
                collinearBoundaryCycleProofs=cycle_proofs)
            if ancestry is not None:proof['sourceFaceParents']=ancestry
            return out,proof
        additions=defaultdict(list);parent_additions=defaultdict(list)
        for (sid,ordinal),controls in pending.items():
            children=[out[sid][ordinal].copy()]
            for control in sorted(controls):
                if any(control in t for t in children):continue
                options=[];eligible_edge=False
                for i,t in enumerate(children):
                    for k in range(3):
                        a,b,c=map(int,np.roll(t,-k));d=v[b]-v[a];length=float(np.linalg.norm(d))
                        if length==0:continue
                        fraction=float((v[control]-v[a])@d/(length*length));distance=float(np.linalg.norm(v[control]-v[a]-fraction*d))
                        if 0.<fraction<1. and distance<=tolerance:
                            eligible_edge=True
                            new=np.array([[a,control,c],[control,b,c]],int)
                            old_normal=np.cross(v[b]-v[a],v[c]-v[a]);p=v[new];nn=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
                            if np.all(nn@old_normal>0):options.append((distance,i,k,new))
                if not options:
                    if eligible_edge:raise ValueError('Boundary subdivision reversed a source fragment')
                    continue
                # Thin fragments can put multiple edges within the numerical
                # line bound. Use the nearest valid edge across all children;
                # source order must not select a different, reversing edge.
                distance,i,_,new=min(options,key=lambda option:option[:3])
                children[i]=new[0];children.append(new[1]);edits.append(dict(surface=sid,parentFace=ordinal,existingControl=control,sourceLineResidualMeters=distance))
            out[sid][ordinal]=children[0];additions[sid].extend(children[1:])
            if ancestry is not None:parent_additions[sid].extend([ancestry[sid][ordinal]]*(len(children)-1))
        for sid,fs in additions.items():
            if fs:
                out[sid]=np.r_[out[sid],np.asarray(fs,int)]
                if ancestry is not None:ancestry[sid]=np.r_[ancestry[sid],np.asarray(parent_additions[sid],int)]
    raise ValueError('Boundary conformity did not close within the explicit finite refinement budget')


def reduce_encoded_surface_chains(positions,surfaces):
    """Cancel only exactly zero/opposed triangles on each canonical chart.

    Open surface charts cannot use a whole-cell boundary-loop repair. This
    reduction does not fill a loop or silently drop a repeated positive face.
    """
    v=np.asarray(positions,float);out={};proof={}
    def key(t):
        q=tuple(map(int,t));return min(q,q[1:]+q[:1],q[2:]+q[:2])
    for sid,faces in surfaces.items():
        f=np.asarray(faces,int);p=v[f];zero=np.all(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])==0,axis=1)
        buckets=defaultdict(list)
        for i in np.flatnonzero(~zero):buckets[key(f[i])].append(int(i))
        seen=set();selected=[];cancelled=0
        for k,ids in buckets.items():
            if k in seen:continue
            reverse=key((k[0],k[2],k[1]));other=buckets.get(reverse,[]);seen.update((k,reverse))
            coefficient=len(ids)-len(other);cancelled+=2*min(len(ids),len(other))
            if abs(coefficient)>1:raise ValueError('A canonical chart has repeated positive-area coverage')
            if coefficient:selected.append(ids[0]if coefficient>0 else other[0])
        out[sid]=f[selected]
        if not len(out[sid]):raise ValueError('A whole canonical chart vanished at the selected storage precision')
        proof[sid]=dict(exactZeroAreaTrianglesRemoved=int(zero.sum()),opposedCoincidentTrianglesCancelled=cancelled,coordinatesMoved=False)
    return out,proof
