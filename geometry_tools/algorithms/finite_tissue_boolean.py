"""Provenance-bound subtraction of a closed occupied volume from finite cells.

No anatomical acceptance is implied. Classifier caps are construction data, not
anatomical surfaces. Original controls remain a literal prefix; only newly
computed intersection controls may be canonicalized within the declared tiny
numerical tolerance. Reciprocal tissue surfaces are stored once. Incompatible
boundary subdivisions fail rather than being filled or covered by a disk.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path
import itertools
import json
import numpy as np
from scipy.spatial import cKDTree
import manifold3d as md


class FiniteBooleanError(ValueError):
    pass


def oriented_key(triangle):
    a = tuple(map(int, triangle))
    return min(a, a[1:]+a[:1], a[2:]+a[:2])


def boundary_signature(faces):
    """Directed boundary after cancelling opposite internal edges."""
    edges = Counter((int(t[k]), int(t[(k+1)%3])) for t in faces for k in range(3))
    return Counter({e:n-edges[(e[1],e[0])] for e,n in edges.items() if n>edges[(e[1],e[0])]})


def same_oriented_triangles(first,second):
    """Exact oriented triangle multiset equality, independent of row order.

    Cyclic corner rotations preserve orientation; reversed faces do not.
    This is stronger than comparing a Boolean engine's re-imported volume.
    """
    a,b=np.asarray(first,np.int64),np.asarray(second,np.int64)
    if a.shape!=b.shape:return False
    def canonical(f):
        shift=np.argmin(f,axis=1)
        f=np.take_along_axis(f,(shift[:,None]+np.arange(3))%3,axis=1)
        return f[np.lexsort((f[:,2],f[:,1],f[:,0]))]
    return np.array_equal(canonical(a),canonical(b))


def make_solid(positions, faces, original_id=None, properties=None):
    """Exact-coordinate welding only; face IDs remain local input face IDs.

    When attributes are supplied, chart vertices remain separate and merge
    vectors establish only the geometric topology. Attributes are never
    recalculated by Manifold, including original integer normal words.
    """
    positions=np.asarray(positions,np.float64);faces=np.asarray(faces,np.uint64)
    if original_id is None:original_id=md.Manifold.reserve_ids(1)
    if properties is None:
        used=np.unique(faces);vp,inverse=np.unique(positions[used],axis=0,return_inverse=True)
        lookup=np.full(len(positions),-1,dtype=np.int64);lookup[used]=inverse
        triangles=lookup[faces].astype(np.uint64);kwargs={}
    else:
        corner=positions[faces].reshape(-1,3)
        attrs=np.asarray(properties,np.float64).reshape(len(corner),-1)
        vp=np.c_[corner,attrs];triangles=np.arange(len(corner),dtype=np.uint64).reshape(-1,3)
        _,first,inverse=np.unique(corner,axis=0,return_index=True,return_inverse=True)
        target=first[inverse];source=np.flatnonzero(target!=np.arange(len(corner)))
        kwargs={'merge_from_vert':source.astype(np.uint64),'merge_to_vert':target[source].astype(np.uint64)}
    mesh=md.Mesh64(np.ascontiguousarray(vp),np.ascontiguousarray(triangles),
                   run_index=np.array([0,faces.size],np.uint64),run_original_id=np.array([original_id],np.uint32),
                   face_id=np.arange(len(faces),dtype=np.uint64),**kwargs)
    solid=md.Manifold(mesh)
    if solid.status()!=md.Error.NoError:raise FiniteBooleanError(f'Input is not a closed oriented cell: {solid.status()}')
    return solid,int(original_id)


def provenance(mesh):
    ids=np.empty(len(mesh.tri_verts),np.uint32)
    runs=np.asarray(mesh.run_index)//3
    originals=np.asarray(mesh.run_original_id)
    for i,original in enumerate(originals):ids[runs[i]:runs[i+1]]=original
    return ids,np.asarray(mesh.face_id,dtype=np.int64)


class CommonPositionBank:
    """Literal source prefix + bounded, explicit intersection-only welding."""
    def __init__(self,source,tolerance=2e-12):
        self.source=np.asarray(source,np.float64).copy();self.rows=self.source.tolist()
        self.tolerance=float(tolerance);self.tree=cKDTree(self.source)
        self.exact={tuple(v):i for i,v in reversed(list(enumerate(self.source)))}
        self.buckets=defaultdict(list);self.maximum_snap=0.;self.new_snaps=0
    def insert(self,p):
        p=np.asarray(p,float);key=tuple(p)
        if key in self.exact:return self.exact[key]
        if self.tolerance==0:
            i=len(self.rows);self.rows.append(p.tolist());self.exact[key]=i;return i
        distance,index=self.tree.query(p)
        if distance<=self.tolerance:
            self.maximum_snap=max(self.maximum_snap,float(distance));return int(index)
        bucket=tuple(np.floor(p/self.tolerance).astype(np.int64))
        for delta in itertools.product((-1,0,1),repeat=3):
            for i in self.buckets[tuple(bucket[k]+delta[k] for k in range(3))]:
                distance=np.linalg.norm(p-self.rows[i])
                if distance<=self.tolerance:
                    self.maximum_snap=max(self.maximum_snap,float(distance));self.new_snaps+=1;return i
        i=len(self.rows);self.rows.append(p.tolist());self.exact[key]=i;self.buckets[bucket].append(i);return i
    @property
    def positions(self):return np.asarray(self.rows,np.float64)


def barycentric(points,triangle):
    """Stable local-coordinate parent chart; zero-area parents are rejected."""
    triangle=np.asarray(triangle,float);points=np.asarray(points,float)
    e=np.stack((triangle[1]-triangle[0],triangle[2]-triangle[0]),axis=1)
    uv=np.linalg.lstsq(e,(points-triangle[0]).T,rcond=None)[0].T
    return np.c_[1-uv.sum(1),uv]


def load_cells(archive,manifest):
    z=np.load(archive);meta=json.loads(Path(manifest).read_text());surfaces={r['id']:z[r['array']] for r in meta['surfaces']}
    rows={r['id']:r for r in meta['surfaces']};cells={};parents={}
    for name,uses in meta['cells'].items():
        triangles=[];provenance_rows=[]
        for use in uses:
            fs=surfaces[use['surface']];side=use['side']
            triangles.append(fs if side==1 else fs[:,[0,2,1]])
            provenance_rows.extend((use['surface'],i,side) for i in range(len(fs)))
        cells[name]=np.concatenate(triangles);parents[name]=provenance_rows
    return z,meta,surfaces,rows,cells,parents


def classifier_caps_outside(occupied_positions,occupied_faces,cap_faces,host_positions,host_outer_faces,tolerance=2e-12):
    """Complete finite cap/host separation, not cap-center sampling.

    SAT includes in-plane axes for coplanar triangles; point winding establishes
    the outside component once every cap triangle is disjoint from the host.
    Touching or ambiguous caps are rejected, never silently dropped.
    """
    from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains
    caps=np.asarray(cap_faces,np.int64).ravel()
    if not len(caps):return {'capFaces':0,'verticesInside':0,'ambiguousVertices':0,'triangleContacts':0}
    if np.any(caps<0) or np.any(caps>=len(occupied_faces)):raise FiniteBooleanError('Classifier cap face ID outside input')
    cap_tri=np.asarray(occupied_positions,float)[np.asarray(occupied_faces)[caps]]
    host_tri=np.asarray(host_positions,float)[np.asarray(host_outer_faces)]
    points=np.unique(cap_tri.reshape(-1,3),axis=0);inside,ambiguous=ray_contains(points,host_tri)
    if np.any(inside|ambiguous):raise FiniteBooleanError('Classifier cap has inside or ambiguous vertices')
    low,high=host_tri.min(1),host_tri.max(1);tested=0
    for a in cap_tri:
        ids=np.flatnonzero(np.all(high>=a.min(0)-tolerance,axis=1)&np.all(low<=a.max(0)+tolerance,axis=1))
        ea=np.roll(a,-1,axis=0)-a;na=np.cross(ea[0],ea[1])
        if np.linalg.norm(na)==0:raise FiniteBooleanError('Degenerate classifier cap')
        for i in ids:
            tested+=1;b=host_tri[i];eb=np.roll(b,-1,axis=0)-b;nb=np.cross(eb[0],eb[1])
            axes=np.r_[na[None],nb[None],np.cross(ea[:,None],eb[None,:]).reshape(-1,3),np.cross(na,ea),np.cross(nb,eb)]
            length=np.linalg.norm(axes,axis=1);axes=axes[length>1e-25]/length[length>1e-25,None]
            # Translate projections locally to avoid metre-origin cancellation.
            ap=(a-a[0])@axes.T;bp=(b-a[0])@axes.T
            separated=(ap.max(0)<bp.min(0)-tolerance)|(bp.max(0)<ap.min(0)-tolerance)
            if not np.any(separated):raise FiniteBooleanError(f'Classifier cap touches/crosses host outer face {i}')
    return {'capFaces':len(caps),'verticesInside':0,'ambiguousVertices':0,'triangleContacts':0,'finiteTriangleCandidateTests':tested,'method':'Complete triangle SAT + oriented finite-ray component test'}


def subtract_occupied(positions,surfaces,surface_rows,cells,parents,cell_uses,
                      occupied_positions,occupied_faces,classifier_cap_faces=(),
                      intersection_tolerance=2e-12,face_ancestry=None,occupied_face_owners=None,partition_free_outer=False,
                      construction_resolution_meters=None,reconcile_shared_subdivisions=False,
                      diagnostic_directory=None,exact_source_replay=False):
    """Return complementary cells and one canonical contact surface per cell.

    The returned contact surfaces are outward from tissue (inward to occupied
    volume). ``vascularContacts`` records the reverse use for the vascular wall
    assembler; it is NOT an additional endocrine cell or an extra tissue coat.
    Normal charts are represented by exact source-parent face/barycentric data.
    A later packer must emit each chart once and apply reciprocal normal signs.
    Opt-in subdivision reconciliation accepts only the same directed boundary
    with additional collinear controls. It propagates those controls through
    every incident chart; gaps, orientation changes and different curves fail.
    """
    occupied_positions=np.asarray(occupied_positions,float);occupied_faces=np.asarray(occupied_faces,np.int64)
    if construction_resolution_meters is not None:
        construction_resolution_meters=float(construction_resolution_meters)
        if not np.isfinite(construction_resolution_meters)or construction_resolution_meters<=0:
            raise ValueError('An explicit positive measured construction resolution is required')
    caps=set(map(int,np.asarray(classifier_cap_faces).ravel()))
    outer_faces=np.concatenate([surfaces[sid] for sid,row in surface_rows.items() if row["kind"]=="outer"])
    cap_proof=classifier_caps_outside(occupied_positions,occupied_faces,list(caps),positions,outer_faces,intersection_tolerance)
    vascular_owners=np.full(len(occupied_faces),"vascular:occupied-wall",dtype=object) if occupied_face_owners is None else np.asarray(occupied_face_owners,dtype=object)
    if len(vascular_owners)!=len(occupied_faces):raise FiniteBooleanError("One vascular owner is required per occupied outer face")
    bank=CommonPositionBank(np.r_[positions,occupied_positions],0. if exact_source_replay else intersection_tolerance)
    cutter,cutter_id=make_solid(occupied_positions,occupied_faces)
    if construction_resolution_meters is not None:cutter=cutter.set_tolerance(construction_resolution_meters)
    host_ids={};fragments=defaultdict(dict);contacts={};volumes={};raw_output_faces={};source_replays=[]
    for name,fs in cells.items():
        host,host_id=make_solid(positions,fs);host_ids[name]=host_id
        if construction_resolution_meters is not None:host=host.set_tolerance(construction_resolution_meters)
        result=host-cutter
        if result.status()!=md.Error.NoError:raise FiniteBooleanError(f'{name}: difference {result.status()}')
        out=result.to_mesh64();origins,faceids=provenance(out)
        ov=np.asarray(out.vert_properties)[:,:3]
        if exact_source_replay:
            from .exact_source_intersections import replay_boolean_vertices
            ov,replay=replay_boolean_vertices(ov,np.asarray(out.tri_verts),origins,faceids,
                host_id,positions,fs,cutter_id,occupied_positions,occupied_faces)
            source_replays.append(dict(cell=name,**replay))
        mapped=np.array([bank.insert(v)for v in ov],np.int64)
        output=mapped[np.asarray(out.tri_verts)];raw_output_faces[name]=output
        local=defaultdict(list);contact=[];contact_parent=[]
        for tri,origin,faceid in zip(output,origins,faceids):
            if origin==host_id:
                sid,ordinal,side=parents[name][int(faceid)]
                local[(sid,ordinal)].append(tri if side==1 else tri[[0,2,1]])
            elif origin==cutter_id:
                if int(faceid) in caps:raise FiniteBooleanError(f'Classifier cap {faceid} contributed to {name}; it must lie outside the host')
                contact.append(tri);contact_parent.append(int(faceid))
            else:raise FiniteBooleanError(f'Unknown Boolean provenance {origin}')
        for key,tris in local.items():fragments[key][name]=np.asarray(tris,np.int64)
        contacts[name]=(np.asarray(contact,np.int64).reshape(-1,3),np.asarray(contact_parent,np.int64))
        volumes[name]={'before':float(host.volume()),'after':float(result.volume())}
    output_surfaces={};parent_faces={};output_rows=[];evidence=[]
    shared_edge_plans=defaultdict(set);boundary_refinements=[]
    fragment_positions=bank.positions
    for sid,source in surfaces.items():
        owners=[name for name,uses in cell_uses.items() if any(u['surface']==sid for u in uses)]
        combined=[];parent=[]
        for ordinal in range(len(source)):
            versions=fragments.get((sid,ordinal),{})
            if set(versions)!=set(owners):
                # Complete removal is legal, partial absence at one reciprocal
                # owner is not. This also catches a classification disagreement.
                if versions:raise FiniteBooleanError(f'Owner mismatch at {sid}/{ordinal}: {list(versions)} vs {owners}')
                continue
            canonical=versions[owners[0]]
            for owner in owners[1:]:
                other=versions[owner]
                if boundary_signature(canonical)!=boundary_signature(other):
                    if not reconcile_shared_subdivisions:raise FiniteBooleanError(f'Incompatible reciprocal boundary at {sid}/{ordinal}')
                    from .shared_triangle_edges import reconcile_boundaries
                    try:
                        deviation=reconcile_boundaries(fragment_positions,canonical,other,shared_edge_plans,tolerance=intersection_tolerance)
                    except ValueError as error:
                        raise FiniteBooleanError(f'Incompatible reciprocal boundary at {sid}/{ordinal}: {error}') from error
                    boundary_refinements.append(dict(surface=sid,sourceFace=ordinal,owners=[owners[0],owner],maximumPointDeviationMeters=float(deviation)))
                pa=fragment_positions
                area=lambda f:np.linalg.norm(np.cross(pa[f[:,1]]-pa[f[:,0]],pa[f[:,2]]-pa[f[:,0]]),axis=1).sum()/2
                if abs(area(canonical)-area(other))>max(1e-18,area(canonical)*1e-9):raise FiniteBooleanError(f'Reciprocal area mismatch at {sid}/{ordinal}')
            combined.extend(canonical);parent.extend([ordinal]*len(canonical))
        fs=np.asarray(combined,np.int64).reshape(-1,3)
        output_surfaces[sid]=fs;parent_faces[sid]=np.asarray(parent,np.int64)
        output_rows.append({**surface_rows[sid],'array':f'surface{len(output_rows)}Faces','parentSurface':sid,'parentNormalSign':1})
        evidence.append({'surface':sid,'sourceFaces':len(source),'outputFaces':len(fs),'owners':owners,'reciprocalBoundaryMatched':True})
    new_uses={name:[dict(u) for u in uses if len(output_surfaces[u['surface']])] for name,uses in cell_uses.items()}
    vascular_contacts=[]
    for name,(fs,parent) in contacts.items():
        if not len(fs):continue
        for owner in sorted(set(vascular_owners[parent])):
            take=np.flatnonzero(vascular_owners[parent]==owner);triangles=fs[take];faceids=parent[take]
            sid='vascular-contact:'+name+':'+str(owner)
            output_surfaces[sid]=triangles;parent_faces[sid]=faceids
            output_rows.append({'id':sid,'kind':'tissue-vascular-contact','array':f'surface{len(output_rows)}Faces','parentSurface':'occupied-outer','parentNormalSign':-1,'vascularOwner':str(owner)})
            new_uses[name].append({'surface':sid,'side':1,'neighbor':str(owner)})
            vascular_contacts.append({'surface':sid,'side':-1,'neighbor':name,'vascularOwner':str(owner),
                                      'sourceFaceIds':faceids.tolist(),
                                      'sourceAncestry':None if face_ancestry is None else [face_ancestry[int(i)] for i in faceids]})
    free_surfaces=[]
    if partition_free_outer:
        whole,whole_id=make_solid(positions,outer_faces)
        if construction_resolution_meters is not None:whole=whole.set_tolerance(construction_resolution_meters)
        outside=cutter-whole
        if outside.status()!=md.Error.NoError:raise FiniteBooleanError('Outside vascular classifier failed')
        mesh=outside.to_mesh64();origin,faceids=provenance(mesh)
        output_positions=np.asarray(mesh.vert_properties)[:,:3]
        if exact_source_replay:
            output_positions,replay=replay_boolean_vertices(output_positions,np.asarray(mesh.tri_verts),origin,faceids,
                cutter_id,occupied_positions,occupied_faces,whole_id,positions,outer_faces)
            source_replays.append(dict(cell='occupied-free-boundary',**replay))
        mapped=np.array([bank.insert(v) for v in output_positions],np.int64)
        fs=mapped[np.asarray(mesh.tri_verts)]
        take=(origin==cutter_id)&~np.isin(faceids,list(caps))
        for owner in sorted(set(vascular_owners[faceids[take]])):
            ids=np.flatnonzero(take&(vascular_owners[np.clip(faceids,0,len(vascular_owners)-1)]==owner))
            sid='vascular-free:'+str(owner);output_surfaces[sid]=fs[ids];parent_faces[sid]=faceids[ids]
            output_rows.append({'id':sid,'kind':'vascular-free-exterior','array':f'surface{len(output_rows)}Faces','parentSurface':'occupied-outer','parentNormalSign':1,'vascularOwner':str(owner)})
            free_surfaces.append({'surface':sid,'side':1,'vascularOwner':str(owner)})
    shared_refinement=None
    if shared_edge_plans:
        from .shared_triangle_edges import refine_surfaces
        output_surfaces,refinement_parents,changed=refine_surfaces(bank,output_surfaces,shared_edge_plans)
        for sid,indices in refinement_parents.items():parent_faces[sid]=parent_faces[sid][indices]
        for contact in vascular_contacts:
            ids=parent_faces[contact['surface']]
            contact['sourceFaceIds']=ids.tolist()
            contact['sourceAncestry']=None if face_ancestry is None else [face_ancestry[int(i)] for i in ids]
        for item in evidence:item['outputFaces']=len(output_surfaces[item['surface']])
        shared_refinement=dict(method='Common collinear edge subdivisions propagated through every incident output chart',
            edges=len(shared_edge_plans),changedFaces=changed,boundaryChecks=boundary_refinements,
            maximumPointDeviationMeters=max(r['maximumPointDeviationMeters']for r in boundary_refinements),
            originalControlsMoved=False,distinctCurvesJoined=False)
    outer_replay=None
    if partition_free_outer:
        boundary=[output_surfaces[u['surface']][:,[0,2,1]] for u in vascular_contacts]
        boundary.extend(output_surfaces[u['surface']] for u in free_surfaces)
        # These disks belong ONLY to the computational replay certificate.
        # They are absent from all returned anatomical surfaces/cell uses.
        cap_controls=np.array([bank.insert(p) for p in occupied_positions],np.int64)
        proof_caps=cap_controls[occupied_faces[np.array(sorted(caps),np.int64)]]
        boundary.append(proof_caps)
        replay_faces=np.concatenate(boundary)
        replay,_=make_solid(bank.positions,replay_faces)
        if abs(replay.volume()-cutter.volume())>max(1e-19,abs(cutter.volume())*1e-8):raise FiniteBooleanError('Shared contact/free outer replay changed occupied volume')
        area_original=np.linalg.norm(np.cross(occupied_positions[occupied_faces[:,1]]-occupied_positions[occupied_faces[:,0]],occupied_positions[occupied_faces[:,2]]-occupied_positions[occupied_faces[:,0]]),axis=1)/2
        area_replayed=np.zeros(len(occupied_faces));fragment_perimeter=np.zeros(len(occupied_faces));bp=bank.positions
        new_occupied_surfaces={u['surface']for u in vascular_contacts+free_surfaces}
        for sid,fs in output_surfaces.items():
            if sid not in new_occupied_surfaces:continue
            ar=np.linalg.norm(np.cross(bp[fs[:,1]]-bp[fs[:,0]],bp[fs[:,2]]-bp[fs[:,0]]),axis=1)/2
            np.add.at(area_replayed,parent_faces[sid],ar)
            perimeter=np.linalg.norm(bp[fs]-np.roll(bp[fs],-1,axis=1),axis=2).sum(1)
            np.add.at(fragment_perimeter,parent_faces[sid],perimeter)
        area_replayed[list(caps)]=area_original[list(caps)]
        residual=np.max(np.abs(area_replayed-area_original))
        source_perimeter=np.linalg.norm(occupied_positions[occupied_faces]-np.roll(occupied_positions[occupied_faces],-1,axis=1),axis=2).sum(1)
        # The area perturbation is bounded by displacement × perimeter, not a
        # fixed arbitrary relative tolerance on extremely narrow triangles.
        coordinate_bound=intersection_tolerance+64*np.finfo(float).eps*np.max(np.abs(occupied_positions))
        area_bound=coordinate_bound*(source_perimeter+fragment_perimeter)+1e-18
        bad=np.flatnonzero(np.abs(area_replayed-area_original)>area_bound)
        if len(bad):raise FiniteBooleanError(f'Occupied outer parent-face coverage failed: {len(bad)} faces; maximum area residual {residual:.12g} m²; examples '+str([(int(i),float(area_original[i]),float(area_replayed[i])) for i in bad[:8]]))
        outer_replay={'closedClassifierReplay':True,'classifierCapsEmittedAnatomically':False,
                      'sourceVolumeMetersCubed':float(cutter.volume()),'replayedVolumeMetersCubed':float(replay.volume()),
                      'maximumParentTriangleAreaResidualMetersSquared':float(residual),
                      'areaResidualBoundMethod':'Declared coordinate perturbation times source and fragment perimeters',
                      'maximumAreaResidualExcessAboveBoundMetersSquared':float(np.maximum(np.abs(area_replayed-area_original)-area_bound,0).max())}
    # Reconstructed cells use the canonical shared triangulation, not each
    # independent Boolean's potentially different internal diagonal choices.
    reconstructed={name:np.concatenate([output_surfaces[u['surface']] if u['side']==1 else output_surfaces[u['surface']][:,[0,2,1]] for u in uses]) for name,uses in new_uses.items()}
    if diagnostic_directory is not None:
        from geometry_tools.common import atomic_json,atomic_npz
        diagnostic_directory=Path(diagnostic_directory)
        diagnostic_directory.mkdir(exist_ok=False)
        arrays={'positions':bank.positions}
        diagnostic_names=list(reconstructed)
        for i,name in enumerate(diagnostic_names):
            arrays[f'canonical_{i}']=reconstructed[name]
            arrays[f'independent_{i}']=raw_output_faces[name]
        atomic_npz(diagnostic_directory/'replay.npz',**arrays)
        atomic_json(diagnostic_directory/'replay.json',{'names':diagnostic_names,'volumes':volumes})
    for name,fs in reconstructed.items():
        # Re-importing a very thin, already closed Boolean result can cause
        # the engine to simplify its boundary. When the canonical and raw
        # output have exactly the same oriented triangles, equality itself
        # proves replay: a second import is neither needed nor stronger.
        if same_oriented_triangles(fs,raw_output_faces[name]):
            volumes[name]['canonical']=volumes[name]['after']
            volumes[name]['canonicalReplay']='exact oriented triangle multiset'
            continue
        solid,_=make_solid(bank.positions,fs)
        if abs(solid.volume()-volumes[name]['after'])>max(1e-19,abs(solid.volume())*1e-8):raise FiniteBooleanError(f'Canonical reconstruction volume changed: {name}; canonical={solid.volume():.17g} m³; independent={volumes[name]["after"]:.17g} m³; source={volumes[name]["before"]:.17g} m³')
        volumes[name]['canonical']=float(solid.volume())
    # Expose source-parent barycentric charts. Untouched source corners can be
    # detected literally and use original packed words without any interpolation.
    charts={};final_positions=bank.positions
    new_occupied_surfaces={u['surface']for u in vascular_contacts+free_surfaces}
    for sid,fs in output_surfaces.items():
        parent=parent_faces[sid];new_occupied=sid in new_occupied_surfaces
        original=occupied_faces if new_occupied else surfaces[sid]
        parent_positions=occupied_positions if new_occupied else positions
        weights=np.array([barycentric(final_positions[tri],parent_positions[original[int(k)]]) for tri,k in zip(fs,parent)])
        charts[sid]={'parentFaces':parent,'barycentrics':weights.reshape(-1,3,3)}
    return {'positions':bank.positions,'surfaces':output_surfaces,'surfaceRows':output_rows,'cells':new_uses,
            'cellFaces':reconstructed,'charts':charts,'vascularContacts':vascular_contacts,'vascularFreeSurfaces':free_surfaces,
            'evidence':{'classifierCapContributions':0,'classifierCapsOutsideProof':cap_proof,'vascularOuterReplay':outer_replay,'interfaceEvidence':evidence,'cellVolumes':volumes,
                        'exactSourceReplays':source_replays,
                        'sharedSubdivisionReconciliation':shared_refinement,
                        'literalSourcePrefix':np.array_equal(bank.positions[:len(positions)],positions),
                        'intersectionCanonicalToleranceMeters':intersection_tolerance,'maximumCanonicalSnapMeters':bank.maximum_snap,
                        'explicitConstructionResolutionMeters':construction_resolution_meters,
                        'intersectionSnapCount':bank.new_snaps,'sourceControls':len(positions),'outputControls':len(bank.positions)},
            'limitations':'Construction geometry only. Vessel wall/lumen assembly, full contact self-intersection tests, render and anatomical acceptance remain separate gates.'}


def shared_corner_normal_words(result,source_positions,source_surfaces,source_outer_normal_words):
    """One canonical per-corner chart, preserving literal exterior word triples.

    Internal and vascular charts are flat. Newly introduced exterior controls
    interpolate only their original source triangle's normals; original corners
    retain their raw words without normalizing them again.
    """
    source_positions=np.asarray(source_positions,float);words=np.asarray(source_outer_normal_words,np.int16)
    out={};literal_counts={};minimum_cosine={}
    for row in result['surfaceRows']:
        sid=row['id'];fs=result['surfaces'][sid];tri=result['positions'][fs]
        area=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);unit=area/np.linalg.norm(area,axis=1)[:,None]
        normals=np.repeat(np.rint(unit*32767).astype(np.int16)[:,None,:],3,axis=1);literal=0
        if row['kind']=='outer':
            parent=result['charts'][sid]['parentFaces'];bary=result['charts'][sid]['barycentrics'];source=source_surfaces[sid][parent]
            values=np.einsum('tck,tkj->tcj',bary,words[source].astype(float))
            values/=np.linalg.norm(values,axis=2)[:,:,None]
            normals=np.rint(values*32767).astype(np.int16)
            for i,t in enumerate(tri):
                for k,p in enumerate(t):
                    exact=np.flatnonzero(np.all(source_positions[source[i]]==p,axis=1))
                    if len(exact):normals[i,k]=words[source[i,int(exact[0])]];literal+=1
        shade=normals.astype(float);dots=np.einsum('tki,ti->tk',shade,unit)/np.linalg.norm(shade,axis=2)
        if len(dots) and np.min(dots)<=0:raise FiniteBooleanError(f'Nonoutward shared chart at {sid}')
        out[sid]=normals;literal_counts[sid]=literal;minimum_cosine[sid]=float(dots.min()) if len(dots) else None
    return out,{'literalOuterCorners':literal_counts,'minimumCornerNormalFaceCosine':minimum_cosine,'nonOutwardCorners':0}
