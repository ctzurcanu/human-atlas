"""Conform all surfaces to point identifications already forced by storage.

No packed coordinate is changed. A triangle with two equal stored points has
zero area and cannot be emitted; every incident surface uses the same global
identification. Nondegenerate inversions are left to the existing strict gate.
"""
import numpy as np

def conform_quantized_surfaces(positions,surfaces,origin):
    v=np.asarray(positions,float);o=np.asarray(origin,float);q=(v-o).astype('<f4')
    return conform_stored_surfaces(q,surfaces,source_positions=v)

def conform_stored_surfaces(stored_positions,surfaces,*,source_positions=None):
    """Conform literal decoded controls, including multiple chart origins.

    No additional Float32 cast is allowed here: casting world positions would
    destroy the extra precision provided by each chart's local origin.
    """
    q=np.asarray(stored_positions)
    v=np.asarray(source_positions if source_positions is not None else q,float)
    if q.ndim!=2 or q.shape[1]!=3 or v.shape!=q.shape or not np.isfinite(q).all():
        raise ValueError('Explicit finite stored and source control banks are required')
    _,first,inverse=np.unique(q,axis=0,return_index=True,return_inverse=True);canonical=first[inverse];out={};parents={};dropped={}
    for sid,f in surfaces.items():
        mapped=canonical[np.asarray(f,int)];keep=(mapped[:,0]!=mapped[:,1])&(mapped[:,1]!=mapped[:,2])&(mapped[:,2]!=mapped[:,0]);out[sid]=mapped[keep];parents[sid]=np.flatnonzero(keep);dropped[sid]=int((~keep).sum())
        if not len(out[sid]):raise ValueError('A whole tissue surface vanished at the chosen coordinate precision')
    assert np.array_equal(q[canonical],q),'Coordinate storage changed during conformity'
    return out,parents,dict(storedCoordinatesChanged=False,zeroAreaFacesRemoved={k:n for k,n in dropped.items()if n},maximumEquivalentSourcePointSeparationMeters=float(np.linalg.norm(v-v[canonical],axis=1).max()))


def contract_declared_intersection_edges(positions, surfaces, stored_positions, edges, *, protected_prefix):
    """Identify only explicitly diagnosed construction edges below storage error.

    Both controls must be newly computed intersections, joined by an existing
    edge. Their separation must be no larger than their measured encoding
    errors combined. All incident tissue surfaces use the same identification;
    no source control moves and no unrelated nearby surfaces are welded.
    This is a candidate operation: topology, parent-plane, shared-interface,
    occupied-volume and whole-surface checks must still be rerun.
    """
    p, q = np.asarray(positions, float), np.asarray(stored_positions, float)
    if p.shape != q.shape or p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(np.r_[p.ravel(), q.ravel()]).all():
        raise ValueError('Corresponding finite source and stored controls are required')
    output = {sid: np.array(f, int, copy=True) for sid, f in surfaces.items()}
    parents = {sid: np.arange(len(f)) for sid, f in output.items()}
    changes = []
    for keep, drop in edges:
        keep, drop = int(keep), int(drop)
        if keep == drop or min(keep, drop) < protected_prefix or max(keep, drop) >= len(p):
            raise ValueError('Only two distinct new intersection controls may be contracted')
        if any(keep in [e['keep'], e['drop']] or drop in [e['keep'], e['drop']] for e in changes):
            raise ValueError('Dependent edge contractions require a newly diagnosed candidate')
        incident = [sid for sid, f in output.items() if np.any(np.sum(np.isin(f, [keep, drop]), axis=1) == 2)]
        if not incident:
            raise ValueError('Close controls are not one explicitly diagnosed construction edge')
        length = float(np.linalg.norm(p[keep]-p[drop]))
        bound = float(np.linalg.norm(q[[keep, drop]]-p[[keep, drop]], axis=1).sum())
        if not 0 < length <= bound:
            raise ValueError('Construction edge exceeds its actual combined storage error')
        removed = {}
        for sid, f in output.items():
            mapped = np.where(f == drop, keep, f)
            valid = (mapped[:, 0] != mapped[:, 1]) & (mapped[:, 1] != mapped[:, 2]) & (mapped[:, 2] != mapped[:, 0])
            removed[sid] = int((~valid).sum())
            if not valid.any():
                raise ValueError('A whole tissue interface would disappear')
            output[sid], parents[sid] = mapped[valid], parents[sid][valid]
        changes.append(dict(keep=keep, drop=drop, source_edge_length_m=length,
                            measured_combined_encoding_error_m=bound,
                            incident_surfaces=incident, removed_collapsed_faces={s: n for s, n in removed.items() if n}))
    return output, parents, dict(changes=changes, literal_source_controls_moved=False,
                                 unrelated_points_identified=False, full_geometry_acceptance_required=True)

def contract_quantized_intersection_edges(positions,surfaces,origin,mutable_start,mutable_end,maximum_edge_m=1e-8,maximum_contractions=32):
    """Contract only nanometre Boolean slivers which invert during storage.

    Literal source controls are protected. The same edge identification is
    propagated through every surface; topology and crossings still require
    independent checks after contraction. No general source simplification.
    """
    if not isinstance(maximum_contractions,int)or maximum_contractions<1:raise ValueError('A positive finite contraction budget is required')
    v=np.asarray(positions,float);o=np.asarray(origin,float);out={sid:f.copy()for sid,f in surfaces.items()};parents={sid:np.arange(len(f))for sid,f in out.items()};edits=[]
    for iteration in range(maximum_contractions):
        q=(v-o).astype('<f4').astype(float)+o;bad_faces=[]
        for fs in out.values():
            t=v[fs];u=q[fs];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);m=np.cross(u[:,1]-u[:,0],u[:,2]-u[:,0]);bad=(n*m).sum(1)<=0;bad_faces.extend(fs[bad])
        if not bad_faces:return out,parents,dict(intersectionOnlyEdgeContractions=edits,maximumSourceEdgeMeters=max((e['lengthMeters']for e in edits),default=0),literalSourceControlsMoved=False)
        choices=[]
        for face in bad_faces:
            for a,b in zip(face,np.roll(face,-1)):
                length=float(np.linalg.norm(v[a]-v[b]));ma=mutable_start<=a<mutable_end;mb=mutable_start<=b<mutable_end
                if 0<length<=maximum_edge_m and(ma or mb):
                    keep,drop=(a,b)if mb else(b,a);choices.append((length,int(keep),int(drop)))
        if not choices:
            return out,parents,dict(intersectionOnlyEdgeContractions=edits,maximumSourceEdgeMeters=max((e['lengthMeters']for e in edits),default=0),literalSourceControlsMoved=False,remainingFacesForStrictDiagonalRepair=len(bad_faces))
        length,keep,drop=min(choices);edits.append(dict(keep=keep,drop=drop,lengthMeters=length,maximumPackedCoordinateChangeMeters=float(np.linalg.norm(q[keep]-q[drop]))))
        for sid,fs in out.items():
            fs=np.where(fs==drop,keep,fs);good=(fs[:,0]!=fs[:,1])&(fs[:,1]!=fs[:,2])&(fs[:,2]!=fs[:,0]);out[sid]=fs[good];parents[sid]=parents[sid][good]
            if not len(out[sid]):raise ValueError('A complete surface vanished during intersection contraction')
    raise ValueError('Intersection-only contraction did not converge')
