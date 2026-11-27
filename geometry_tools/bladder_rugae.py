"""Finite relief of the existing bladder lining, with a fixed trigone/outer wall.

Ridges deform the owned mucosal tissue, not a second overlapping fold mesh.
Heights/spacing are illustrative review inputs. No clinical opening or filling
state is inferred. The result still needs whole-fluid, source-interface,
packing, rendered and anatomical review before installation.
"""
from pathlib import Path
import json
import numpy as np

from .algorithms.finite_tissue_boolean import CommonPositionBank, barycentric
from .algorithms.shared_triangle_edges import bisect_surfaces
from .algorithms.surface_normals import grouped_normal_words, source_normal_words
from .common import atomic_json, atomic_npz, digest_file
from .model import AtlasReader, Mesh
from .spatial import BVH
from .tissue import pack


def smoothstep(value):
    t = np.clip(value, 0., 1.)
    return t*t*t*(10.+t*(-15.+6.*t))


def fold_height(points, distance_to_fixed_lining, *, center, height_m=.00035,
                spacing_m=.004, transition_m=.004, neck_y_m=.899):
    points = np.asarray(points, float)
    distance = np.asarray(distance_to_fixed_lining, float)
    center = np.asarray(center, float)
    if points.ndim != 2 or points.shape[1] != 3 or center.shape != (3,) or distance.shape != (len(points),):
        raise ValueError('Corresponding physical surface points and fixed-support distances are required')
    if not np.isfinite(np.r_[points.ravel(), center, distance]).all() or np.any(distance < 0):
        raise ValueError('Finite source coordinates and nonnegative distances are required')
    if not 0 < height_m <= .0006 or not .003 <= spacing_m <= .007 or transition_m < .003:
        raise ValueError('Explicit bounded mucosal relief inputs are required')
    p = (points-center)/spacing_m
    first = 2*np.pi*(p[:, 0]+.3*p[:, 1])+.65*np.sin(2*np.pi*p[:, 1]/3)+.4*np.sin(2*np.pi*p[:, 2]/4)
    second = 2*np.pi*(.75*p[:, 2]-.4*p[:, 1])+.4*np.sin(2*np.pi*p[:, 0]/4)
    a = np.maximum(0., np.cos(first))**2
    b = .65*np.maximum(0., np.cos(second))**2
    ridge = 1.-(1.-a)*(1.-b)
    return height_m*ridge*smoothstep(distance/transition_m)*smoothstep((points[:, 1]-neck_y_m)/transition_m)


def fixed_source_fields(source_positions, encoded_positions, surfaces, moving, source, *, parent_hints=None, parent_faces_out=None):
    """Restore static facets from the installed field, never a failed trial.

    Closest centroids only propose parent facets. Every actual corner must
    independently lie on its proposed source triangle at Float64 roundoff.
    """
    from .distance import point_triangle
    with np.load(source['source'], allow_pickle=False) as z:
        old_p=z['encoded_positions'].copy()
        old_faces={sid:z[key].copy()for sid,key in source['surfaces'].items()}
        old_words={sid:z[key].copy()for sid,key in source['corner_words'].items()}
    fields={sid:np.zeros((len(f),3,3),dtype='<i2')for sid,f in surfaces.items()}
    evidence=[]
    bound=128*np.finfo(float).eps*max(float(np.abs(old_p).max()),np.finfo(float).tiny)
    for sid,faces in surfaces.items():
        take=np.flatnonzero(~np.any(moving[faces],axis=1))
        if not len(take):continue
        fs=faces[take];native=old_faces[sid]
        lookup={tuple(sorted(map(int,f))):i for i,f in enumerate(native)}
        parents=np.array([lookup.get(tuple(sorted(map(int,f))),-1)for f in fs])
        if parent_hints is not None and sid in parent_hints:
            hints=np.asarray(parent_hints[sid],int)
            if hints.shape!=(len(faces),)or np.any((hints < -1)|(hints>=len(native))):
                raise ValueError('Invalid source-facet lineage hints')
            use=hints[take]>=0;parents[use]=hints[take][use]
        missing=np.flatnonzero(parents<0)
        if len(missing):
            query,_=BVH(Mesh(old_p,native,source['frame'])).closest_points(source_positions[fs[missing]].mean(1),backend='auto')
            parents[missing]=query['face_ids']
        triangles=old_p[native[parents]]
        distance,_,_=point_triangle(source_positions[fs].reshape(-1,3),np.repeat(triangles,3,axis=0))
        maximum=float(np.sqrt(distance.max(initial=0)))
        if maximum>bound:
            raise ValueError('A static refined facet has no complete native source witness: '+sid+': '+str(maximum))
        if parent_faces_out is not None:
            parent_faces_out[sid]=np.full(len(faces),-1,int);parent_faces_out[sid][take]=parents
        restored,proof=source_normal_words(encoded_positions,{sid:fs},{sid:parents},old_p,native,old_words[sid])
        words=restored[sid]
        # Source corners retain their exact stored words, including chart seams.
        literal=source_positions[fs][:,:,None,:]==triangles[:,None,:,:]
        literal=np.all(literal,axis=3)
        for corner in range(3):
            actual,parent_corner=np.nonzero(literal[:,corner,:])
            words[actual,corner]=old_words[sid][parents[actual],parent_corner]
        tri=encoded_positions[fs];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        if np.any(np.einsum('tki,ti->tk',words.astype(float),normal)<=0):
            raise ValueError('Literal native static normal became inward: '+sid)
        fields[sid][take]=words
        evidence.append(dict(surface=sid,static_faces=len(take),maximum_source_corner_distance_m=maximum,
                             numerical_distance_bound_m=bound,source_field=proof))
    return fields,evidence


def restore_saved_fields(candidate_directory, output):
    """Checkpoint a saved candidate with corrected source-bound static fields."""
    folder,output=map(lambda p:Path(p).resolve(),(candidate_directory,output))
    if output.exists():raise FileExistsError(output)
    previous=json.loads((folder/'construction.json').read_text())
    if digest_file(previous['catalogue'])!=previous['catalogue_sha256']:
        raise ValueError('The installed source atlas changed')
    if digest_file(folder/'construction.npz')!=previous['source_geometry_sha256']:
        raise ValueError('The saved candidate changed')
    source=json.loads(Path(previous['canonical_delivery']).read_text())
    with np.load(folder/'construction.npz',allow_pickle=False)as z:arrays={key:z[key].copy()for key in z.files}
    surfaces={sid:arrays[key]for sid,key in previous['surfaces'].items()}
    moving=arrays['height']>0
    fields,evidence=fixed_source_fields(arrays['source_positions'],arrays['encoded_positions'],surfaces,moving,source)
    for sid,faces in surfaces.items():
        changed=np.any(moving[faces],axis=1)
        fields[sid][changed]=arrays[previous['corner_words'][sid]][changed]
        arrays[previous['corner_words'][sid]]=fields[sid]
    output.mkdir(parents=True)
    atomic_npz(output/'construction.npz',**arrays)
    plan=dict(previous,source_geometry=str(output/'construction.npz'),source_geometry_sha256=digest_file(output/'construction.npz'),
              resumed_construction_sha256=digest_file(folder/'construction.npz'),static_normal_source=evidence,
              geometric_positions_unchanged_from_saved_candidate=True,installed=False,anatomical_acceptance=False)
    atomic_json(output/'construction.json',plan)
    return dict(state='prepared',folder=str(output),static_surfaces=len(evidence),installed=False)


def pack_saved(candidate_directory):
    """Validate one saved construction; never rerun its geometry solve."""
    folder=Path(candidate_directory).resolve();plan=json.loads((folder/'construction.json').read_text())
    if digest_file(plan['catalogue'])!=plan['catalogue_sha256']or digest_file(folder/'construction.npz')!=plan['source_geometry_sha256']:
        raise ValueError('Source or candidate changed')
    source=json.loads(Path(plan['canonical_delivery']).read_text())
    with np.load(folder/'construction.npz',allow_pickle=False)as z:
        p=z['positions'].copy();surfaces={sid:z[key].copy()for sid,key in plan['surfaces'].items()};words={sid:z[key].copy()for sid,key in plan['corner_words'].items()}
    result=pack(p,surfaces,plan['cells'],{pid:pid for pid in plan['cells']},source['frame'],folder/'packed',
                origin=source['origins']['bladder'],corner_words=words,max_pairs=1000000)
    report=dict(**plan,packed=result,passed=result['passed'])
    atomic_json(folder/'candidate.json',report)
    return report


def ambient_relief_height(points, clinical_distance, covering_distance, *, center,
                          height_m=.00035, spacing_m=.004, covering_transition_m=.0012):
    """Separate clinical protection from finite wall-depth modulation."""
    covering_distance=np.asarray(covering_distance,float)
    if covering_distance.shape!=(len(points),)or np.any(covering_distance<0)or not np.isfinite(covering_distance).all():
        raise ValueError('Corresponding finite distances to the fixed outer band required')
    if not .001<=covering_transition_m<=.0035:
        raise ValueError('Explicit bounded tissue-band falloff required')
    return fold_height(points,clinical_distance,center=center,height_m=height_m,spacing_m=spacing_m,transition_m=.0035)*smoothstep(covering_distance/covering_transition_m)


def ambient_band_field(candidate_directory, source_carriers, output, *, phase, maximum_batches=2):
    """Evaluate a continuous physical-coordinate field, without source rays.

    Distance to the complete fixed band boundary is continuous even where a
    radial first hit switches to another source facet. This construction still
    needs actual discrete cell/fluid checks; no analytical acceptance is assumed.
    """
    folder,source_carriers,output=map(lambda p:Path(p).resolve(),(candidate_directory,source_carriers,output))
    if phase not in ['clinical','covering','geometry','normals']:
        raise ValueError('Choose one explicit saved ambient-field stage')
    previous=json.loads((folder/'construction.json').read_text())
    if digest_file(previous['catalogue'])!=previous['catalogue_sha256']or digest_file(previous['source_geometry'])!=previous['source_geometry_sha256']:
        raise ValueError('Source atlas or construction changed')
    source=json.loads(Path(previous['canonical_delivery']).read_text());guide=json.loads((source_carriers/'guide.json').read_text())
    if not guide['passed']or guide['source_catalogue_sha256']!=previous['catalogue_sha256']:
        raise ValueError('Current verified tissue-band guide required')
    binding=dict(construction=str(folder/'construction.json'),construction_sha256=digest_file(folder/'construction.json'),
                 guide_sha256=digest_file(source_carriers/'guide.json'),clinical_transition_m=.0035,covering_transition_m=.0012)
    output.mkdir(parents=True,exist_ok=True)
    if(output/'input-binding.json').exists():
        if json.loads((output/'input-binding.json').read_text())!=binding:raise ValueError('Saved ambient-field inputs changed')
    else:atomic_json(output/'input-binding.json',binding)
    band=set(guide['selected_part_ids']);cells=previous['cells']
    with np.load(previous['source_geometry'],allow_pickle=False)as z:
        p=z['source_positions'].copy();surfaces={sid:z[key].copy()for sid,key in previous['surfaces'].items()}
        hints={sid:z[key].copy()for sid,key in previous.get('normal_source_parents',{}).items()}
        seed={key:z[key].copy()for key in ['queried_controls','clinical_distance','covering_distance']if key in z}
    lining=set(previous['deformed_surfaces']+previous['protected_surfaces'])
    covering_surfaces={u['surface']for pid in band for u in source['cells'][pid]if u.get('neighbor')not in band and u['surface']not in lining}
    clinical_surfaces=set(previous['protected_surfaces']);fixed_surfaces=covering_surfaces|clinical_surfaces
    fixed=np.unique(np.concatenate([surfaces[sid]for sid in fixed_surfaces]))
    ids=np.unique(np.concatenate([surfaces[u['surface']]for pid in band for u in cells[pid]]))
    ids=np.setdiff1d(ids[p[ids,1]>previous['neck_fixed_below_y_m']],fixed)
    if previous.get('ambient_boundary_reference'):
        reference=Path(previous['ambient_boundary_reference'])
        if digest_file(reference)!=previous['ambient_boundary_reference_sha256']:raise ValueError('The continuous field reference changed')
        ref=json.loads(reference.read_text())
        if digest_file(ref['source_geometry'])!=ref['source_geometry_sha256']:raise ValueError('The reference boundary geometry changed')
        with np.load(ref['source_geometry'],allow_pickle=False)as z:
            refp=z['source_positions'].copy();ref_surfaces={sid:z[key].copy()for sid,key in ref['surfaces'].items()}
        protected=Mesh(refp,np.concatenate([ref_surfaces[sid]for sid in clinical_surfaces]),source['frame'])
        covering=Mesh(refp,np.concatenate([ref_surfaces[sid]for sid in covering_surfaces]),source['frame'])
        if seed:
            seed_path=Path(previous['saved_query_seed'])
            if digest_file(seed_path)!=previous['saved_query_seed_sha256']:raise ValueError('Saved query seed changed')
            seed_meta=json.loads(seed_path.read_text())
            if digest_file(seed_meta['source_geometry'])!=seed_meta['source_geometry_sha256']:raise ValueError('Saved query geometry changed')
            with np.load(seed_meta['source_geometry'],allow_pickle=False)as z:
                if any(not np.array_equal(seed[key],z[key])for key in seed):raise ValueError('A copied query seed differs from its source')
                if not np.array_equal(p[seed['queried_controls']],z['source_positions'][seed['queried_controls']]):
                    raise ValueError('A reused query control changed position')
    else:
        protected=Mesh(p,np.concatenate([surfaces[sid]for sid in clinical_surfaces]),source['frame'])
        covering=Mesh(p,np.concatenate([surfaces[sid]for sid in covering_surfaces]),source['frame'])
        seed={}
    if phase in ['clinical','covering']:
        mesh=protected if phase=='clinical'else covering
        path=output/(phase+'-query.npz');proof=output/(phase+'-query.json')
        if path.exists()or proof.exists():raise FileExistsError('Reuse the saved query instead of repeating it')
        if not 1<=maximum_batches<=3:raise ValueError('Each query stage may calculate at most three saved batches')
        tree=BVH(mesh);saved=[];done=0;batch_size=8192
        for start in range(0,len(ids),batch_size):
            batch=output/(phase+'-batch-'+str(start)+'.npz');batch_proof=batch.with_suffix('.json');stop=min(start+batch_size,len(ids))
            if batch.exists()and batch_proof.exists():
                meta=json.loads(batch_proof.read_text())
                if meta['binding']!=binding or digest_file(batch)!=meta['query_sha256']:raise ValueError('Saved distance batch changed')
            else:
                if done>=maximum_batches:
                    return dict(state='distance_batches_saved',phase=phase,saved_batches=len(saved),total_batches=(len(ids)+batch_size-1)//batch_size,installed=False)
                selected=ids[start:stop];distances=np.empty(len(selected));reuse=np.zeros(len(selected),bool)
                if seed:
                    seed_lookup={int(k):i for i,k in enumerate(seed['queried_controls'])}
                    matched=np.array([seed_lookup.get(int(k),-1)for k in selected]);reuse=matched>=0
                    distances[reuse]=seed[phase+'_distance'][matched[reuse]]
                evidence=dict(queries=0,nearest_triangle_verified=True,float64_source_witnesses=True)
                if np.any(~reuse):
                    query,evidence=tree.closest_points(p[selected[~reuse]],backend='auto');distances[~reuse]=query['distance_m']
                evidence.update(reused_source_bound_queries=int(reuse.sum()),new_queries=int((~reuse).sum()),
                                original_field_reference=previous.get('ambient_boundary_reference'),
                                saved_query_seed=previous.get('saved_query_seed'),saved_query_seed_sha256=previous.get('saved_query_seed_sha256'))
                atomic_npz(batch,ids=selected,distance_m=distances)
                meta=dict(evidence=evidence,query_sha256=digest_file(batch),binding=binding,passed=True)
                atomic_json(batch_proof,meta);done+=int(np.any(~reuse))
                print(json.dumps(dict(phase=phase,completed_query_controls=stop,total_query_controls=len(ids))),flush=True)
            with np.load(batch,allow_pickle=False)as z:
                if not np.array_equal(z['ids'],ids[start:stop]):raise ValueError('Saved distance batch has different controls')
                saved.append(z['distance_m'].copy())
        evidence=dict(method='source-bound independently certified distance batches',queries=len(ids),batch_size=batch_size,
                      nearest_triangle_verified=True,float64_source_witnesses=True,anatomical_acceptance=False)
        atomic_npz(path,ids=ids,distance_m=np.concatenate(saved))
        atomic_json(proof,dict(evidence=evidence,query_sha256=digest_file(path),binding=binding,passed=True))
        return dict(state='distance_saved',phase=phase,queries=len(ids),installed=False)
    distances={};proofs={}
    for role in ['clinical','covering']:
        path=output/(role+'-query.npz');proof=json.loads((output/(role+'-query.json')).read_text())
        if proof['binding']!=binding or digest_file(path)!=proof['query_sha256']:raise ValueError('A saved distance query changed')
        with np.load(path,allow_pickle=False)as z:
            if not np.array_equal(z['ids'],ids):raise ValueError('Saved distances refer to different controls')
            distances[role]=z['distance_m'].copy()
        proofs[role]=proof['evidence']
    query={'distance_m':distances['clinical']};evidence=proofs['clinical']
    outer_query={'distance_m':distances['covering']};outer_evidence=proofs['covering']
    height=np.zeros(len(p));height[ids]=ambient_relief_height(p[ids],query['distance_m'],outer_query['distance_m'],center=previous['center_m'],
        height_m=previous['maximum_height_m'],spacing_m=previous['nominal_spacing_m'])
    direction=np.asarray(previous['center_m'])-p;direction/=np.maximum(np.linalg.norm(direction,axis=1)[:,None],1e-20)
    moved=p+height[:,None]*direction;origin=np.asarray(source['origins']['bladder'],float)
    encoded=(moved-origin).astype('<f4').astype(float)+origin;moving=height>0
    geometry=dict(source_positions=p,positions=moved,encoded_positions=encoded,height=height,fixed_controls=fixed,
                  queried_controls=ids,clinical_distance=query['distance_m'],covering_distance=outer_query['distance_m'])
    if phase=='geometry':
        atomic_npz(output/'field.npz',**geometry)
        atomic_json(output/'field.json',dict(binding=binding,field_sha256=digest_file(output/'field.npz'),
                    maximum_height_m=float(height.max()),moved_controls=int(moving.sum()),installed=False))
        return dict(state='field_saved',phase=phase,moved_controls=int(moving.sum()),maximum_height_m=float(height.max()),installed=False)
    field_proof=json.loads((output/'field.json').read_text())
    if field_proof['binding']!=binding or digest_file(output/'field.npz')!=field_proof['field_sha256']:
        raise ValueError('Saved field changed before normal generation')
    with np.load(output/'field.npz',allow_pickle=False)as z:
        if any(not np.array_equal(z[key],value)for key,value in geometry.items()):raise ValueError('Saved field does not replay')
    affected={sid for sid,fs in surfaces.items()if np.any(moving[fs])}
    generated,normal_evidence=grouped_normal_words(encoded,{sid:surfaces[sid]for sid in affected},
        {sid:np.full(len(surfaces[sid]),'fluid-lining'if sid in lining else sid,dtype=object)for sid in affected})
    parents={};fields,static=fixed_source_fields(p,encoded,surfaces,moving,source,parent_hints=hints,parent_faces_out=parents)
    for sid,words in generated.items():
        changed=np.any(moving[surfaces[sid]],axis=1);fields[sid][changed]=words[changed]
    arrays=geometry
    keys={};word_keys={};parent_keys={}
    for i,(sid,fs)in enumerate(surfaces.items()):
        keys[sid]='surface'+str(i);arrays[keys[sid]]=fs;word_keys[sid]='words'+str(i);arrays[word_keys[sid]]=fields[sid]
        if sid in parents:parent_keys[sid]='native_parent'+str(i);arrays[parent_keys[sid]]=parents[sid]
    atomic_npz(output/'construction.npz',**arrays)
    plan=dict(previous,schema='bladder-lining-ambient-relief/v1',source_geometry=str(output/'construction.npz'),
        source_geometry_sha256=digest_file(output/'construction.npz'),surfaces=keys,corner_words=word_keys,normal_source_parents=parent_keys,
        field_method='continuous physical-coordinate ridge field with separate C2 clinical and outer-band distance gates',
        radial_source_queries=[],source_ray_monotonicity_bound=None,ambient_distance_query=evidence,outer_band_distance_query=outer_evidence,
        clinical_transition_m=.0035,covering_transition_m=.0012,
        static_normal_source=static,normal_field=normal_evidence,source_ray_bound_is_not_a_discrete_mesh_proof=True,
        moved_controls=int(moving.sum()),actual_maximum_displacement_m=float(height.max()),installed=False,anatomical_acceptance=False)
    atomic_json(output/'construction.json',plan)
    return dict(state='prepared',moved_controls=plan['moved_controls'],maximum_height_m=plan['actual_maximum_displacement_m'],installed=False)


def construct(delivery_directory, catalogue, output, *, expected_catalogue_sha256,
              maximum_edge_m=.0008, height_m=.00035, spacing_m=.004):
    folder, catalogue, output = map(lambda p: Path(p).resolve(), (delivery_directory, catalogue, output))
    if output.exists():
        raise FileExistsError(output)
    if digest_file(catalogue) != expected_catalogue_sha256:
        raise ValueError('The installed bladder revision changed')
    if not .0006 <= maximum_edge_m <= .0015:
        raise ValueError('Bounded local lining refinement spacing required')
    source = json.loads((folder/'delivery.json').read_text())
    atlas = json.loads(catalogue.read_text())
    records = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    for pid in source['groups']['bladder']['parts']:
        actual = reader.read_part(atlas, records[pid], frame=source['frame'])
        original = Mesh.load(folder/'bladder'/(pid.replace(':', '_')+'.npz'))
        if actual.geometry_sha != original.geometry_sha:
            raise ValueError('The canonical tissue source is not the installed cell: '+pid)
    with np.load(source['source'], allow_pickle=False) as z:
        original_positions = z['encoded_positions'].copy()
        old_surfaces = {sid: z[key].copy() for sid, key in source['surfaces'].items()}
        old_words = {sid: z[key].copy() for sid, key in source['corner_words'].items()}
    owners = {sid: [] for sid in old_surfaces}
    for pid, uses in source['cells'].items():
        for use in uses:
            owners[use['surface']].append(pid)
    free_lining = {sid for sid, pids in owners.items() if len(pids) == 1 and records[pids[0]].get('tissueWallLayer', '').startswith('mucosa')
                   and not next(u for u in source['cells'][pids[0]] if u['surface'] == sid).get('neighbor')}
    target = {sid for sid in free_lining if owners[sid][0] in ['BLADDER:male:body:mucosa', 'BLADDER:male:fundus:mucosa']}
    protected_surfaces = free_lining-target
    if len(target) != 2 or not protected_surfaces:
        raise ValueError('Explicit native body, fundus and protected trigonal lining are required')
    output.mkdir(parents=True)
    center = np.array([0., .919, .002])
    plan = dict(schema='bladder-lining-relief/v1', catalogue=str(catalogue), catalogue_sha256=expected_catalogue_sha256,
                canonical_delivery=str(folder/'delivery.json'), canonical_delivery_sha256=digest_file(folder/'delivery.json'),
                maximum_edge_m=maximum_edge_m, maximum_height_m=height_m, nominal_spacing_m=spacing_m,
                fixed_support_falloff_m=.004, neck_fixed_below_y_m=.899, center_m=center.tolist(),
                deformed_surfaces=sorted(target), protected_surfaces=sorted(protected_surfaces),
                dimensions_are_illustrative=True, original_outer_wall_fixed=True, anatomical_acceptance=False, installed=False)
    atomic_json(output/'plan.json', plan)
    bank = CommonPositionBank(original_positions)
    surfaces = {sid: f.copy() for sid, f in old_surfaces.items()}
    ancestry = {sid: np.arange(len(f)) for sid, f in surfaces.items()}
    rounds = []
    for iteration in range(10):
        p = bank.positions
        selected = np.concatenate([surfaces[sid] for sid in target])
        selected = selected[p[selected, 1].max(1) > .899]
        marked = set()
        for a, b in [(0, 1), (1, 2), (2, 0)]:
            e = selected[:, [a, b]]
            long = np.linalg.norm(p[e[:, 0]]-p[e[:, 1]], axis=1) > maximum_edge_m
            marked.update(tuple(sorted(map(int, edge))) for edge in e[long])
        if not marked:
            break
        midpoints = {edge: bank.insert((p[edge[0]]+p[edge[1]])*.5) for edge in sorted(marked)}
        surfaces, parents, changed = bisect_surfaces(surfaces, midpoints)
        ancestry = {sid: ancestry[sid][parents[sid]] for sid in surfaces}
        rounds.append(dict(iteration=iteration, split_edges=len(marked), changed_triangles=changed))
        if sum(len(surfaces[sid]) for sid in target) > 100000:
            atomic_json(output/'refinement-rejection.json',dict(passed=False,reason='Local lining refinement exceeded the fixed 100,000-triangle budget',
                observed_lining_triangles=sum(len(surfaces[sid])for sid in target),rounds=rounds,installed=False))
            raise ValueError('Local lining refinement exceeded the explicit 100,000-triangle budget')
    else:
        raise ValueError('Local lining refinement did not converge')
    p = bank.positions
    lining_ids = np.unique(np.concatenate([surfaces[sid] for sid in target]))
    fixed = np.unique(np.concatenate([surfaces[sid] for sid in protected_surfaces]))
    fixed_mesh = Mesh(p, np.concatenate([surfaces[sid] for sid in protected_surfaces]), source['frame'])
    query, query_evidence = BVH(fixed_mesh).closest_points(p[lining_ids], backend='auto')
    height = np.zeros(len(p))
    height[lining_ids] = fold_height(p[lining_ids], query['distance_m'], center=center, height_m=height_m, spacing_m=spacing_m)
    height[fixed] = 0.
    direction = center-p
    direction /= np.maximum(np.linalg.norm(direction, axis=1)[:, None], 1e-20)
    moved = p+direction*height[:, None]
    moving = height > 0
    if not moving.any() or not np.array_equal(moved[fixed], p[fixed]):
        raise ValueError('The lining relief must be nonempty and retain every protected source control')
    origin = np.asarray(source['origins']['bladder'], float)
    encoded = (moved-origin).astype('<f4').astype(float)+origin
    # Every original nonlining point remains in its existing native chart.
    # Only the bladder chart is repacked; other delivery banks are reused.
    used = {u['surface'] for pid, uses in source['cells'].items() if source['chart_for'][pid] == 'bladder' for u in uses}
    fields = {}
    for sid in used:
        fs = surfaces[sid]
        src = old_surfaces[sid][ancestry[sid]]
        weights = np.array([barycentric(t, parent) for t, parent in zip(p[fs], original_positions[src])])
        values = np.einsum('tck,tkj->tcj', weights, old_words[sid][ancestry[sid]].astype(float))
        values /= np.linalg.norm(values, axis=2)[:, :, None]
        values = np.rint(values*32767).astype('<i2')
        for i, (actual, parent) in enumerate(zip(p[fs], original_positions[src])):
            for j, point in enumerate(actual):
                literal = np.flatnonzero(np.all(parent == point, axis=1))
                if len(literal):
                    values[i, j] = old_words[sid][ancestry[sid][i], literal[0]]
        fields[sid] = values
    affected = {sid for sid in used if np.any(moving[surfaces[sid]])}
    geometric, normal_evidence = grouped_normal_words(encoded, {sid: surfaces[sid] for sid in affected},
        {sid: np.full(len(surfaces[sid]), 'fluid-lining' if sid in free_lining else sid, dtype=object) for sid in affected})
    for sid, values in geometric.items():
        changed_faces = np.any(moving[surfaces[sid]], axis=1)
        fields[sid][changed_faces] = values[changed_faces]
    cells = {pid: [dict(u, neighbor=u.get('neighbor') if u.get('neighbor') in source['chart_for'] and source['chart_for'][u['neighbor']] == 'bladder' else None) for u in uses]
             for pid, uses in source['cells'].items() if source['chart_for'][pid] == 'bladder'}
    arrays, keys, word_keys = dict(source_positions=p, positions=moved, encoded_positions=encoded, height=height, fixed_controls=fixed), {}, {}
    outward_failures=[]
    for i, sid in enumerate(used):
        keys[sid] = 'surface'+str(i); arrays[keys[sid]] = surfaces[sid]
        word_keys[sid]='words'+str(i);arrays[word_keys[sid]]=fields[sid]
        tri=encoded[surfaces[sid]];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        bad=np.flatnonzero(np.any(np.einsum('tki,ti->tk',fields[sid].astype(float),normal)<=0,axis=1))
        if len(bad):outward_failures.append(dict(surface=sid,faces=bad.tolist(),zero_area_faces=int(np.sum(np.linalg.norm(normal[bad],axis=1)==0))))
    atomic_npz(output/'construction.npz', **arrays)
    construction=dict(**plan,source_geometry=str(output/'construction.npz'),source_geometry_sha256=digest_file(output/'construction.npz'),
        surfaces=keys,corner_words=word_keys,cells=cells,refinement=rounds,moved_controls=int(moving.sum()),
        actual_maximum_displacement_m=float(height.max()),outward_failures=outward_failures,
        refined_source_controls=len(p),original_source_controls=len(original_positions),distance_query=query_evidence,normal_field=normal_evidence)
    atomic_json(output/'construction.json',construction)
    if outward_failures:raise ValueError('Actual emitted normals/areas fail; inspect construction.json')
    result = pack(moved, {sid: surfaces[sid] for sid in used}, cells, {pid: pid for pid in cells}, source['frame'], output/'packed',
                  origin=origin, corner_words=fields, max_pairs=1000000)
    report = dict(**plan, source_geometry=str(output/'construction.npz'), source_geometry_sha256=digest_file(output/'construction.npz'),
                  surfaces=keys, refined_source_controls=len(p), original_source_controls=len(original_positions), moved_controls=int(moving.sum()),
                  actual_maximum_displacement_m=float(height.max()), refinement=rounds, distance_query=query_evidence, normal_field=normal_evidence,
                  protected_lining_controls_literal=True, source_other_charts_reused=True, packed=result, passed=result['passed'])
    atomic_json(output/'candidate.json', report)
    return report


def verify_fluid(candidate_directory, source_fluid, output):
    """Replace the literal source lining in the existing complete fluid cell.

    The original computational terminal disks are retained for this volume
    proof only. They are not returned as visible anatomy or as tissue caps.
    """
    from .audit import audit
    from .algorithms.finite_tissue_boolean import make_solid
    from .algorithms.constrained_boundary_disk import cancel_opposed_faces
    folder, source_fluid, output = map(lambda p: Path(p).resolve(), (candidate_directory, source_fluid, output))
    if output.exists():
        raise FileExistsError(output)
    candidate = json.loads((folder/'construction.json').read_text())
    if digest_file(folder/'construction.npz') != candidate['source_geometry_sha256']:
        raise ValueError('The saved lining candidate changed')
    source = json.loads(Path(candidate['canonical_delivery']).read_text())
    with np.load(source['source'], allow_pickle=False) as z:
        old_p = z['encoded_positions'].copy()
        names = candidate['deformed_surfaces']+candidate['protected_surfaces']
        source_lining = np.concatenate([old_p[z[source['surfaces'][sid]]] for sid in names])
    old = Mesh.load(source_fluid)
    key = lambda t: min(tuple(map(tuple, np.roll(t, k, axis=0))) for k in range(3))
    expected = {key(t[[0, 2, 1]]) for t in source_lining}
    actual = {key(t) for t in old.triangles}
    missing = expected-actual
    if missing:
        raise ValueError('The complete fluid source does not contain every literal source lining facet: '+str(len(missing)))
    with np.load(folder/'construction.npz', allow_pickle=False) as z:
        p = z['encoded_positions'].copy()
        lining = np.concatenate([z[candidate['surfaces'][sid]][:, [0, 2, 1]] for sid in names])
    bank = CommonPositionBank(p)
    old_map = np.array([bank.insert(v) for v in old.vertices], int)
    retained = old_map[old.faces[[i for i, t in enumerate(old.triangles) if key(t) not in expected]]]
    faces, cancellation = cancel_opposed_faces(np.r_[retained, lining])
    result = Mesh(bank.positions, faces, old.frame)
    check = audit(result, require_closed=True, expected_components=1, expected_genus=[0], max_pairs=1000000)
    output.mkdir(parents=True)
    result.save(output/'fluid.npz')
    proof = dict(source_fluid=str(source_fluid), source_fluid_sha256=digest_file(source_fluid),
                 construction_sha256=digest_file(folder/'construction.npz'), all_source_lining_facets_replaced=len(expected),
                 original_nonlining_triangles_retained=len(retained), check=check, cancellation=cancellation,
                 computational_terminal_disks_not_anatomy=True, anatomical_acceptance=False, passed=check['passed'])
    if check['passed']:
        before, _ = make_solid(old.vertices, old.faces)
        after, _ = make_solid(result.vertices, result.faces)
        outside = float((after-before).volume())
        # Only source Float32 storage changes can protrude beyond the old
        # clinical cavity; the constructed relief itself points into it.
        extent = float(np.ptp(old.bounds, axis=0).max())
        bound = 2048*np.finfo(float).eps*extent**3*(len(old.faces)+len(result.faces))
        proof.update(source_volume_m3=float(before.volume()), current_volume_m3=float(after.volume()),
                     outside_source_fluid_m3=outside, numerical_volume_bound_m3=bound,
                     passed=bool(check['passed'] and after.volume()>0 and after.volume()<before.volume() and outside<=bound))
    atomic_json(output/'verification.json', proof)
    return proof


def refine_failed_cell_source(candidate_directory, output):
    """Bisect only measured crossing facets, across every incident tissue.

    New controls lie on the undeformed source edges. They must subsequently
    receive the same evaluated wall field; merely splitting deformed triangles
    would preserve the defect and is deliberately not emitted as a candidate.
    """
    from .contact import pair_contacts
    folder,output=map(lambda p:Path(p).resolve(),(candidate_directory,output))
    if output.exists():raise FileExistsError(output)
    plan=json.loads((folder/'construction.json').read_text())
    if digest_file(plan['catalogue'])!=plan['catalogue_sha256']or digest_file(plan['source_geometry'])!=plan['source_geometry_sha256']:
        raise ValueError('Current source or saved candidate changed')
    staged=json.loads((folder/'staged/tissue.json').read_text())
    if digest_file(staged['bank'])!=staged['bank_sha256']:raise ValueError('The emitted tissue bank changed')
    failed=[name for name,check in staged['checks']['cells'].items()if not check['passed']]
    if not failed:raise ValueError('Actual measured cell defects are required')
    with np.load(plan['source_geometry'],allow_pickle=False)as z:
        p=z['source_positions'].copy();encoded=z['encoded_positions'].copy()
        surfaces={sid:z[key].copy()for sid,key in plan['surfaces'].items()}
        hints={sid:z[key].copy()for sid,key in plan.get('normal_source_parents',{}).items()}
        seed={key:z[key].copy()for key in ['queried_controls','clinical_distance','covering_distance']if key in z}
    lookup={tuple(v):i for i,v in enumerate(encoded)};marked=set();witnesses=[]
    cache=folder/'refinement-witnesses';cache.mkdir(exist_ok=True)
    for name in sorted(failed):
        row=staged['mesh_files'][name]
        if digest_file(row['path'])!=row['sha256']:raise ValueError('An emitted source cell changed')
        mesh=Mesh.load(row['path']);path=cache/(name.replace(':','_')+'.json')
        if path.exists():
            saved=json.loads(path.read_text())
            if saved['mesh_sha256']!=row['sha256']:raise ValueError('Stale crossing witnesses')
            bad=np.asarray(saved['pairs'],int)
        else:
            check=staged['checks']['cells'][name];summary=check.get('self_intersections',{})
            examples=summary.get('examples',[])
            reusable=(check.get('geometry_sha256')==mesh.geometry_sha and check.get('attributes_sha256')==mesh.attributes_sha
                      and not summary.get('truncated',True)and summary.get('defect_pairs')==len(examples)and len(examples)>0)
            if reusable:bad=np.asarray(examples,int)
            else:
                pairs,summary=pair_contacts(mesh,mesh,0.,selfcheck=True,max_pairs=1000000)
                if summary['truncated']:raise ValueError('Crossing witness search exceeded its fixed pair budget')
                bad=pairs['pairs'][np.isin(pairs['relation'],['crossing','coplanar-overlap'])]
            saved=dict(cell=name,mesh_sha256=row['sha256'],pairs=bad.tolist(),summary=summary,complete_bound_evidence_reused=reusable)
            atomic_json(path,saved)
        if not len(bad):raise ValueError('The actual emitted cell no longer reproduces its defect')
        for tri in mesh.triangles[np.unique(bad)]:
            ids=[lookup[tuple(v)]for v in tri]
            marked.update(tuple(sorted((ids[a],ids[b])))for a,b in[(0,1),(1,2),(2,0)])
        witnesses.append(dict(cell=name,defect_pairs=len(bad),source_mesh_sha256=row['sha256']))
        print(json.dumps(witnesses[-1]),flush=True)
    if len(marked)>2000:raise ValueError('Adaptive local refinement exceeds its explicit 2,000-edge ceiling')
    bank=CommonPositionBank(p)
    midpoints={edge:bank.insert((p[edge[0]]+p[edge[1]])*.5)for edge in sorted(marked)}
    refined,parents,changed=bisect_surfaces(surfaces,midpoints)
    if sum(map(len,refined.values()))>sum(map(len,surfaces.values()))+10000:
        raise ValueError('Adaptive local refinement exceeds its 10,000-triangle increment ceiling')
    output.mkdir(parents=True);arrays=dict(source_positions=bank.positions,**seed);keys={};parent_keys={}
    for i,(sid,faces)in enumerate(refined.items()):
        keys[sid]='surface'+str(i);arrays[keys[sid]]=faces
        if sid in hints:
            parent_keys[sid]='native_parent'+str(i);arrays[parent_keys[sid]]=hints[sid][parents[sid]]
    atomic_npz(output/'construction.npz',**arrays)
    result=dict(plan,source_geometry=str(output/'construction.npz'),source_geometry_sha256=digest_file(output/'construction.npz'),
                surfaces=keys,corner_words={},normal_source_parents=parent_keys,refined_from_saved_geometry_sha256=plan['source_geometry_sha256'],
                refinement_only_source=True,adaptive_witnesses=witnesses,adaptive_split_edges=len(marked),
                adaptive_changed_triangles=changed,source_cell_checks_complete=staged['checks']['all_emitted_cells_checked'],
                installed=False,anatomical_acceptance=False)
    if seed:
        reference=plan.get('ambient_boundary_reference',str(folder/'construction.json'))
        result.update(ambient_boundary_reference=reference,ambient_boundary_reference_sha256=digest_file(reference),
                      saved_query_seed=str(folder/'construction.json'),saved_query_seed_sha256=digest_file(folder/'construction.json'))
    atomic_json(output/'construction.json',result)
    return dict(state='source_refined',split_edges=len(marked),changed_triangles=changed,witnesses=witnesses,installed=False)


def propagate_through_wall(candidate_directory, source_carriers, output, *, validate=True):
    """Carry a resolved radial field through the mucosal/submucosal band."""
    from .open3d_backend import SurfaceIndex
    from .distance import barycentric as triangle_weights
    folder, source_carriers, output = map(lambda p:Path(p).resolve(),(candidate_directory,source_carriers,output))
    if output.exists():
        raise FileExistsError(output)
    previous = json.loads((folder/'construction.json').read_text())
    if digest_file(previous['catalogue']) != previous['catalogue_sha256']:
        raise ValueError('The source atlas changed')
    if digest_file(folder/'construction.npz')!=previous['source_geometry_sha256']:
        raise ValueError('The saved source construction changed')
    source = json.loads(Path(previous['canonical_delivery']).read_text())
    guide=json.loads((source_carriers/'guide.json').read_text())
    if not guide['passed']or guide['source_catalogue_sha256']!=previous['catalogue_sha256']:
        raise ValueError('A current verified finite inner-wall guide is required')
    band=set(guide['selected_part_ids'])
    atlas = json.loads(Path(previous['catalogue']).read_text())
    records = {p['id']: p for p in atlas['parts']}
    with np.load(folder/'construction.npz', allow_pickle=False) as z:
        p = z['source_positions'].copy()
        surfaces = {sid: z[key].copy() for sid, key in previous['surfaces'].items()}
        parent_hints={sid:z[key].copy()for sid,key in previous.get('normal_source_parents',{}).items()}
    with np.load(source['source'], allow_pickle=False) as z:
        old_p = z['encoded_positions'].copy()
        old_surfaces = {sid: z[key].copy() for sid, key in source['surfaces'].items()}
    cells = previous['cells']
    lining = previous['deformed_surfaces']+previous['protected_surfaces']
    source_p=p.copy();source_surfaces={sid:f.copy()for sid,f in surfaces.items()};ancestry={sid:np.arange(len(f))for sid,f in surfaces.items()}
    bank=CommonPositionBank(p);refinement=[]
    eligible_surfaces={u['surface']for pid in band for u in cells[pid]}
    fixed_edges=set()
    for pid in band:
        for use in source['cells'][pid]:
            peer=use.get('neighbor')
            if peer and source['chart_for'].get(peer)!='bladder':
                for face in surfaces[use['surface']]:
                    fixed_edges.update(tuple(sorted(map(int,e)))for e in zip(face,np.roll(face,-1)))
    fixed_profile_controls=np.unique(np.asarray(list(fixed_edges),int))if fixed_edges else np.empty(0,int)
    output.mkdir(parents=True)
    atomic_json(output/'refinement-plan.json',dict(source_construction_sha256=digest_file(folder/'construction.npz'),source_guide_sha256=digest_file(source_carriers/'guide.json'),
        maximum_edge_m=previous['maximum_edge_m'],profile_controls_literal=True,fixed_profile_controls=fixed_profile_controls.tolist(),
        incident_profile_fans_retained=True,maximum_canonical_triangles=500000,installed=False))
    for iteration in range(10):
        p=bank.positions;faces=np.concatenate([surfaces[sid]for sid in eligible_surfaces]);faces=faces[p[faces,1].max(1)>.899];marked=set()
        # A frozen long profile edge cannot satisfy a smaller local target.
        # Keep its entire incident source fan, instead of infinitely bisecting
        # new diagonals towards a boundary that is explicitly immutable.
        faces=faces[~np.any(np.isin(faces,fixed_profile_controls),axis=1)]
        for a,b in [(0,1),(1,2),(2,0)]:
            edges=faces[:,[a,b]];long=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)>previous['maximum_edge_m']
            marked.update(tuple(sorted(map(int,e)))for e in edges[long])
        marked-=fixed_edges
        if not marked:break
        midpoints={e:bank.insert((p[e[0]]+p[e[1]])*.5)for e in sorted(marked)}
        surfaces,parents,changed=bisect_surfaces(surfaces,midpoints)
        parent_hints={sid:hints[parents[sid]]for sid,hints in parent_hints.items()}
        ancestry={sid:ancestry[sid][parents[sid]]for sid in surfaces}
        refinement.append(dict(iteration=iteration,split_edges=len(marked),changed_triangles=changed))
        if sum(map(len,surfaces.values()))>500000:
            atomic_json(output/'refinement-failure.json',dict(passed=False,rounds=refinement,canonical_triangles=sum(map(len,surfaces.values())),reason='Fixed triangle budget exceeded',installed=False))
            raise ValueError('Resolved inner-wall candidate exceeds the fixed 500,000 canonical-triangle budget')
    else:
        atomic_json(output/'refinement-failure.json',dict(passed=False,rounds=refinement,reason='Refinement did not converge within ten rounds',installed=False))
        raise ValueError('Inner-wall refinement did not converge')
    p=bank.positions
    owners={sid:[]for sid in surfaces}
    for pid,uses in cells.items():
        for use in uses:owners[use['surface']].append(pid)
    outside={u['surface']for pid in band for u in source['cells'][pid]if u.get('neighbor')not in band and u['surface']not in lining}
    inner = Mesh.load(source_carriers/'fluid.npz')
    outer = Mesh.load(source_carriers/'outer.npz')
    if inner.frame!=source['frame']or outer.frame!=source['frame']:raise ValueError('Source wall and fluid carriers occupy another frame')
    protected = Mesh(old_p, np.concatenate([old_surfaces[s] for s in previous['protected_surfaces']]), source['frame'])
    center = np.asarray(previous['center_m'], float)
    ids = np.unique(np.concatenate([surfaces[u['surface']]for pid in band for u in cells[pid]]))
    ids = ids[p[ids, 1] > previous['neck_fixed_below_y_m']]
    fixed = np.unique(np.concatenate([surfaces[s] for s in outside|set(previous['protected_surfaces'])]))
    ids = np.setdiff1d(ids, fixed)
    direction = p[ids]-center
    radius = np.linalg.norm(direction, axis=1)
    direction /= radius[:, None]
    angular_gate=np.ones(len(ids));opening_cones=[]
    for group in atlas['finiteUrinaryIntegration']['termGroups']:
        for boundary in group.get('openingBoundaries',[]):
            ring=np.asarray(boundary['points'])-center;length=np.linalg.norm(ring,axis=1);unit=ring/length[:,None]
            axis=unit.mean(0);axis/=np.linalg.norm(axis)
            angle=np.arccos(np.clip(direction@axis,-1,1));extent=float(np.arccos(np.clip(unit@axis,-1,1)).max())
            angular_gate*=smoothstep((angle-extent)/(previous['fixed_support_falloff_m']/length.min()))
            opening_cones.append(dict(name=boundary['name'],axis=axis.tolist(),source_ring_angular_extent_rad=extent))
    active=angular_gate>0;angular_fixed=ids[~active]
    ids,direction,radius,angular_gate=ids[active],direction[active],radius[active],angular_gate[active]
    ray_evidence = []
    def radii(mesh):
        source_ids=np.arange(len(mesh.faces));omitted=[]
        try:index = SurfaceIndex(mesh)
        except ValueError as error:
            if 'Float32 normalization collapses a source triangle'not in str(error):raise
            origin=mesh.bounds.mean(0);scale=float(np.ptp(mesh.bounds,axis=0).max())
            local=((mesh.vertices-origin)/scale).astype('<f4').astype(float);tri=local[mesh.faces]
            keep=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)>0
            source_ids=np.flatnonzero(keep);omitted=np.flatnonzero(~keep).tolist()
            # This is a native proposal copy only. Every original source facet
            # remains in the Float64 fallback and the independent geometry.
            index=SurfaceIndex(Mesh(mesh.vertices,mesh.faces[keep],mesh.frame))
        hits, evidence = index.cast_rays(np.repeat(center[None], len(ids), axis=0), direction)
        proposed=hits['face_ids'].copy();proposed[~hits['hit']]=0
        proposed=source_ids[proposed]
        tri = mesh.triangles[proposed]
        n = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        denominator = np.einsum('ij,ij->i', n, direction)
        t = np.divide(np.einsum('ij,ij->i',tri[:,0]-center,n),denominator,out=np.full(len(ids),np.nan),where=denominator!=0)
        points = center+t[:, None]*direction
        weights, valid = triangle_weights(points, tri)
        # Explicit source-plane intersections are independent of native hit
        # distances. Near-edge native candidates retain their measured residual
        # and never become an anatomical or first-hit certificate.
        overshoot = np.maximum(-weights.min(1), 0.)
        retry=np.flatnonzero(~hits['hit']|~valid|~np.isfinite(t)|(t<=0)|(overshoot>2e-5))
        triangles=mesh.triangles;normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);numerator=np.einsum('ij,ij->i',triangles[:,0]-center,normals)
        for query_id in retry:
            denominator=normals@direction[query_id]
            travel=np.divide(numerator,denominator,out=np.full(len(triangles),np.inf),where=denominator!=0)
            take=np.flatnonzero((travel>0)&np.isfinite(travel))
            point=center+travel[take,None]*direction[query_id]
            weight,good=triangle_weights(point,triangles[take])
            good&=np.all(weight>=-1e-12,axis=1)
            if not good.any():raise ValueError('A Float64 source ray has no hit outside the explicit protected opening cones')
            choices=take[good];hit=int(choices[np.argmin(travel[choices])]);t[query_id]=travel[hit]
            points[query_id]=center+t[query_id]*direction[query_id];overshoot[query_id]=float(max(0.,-weight[good][np.argmin(travel[choices])].min()))
        evidence.update(maximum_source_barycentric_overshoot=float(overshoot.max()),
                        exact_source_plane_intersections=True,complete_float64_fallback_queries=len(retry),nearest_hit_verified=False,
                        native_proposal_omitted_collapsed_source_faces=omitted,all_float64_source_faces_retained=True,
                        source_triangle_order_preserved=not omitted)
        ray_evidence.append(evidence)
        return t, points
    inner_radius, foot = radii(inner)
    outer_radius, _ = radii(outer)
    width = outer_radius-inner_radius
    if np.any(width <= 0):
        raise ValueError('Source ray does not resolve a positive whole-wall thickness')
    distance, distance_evidence = BVH(protected).closest_points(foot, backend='auto')
    amplitude = fold_height(foot, distance['distance_m'], center=center, height_m=previous['maximum_height_m'], spacing_m=previous['nominal_spacing_m'])*angular_gate
    if np.any(amplitude >= inner_radius):
        raise ValueError('The lining relief would collapse the radial fluid interval')
    depth = np.clip((outer_radius-radius)/width, 0., 1.)
    height = np.zeros(len(p))
    height[ids] = amplitude*smoothstep(depth)
    # Retain all literal outer wall and protected free-lining controls;
    # interpolation noise never authorizes moving one of these source walls.
    fixed = np.unique(np.r_[fixed,angular_fixed])
    height[fixed] = 0.
    all_direction = center-p
    all_direction /= np.maximum(np.linalg.norm(all_direction, axis=1)[:, None], 1e-20)
    moved = p+height[:, None]*all_direction
    origin = np.asarray(source['origins']['bladder'], float)
    encoded = (moved-origin).astype('<f4').astype(float)+origin
    moving = height > 0
    affected = {sid for sid, fs in surfaces.items() if np.any(moving[fs])}
    generated, normal_evidence = grouped_normal_words(encoded, {sid: surfaces[sid] for sid in affected},
        {sid: np.full(len(surfaces[sid]), 'fluid-lining' if sid in lining else sid, dtype=object) for sid in affected})
    source_parent_faces={}
    fields,static_normal_evidence=fixed_source_fields(p,encoded,surfaces,moving,source,parent_hints=parent_hints,parent_faces_out=source_parent_faces)
    for sid, values in generated.items():
        changed_faces = np.any(moving[surfaces[sid]], axis=1)
        fields[sid][changed_faces] = values[changed_faces]
    arrays, surface_keys, word_keys = dict(source_positions=p, positions=moved, encoded_positions=encoded, height=height, fixed_controls=fixed), {}, {}
    parent_keys={}
    for i, sid in enumerate(surfaces):
        surface_keys[sid] = 'surface'+str(i); arrays[surface_keys[sid]] = surfaces[sid]
        word_keys[sid] = 'words'+str(i); arrays[word_keys[sid]] = fields[sid]
        if sid in source_parent_faces:
            parent_keys[sid]='native_parent'+str(i);arrays[parent_keys[sid]]=source_parent_faces[sid]
    atomic_npz(output/'construction.npz', **arrays)
    plan = {k: v for k, v in previous.items() if k not in ['surfaces', 'corner_words', 'source_geometry', 'source_geometry_sha256', 'distance_query', 'normal_field', 'outward_failures']}
    plan.update(schema='bladder-lining-radial-relief/v3', source_geometry=str(output/'construction.npz'), source_geometry_sha256=digest_file(output/'construction.npz'),
                source_construction_sha256=digest_file(folder/'construction.npz'), surfaces=surface_keys, corner_words=word_keys,
                radial_source_queries=ray_evidence, protected_distance_query=distance_evidence, normal_field=normal_evidence,static_normal_source=static_normal_evidence,normal_source_parents=parent_keys,
                opening_cones=opening_cones,source_opening_cone_controls_fixed=len(angular_fixed),
                resolved_inner_band_transport=True,band_part_ids=sorted(band),coupled_refinement=refinement,source_band_guide_sha256=digest_file(source_carriers/'guide.json'),
                source_fluid_carrier=str(source_carriers/'fluid.npz'),source_fluid_carrier_sha256=digest_file(source_carriers/'fluid.npz'),
                source_outer_carrier=str(source_carriers/'outer.npz'),source_outer_carrier_sha256=digest_file(source_carriers/'outer.npz'),
                outer_wall_interface_literal=True, refinement_only_source=False,moved_controls=int(moving.sum()), actual_maximum_displacement_m=float(height.max()),
                minimum_resolved_wall_thickness_m=float(width.min()), source_ray_monotonicity_bound=float(np.min((1-amplitude/inner_radius)**2)),
                source_ray_bound_is_not_a_discrete_mesh_proof=True, installed=False, anatomical_acceptance=False)
    atomic_json(output/'construction.json', plan)
    if not validate:
        return dict(state='prepared',source_geometry=plan['source_geometry'],source_geometry_sha256=plan['source_geometry_sha256'],
                    moved_controls=plan['moved_controls'],installed=False)
    result = pack(moved, surfaces, cells, {pid: pid for pid in cells}, source['frame'], output/'packed', origin=origin, corner_words=fields, max_pairs=1000000)
    report = dict(**plan, packed=result, passed=result['passed'])
    atomic_json(output/'candidate.json', report)
    return report
