"""Coats across a joined bladder and ureters, with literal renal end charts.

The two UPJ annuli inherit their tissue roles from the installed reciprocal
interfaces. A common harmonic coordinate joins those controls to the actual
outer/inner sheets; construction tetrahedra are never renderer geometry.
"""
from collections import defaultdict
from pathlib import Path
import itertools
import json

import numpy as np

from .common import atomic_json, atomic_npz, digest_file
from .audit import topology
from .harmonic_wall import graph_values, tetrahedralize
from .implicit_cells import Arrangement
from .model import AtlasReader, Mesh
from .tissue import pack


ROLES = ('adventitia', 'muscularis', 'mucosa')


def repair_coordinate_orientation(coordinates_dir,output):
    """Reuse solved fields; repair only tetrahedral orientation from topology."""
    from .algorithms.tetrahedral_orientation import orient_complex
    coordinates_dir,output=map(lambda p:Path(p).resolve(),(coordinates_dir,output))
    if output.exists():
        raise FileExistsError(output)
    info=json.loads((coordinates_dir/'coordinates.json').read_text())
    for name in ('wall','body_wall','catalogue'):
        if digest_file(info[name])!=info[name+'_sha256']:
            raise ValueError('A source-bound coordinate input changed: '+name)
    with np.load(coordinates_dir/'coordinates.npz',allow_pickle=False)as z:
        arrays={key:z[key].copy()for key in z.files}
    try:
        arrays['tetrahedra'],proof=orient_complex(arrays['positions'],arrays['tetrahedra'],arrays['boundary'])
    except ValueError as error:
        if hasattr(error,'details'):
            output.mkdir(parents=True,exist_ok=True)
            atomic_json(output/'orientation-failure.json',dict(error.details,
                source_coordinates_sha256=digest_file(coordinates_dir/'coordinates.npz'),
                source_wall_sha256=info['wall_sha256']))
        raise
    output.mkdir(parents=True)
    atomic_npz(output/'coordinates.npz',**arrays)
    info.update(orientation_repair=proof,source_coordinates_sha256=digest_file(coordinates_dir/'coordinates.npz'))
    atomic_json(output/'coordinates.json',info)
    return proof


def outlet_membership(p, boundary, tets, non_outlet_points, outlet_points):
    """Transport explicit source ownership through only the new outlet collar.

    A height or side plane is not an anatomical dependency: it can assign a
    piece of the outlet to a ureter. The boundary controls retain their actual
    source identity; only newly constructed collar points interpolate.
    """
    non_outlet_points, outlet_points = set(non_outlet_points), set(outlet_points)
    if non_outlet_points & outlet_points:
        raise ValueError('Outlet and adjoining source ownership overlap')
    lookup = {tuple(q): i for i, q in enumerate(p)}
    controls = {lookup[q]: 0. for q in non_outlet_points if q in lookup}
    controls.update({lookup[q]: 1. for q in outlet_points if q in lookup})
    if set(controls.values()) != {0., 1.}:
        raise ValueError('Both explicit outlet ownership sources must survive')
    edges = np.concatenate([boundary[:, [i, j]] for i, j in [(0, 1), (1, 2), (2, 0)]])
    ids = np.array(list(controls), int)
    surface, surface_evidence = graph_values(p, edges, ids, np.array(list(controls.values())))
    volume_edges = np.concatenate([tets[:, [i, j]] for i, j in itertools.combinations(range(4), 2)])
    boundary_ids = np.unique(boundary)
    membership, evidence = graph_values(p, volume_edges, boundary_ids, surface[boundary_ids])
    evidence.update(boundary_graph=surface_evidence, fixed_source_controls=len(controls),
                    source_roles_explicit=True)
    if np.max(abs(membership[ids]-np.array(list(controls.values())))) > 1e-12:
        raise ValueError('The outlet ownership solve moved a literal source control')
    return membership, evidence


def refine_outlet_membership(coordinates_dir, outlet_wall, output):
    """Add source ownership to a saved tetrahedral construction without remeshing."""
    coordinates_dir, outlet_wall, output = map(lambda q: Path(q).resolve(),
                                              (coordinates_dir, outlet_wall, output))
    if output.exists():
        raise FileExistsError(output)
    info = json.loads((coordinates_dir/'coordinates.json').read_text())
    if not info.get('finite_urethral_outlet'):
        raise ValueError('This coordinate bank has no finite urethral outlet')
    for name in ('wall', 'body_wall', 'catalogue'):
        if digest_file(info[name]) != info[name+'_sha256']:
            raise ValueError('A bound urinary source changed: '+name)
    with np.load(coordinates_dir/'coordinates.npz', allow_pickle=False) as z:
        arrays = {k: z[k].copy() for k in z.files}
    body = Mesh.load(info['body_wall'])
    other = {tuple(q) for q in body.vertices}
    for side in ('left', 'right'):
        other.update(map(tuple, Mesh.load(Path(info['pose_dir'])/(side+'-transported-wall.npz')).vertices))
    outlet = Mesh.load(outlet_wall)
    membership, evidence = outlet_membership(arrays['positions'], arrays['boundary'], arrays['tetrahedra'],
                                             other, map(tuple, outlet.vertices))
    arrays['urethral_membership'] = membership
    output.mkdir(parents=True)
    atomic_npz(output/'coordinates.npz', **arrays)
    info.update(urethral_membership_harmonic=evidence, source_outlet_wall=str(outlet_wall),
                source_outlet_wall_sha256=digest_file(outlet_wall),
                reused_tetrahedra_sha256=digest_file(coordinates_dir/'coordinates.npz'))
    atomic_json(output/'coordinates.json', info)
    return dict(source_controls=evidence['fixed_source_controls'],
                reused_tetrahedra=len(arrays['tetrahedra']), geometry_changed=False)


def root_contract(atlas, wall):
    """Read actual installed end facets, not a planar or radial approximation."""
    reader = AtlasReader()
    source_faces = {tuple(sorted(tuple(q) for q in wall.vertices[t])): i
                    for i, t in enumerate(wall.faces)}
    roles, vertices = {}, defaultdict(set)
    for side in ('left', 'right'):
        for role in ROLES:
            pid = f'RENALDETAIL:male:{side}:ureter-wall:{role}'
            record = next(p for p in atlas['parts'] if p['id'] == pid)
            mesh = reader.read_part(atlas, record, frame=wall.frame)
            groups = [g for g in record['tissueDrawGroups']
                      if g.get('neighbor') == f'RENALDETAIL:male:{side}:upj:{role}']
            if not groups:
                raise ValueError('Missing reciprocal renal-root coat interface: ' + pid)
            for group in groups:
                fs = mesh.faces[group['start']//3:(group['start']+group['count'])//3]
                for tri in mesh.vertices[fs]:
                    key = tuple(sorted(tuple(q) for q in tri))
                    if key not in source_faces:
                        raise ValueError('Joined wall changed a literal renal-root facet')
                    index = source_faces[key]
                    if index in roles:
                        raise ValueError('Renal-root tissue faces overlap')
                    roles[index] = (side, role)
                    for q in tri:
                        vertices[tuple(q)].add(role)
    values = {}
    for point, names in vertices.items():
        if names == {'adventitia', 'muscularis'}:
            value = 1/3
        elif names == {'muscularis', 'mucosa'}:
            value = 2/3
        elif len(names) == 1:
            value = (ROLES.index(next(iter(names))) + .5) / 3
        else:
            raise ValueError('Collapsed or conflicting renal-root coat chart')
        values[point] = value
    return roles, values


def prepare(joined_dir, body_wall, pose_dir, catalogue, output, *, body_carriers=None, outlet_wall=None):
    joined_dir, pose_dir = Path(joined_dir).resolve(), Path(pose_dir).resolve()
    body_wall, catalogue, output = Path(body_wall).resolve(), Path(catalogue).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    report_path=joined_dir/('junctions.json' if (joined_dir/'junctions.json').exists() else 'outlet-with-carriers.json')
    report = json.loads(report_path.read_text())
    if report.get('passed') is not True:
        raise ValueError('Use a passed continuous urinary wall')
    original = Mesh.load(joined_dir/'urinary-wall.npz')
    if original.geometry_sha!=report['whole_check']['geometry_sha256']:
        raise ValueError('The passed whole-wall report does not bind the current geometry')
    wall_path = joined_dir/'verified-carriers/wall.npz'
    wall = Mesh.load(wall_path)
    if wall.surface_sha != original.surface_sha:
        raise ValueError('Carrier compaction changed the continuous wall surface')
    with np.load(joined_dir/'verified-carriers/carriers.npz', allow_pickle=False) as z:
        if not np.array_equal(wall.faces, z['wall_faces']):
            raise ValueError('Joined wall/carrier face sequence changed')
        inner = z['inner_side_faces'].copy()
        inner_ids = z['inner_source_face_ids'].copy()
        rest = z['outer_and_rim_faces'].copy()
        neck = wall.vertices[z['urethral_neck_ring']].mean(0) if 'urethral_neck_ring' in z else None
    if body_carriers:
        body_source=Mesh.load(body_wall)
        with np.load(body_carriers,allow_pickle=False)as z:
            if not np.array_equal(body_source.faces,z['wall_faces']):
                raise ValueError('The source bladder landmark connectivity changed')
            neck=body_source.vertices[z['urethral_neck_ring']].mean(0)
    if neck is None:
        raise ValueError('A source-bound bladder-neck landmark is required')
    atlas = json.loads(catalogue.read_text())
    root_faces, root_values = root_contract(atlas, wall)
    p, tets, _, _, mesher = tetrahedralize(wall.vertices, wall.faces)
    lookup = {tuple(q): i for i, q in enumerate(p)}
    boundary = np.array([[lookup[tuple(wall.vertices[i])] for i in t] for t in wall.faces], int)
    lining = {tuple(q) for q in wall.vertices[np.unique(inner)]}
    # The physical end annuli are free except their explicit tissue chart.
    # Inner lips stay phi=1 even though their vertices also touch an annulus.
    annuli = set(root_faces)
    outside = {tuple(q) for q in wall.vertices[np.unique(rest)]} - lining
    root_points = set(root_values)
    phi_constraints = {q: 0. for q in outside-root_points}
    phi_constraints.update(root_values)
    phi_constraints.update({q: 1. for q in lining})
    # Actual outer rim controls, unlike intermediate annular controls, lie on
    # a non-annular outer source facet and must remain phi=0.
    outer_sheet = np.setdiff1d(np.arange(len(wall.faces)), np.r_[inner_ids, list(annuli)])
    exterior = {tuple(q) for q in wall.vertices[np.unique(wall.faces[outer_sheet])]} - lining
    phi_constraints.update({q: 0. for q in exterior})
    edges = np.concatenate([boundary[:, [i, j]] for i, j in [(0, 1), (1, 2), (2, 0)]])
    fixed = np.array([lookup[q] for q in phi_constraints], int)
    surface, surface_evidence = graph_values(p, edges, fixed, np.array(list(phi_constraints.values())))
    volume_edges = np.concatenate([tets[:, [i, j]] for i, j in itertools.combinations(range(4), 2)])
    boundary_ids = np.unique(boundary)
    phi, harmonic = graph_values(p, volume_edges, boundary_ids, surface[boundary_ids])

    # Bladder/ureter transition follows source lineage. Only new collar
    # controls interpolate; source body and ureter controls keep their roles.
    body = Mesh.load(body_wall)
    body_points = {tuple(q) for q in body.vertices}
    tube_points = set()
    for side in ('left', 'right'):
        tube_points.update(tuple(q) for q in Mesh.load(pose_dir/(side+'-transported-wall.npz')).vertices)
    outlet_assignments={};outlet_limits=None
    if outlet_wall:
        outlet_wall=Path(outlet_wall).resolve();outlet=Mesh.load(outlet_wall)
        with np.load(outlet_wall.parent/'source-lineage.npz',allow_pickle=False)as z:
            stations=z['source_station_ids'].copy();arc=z['source_station_arclength_m'].copy();axis=z['native_course_axis'].copy()
        if len(outlet.vertices)!=2*len(stations):
            raise ValueError('The finite urethral wall lost its native course correspondence')
        values=np.tile(arc[stations]/arc[-1],2)
        outlet_assignments={tuple(q):float(t)for q,t in zip(outlet.vertices,values)}
        tube_points.update(outlet_assignments)
        outlet_limits=[]
        for level in (-.008667132817208767,.0070313639007508755):
            crossings=np.flatnonzero((axis[:-1,2]-level)*(axis[1:,2]-level)<=0)
            if not len(crossings):
                raise ValueError('A retained native urethral region boundary is absent')
            i=int(crossings[0]);t=(level-axis[i,2])/(axis[i+1,2]-axis[i,2])
            outlet_limits.append(float(((1-t)*arc[i]+t*arc[i+1])/arc[-1]))
    if body_points & tube_points:
        raise ValueError('Bladder and ureter source controls are ambiguously shared')
    assignments = {lookup[q]: 1. for q in body_points if q in lookup}
    assignments.update({lookup[q]: 0. for q in tube_points if q in lookup})
    if not assignments or 0. not in assignments.values() or 1. not in assignments.values():
        raise ValueError('Both tissue sources must constrain the collar transition')
    ids = np.array(list(assignments), int)
    surface_b, transition_surface = graph_values(p, edges, ids, np.array(list(assignments.values())))
    blend, transition = graph_values(p, volume_edges, boundary_ids, surface_b[boundary_ids])
    course=np.zeros(len(p));course_evidence=None
    membership=np.zeros(len(p));membership_evidence=None
    if outlet_assignments:
        controls={lookup[q]:0. for q in body_points|tube_points if q in lookup}
        controls.update({lookup[q]:t for q,t in outlet_assignments.items() if q in lookup})
        ids=np.array(list(controls),int)
        surface_q,sq=graph_values(p,edges,ids,np.array(list(controls.values())))
        course,course_evidence=graph_values(p,volume_edges,boundary_ids,surface_q[boundary_ids])
        course_evidence['boundary_graph']=sq
        membership,membership_evidence=outlet_membership(p,boundary,tets,
            (body_points|tube_points)-set(outlet_assignments),outlet_assignments)
    output.mkdir(parents=True)
    atomic_npz(output/'coordinates.npz', positions=p, tetrahedra=tets, boundary=boundary,
        phi=phi, bladder_blend=blend, source_inner_face_ids=inner_ids,
        root_face_ids=np.array(sorted(root_faces), int), neck=neck,urethral_course_fraction=course,
        urethral_membership=membership)
    result = dict(schema='continuous-urinary-wall-coordinates/v1',
        wall=str(wall_path), wall_sha256=digest_file(wall_path),
        continuous_wall_sha256=digest_file(joined_dir/'urinary-wall.npz'),
        body_wall=str(body_wall), body_wall_sha256=digest_file(body_wall),
        catalogue=str(catalogue), catalogue_sha256=digest_file(catalogue),
        pose_dir=str(pose_dir), source_root_faces={str(k): list(v) for k, v in root_faces.items()},
        finite_urethral_outlet=bool(outlet_assignments),native_urethral_region_limits=outlet_limits,
        urethral_course_harmonic=course_evidence,
        urethral_membership_harmonic=membership_evidence,
        source_outlet_wall=str(outlet_wall)if outlet_wall else None,
        source_outlet_wall_sha256=digest_file(outlet_wall)if outlet_wall else None,
        mesher=mesher, harmonic=harmonic, surface_harmonic=surface_evidence,
        transition_harmonic=transition, transition_surface_harmonic=transition_surface,
        source_tissue_roles_explicit=True, source_renal_end_facets_literal=True,
        clinical_thickness_measured=False, installed=False)
    atomic_json(output/'coordinates.json', result)
    return result


def construct(coordinates_dir, source_carriers, output, *, arrangement_only=False, feature_plan=None, opening_surround_width_m=None):
    coordinates_dir, output = Path(coordinates_dir).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    info = json.loads((coordinates_dir/'coordinates.json').read_text())
    if digest_file(info['wall']) != info['wall_sha256']:
        raise ValueError('The continuous source wall changed')
    wall = Mesh.load(info['wall'])
    with np.load(coordinates_dir/'coordinates.npz', allow_pickle=False) as z:
        p, tets, boundary = z['positions'].copy(), z['tetrahedra'].copy(), z['boundary'].copy()
        phi, b, inner_ids, neck = z['phi'].copy(), z['bladder_blend'].copy(), z['source_inner_face_ids'].copy(), z['neck'].copy()
        course=z['urethral_course_fraction'].copy() if 'urethral_course_fraction'in z else np.zeros(len(p))
        membership=z['urethral_membership'].copy() if 'urethral_membership'in z else None
    if info.get('finite_urethral_outlet') and membership is None:
        raise ValueError('Refine the saved coordinates with explicit outlet ownership before partitioning')
    body = Mesh.load(info['body_wall'])
    if digest_file(info['body_wall']) != info['body_wall_sha256']:
        raise ValueError('The source bladder changed')
    with np.load(source_carriers, allow_pickle=False) as z:
        left, right = [body.vertices[z[side+'_ureteric_ring']].mean(0) for side in ('left', 'right')]
    basis = np.column_stack([left-neck, right-neck])
    uv = (p-neck) @ np.linalg.pinv(basis).T
    bary = np.c_[1.-uv.sum(1), uv]
    normal = np.cross(left-neck, right-neck); normal /= np.linalg.norm(normal)
    if normal[2] < 0:
        normal *= -1
    distance = (p-neck) @ normal
    def regional(value, outside=1.):
        # A constant sign outside the reservoir avoids needlessly splitting
        # the immutable renal annuli with distant bladder regional planes.
        return b*value + (1-b)*outside
    coat_outer = (1-b)/3 + b*.08
    coat_muscle = (1-b)*(2/3) + b*.78
    fields = np.column_stack([
        np.full(len(p), -1.), phi, b-.5,
        phi-coat_outer, phi-coat_muscle, regional(phi-.94),
        regional(phi-.015),
        *[regional(phi-c) for c in (.30, .55, .45, .60)],
        regional((p[:, 1]-neck[1]-.008)/.01),
        regional((p[:, 1]-(left[1]+right[1])/2-.004)/.01),
        regional((.010-distance)/.01),
        regional((max(left[1], right[1])+.004-p[:, 1])/.01),
        *[regional(bary[:, i]+.08) for i in range(3)],
        regional((.008-distance)/.01), regional((distance+.012)/.01),
        p[:, 0]+.2*course,
    ])
    has_outlet=info.get('finite_urethral_outlet',False)
    if has_outlet:
        limits=info['native_urethral_region_limits']
        # Compare outlet membership with the other conduit membership. The
        # reservoir itself is still selected by b, independently of height.
        fields=np.c_[fields,course-limits[0],course-limits[1],2*membership+b-1]
    feature_indices=None
    if feature_plan:
        from .bladder_features import TrigonalRelief
        feature_plan=Path(feature_plan).resolve();features=json.loads(feature_plan.read_text())
        feature_source=Mesh.load(feature_plan.parent/'urinary-wall.npz')
        if features['catalogue_sha256']!=info['catalogue_sha256']or feature_source.surface_sha!=wall.surface_sha:
            raise ValueError('The feature plan does not bind this whole bladder wall and source catalogue')
        with np.load(source_carriers,allow_pickle=False)as z:
            inner_source=z['inner_side_faces'].copy();outer_source=z['outer_and_rim_faces'].copy()
        settings=features['relief']
        relief=TrigonalRelief(body,inner_source,outer_source,
            {name:np.array(q,float)for name,q in features['actual_ostia'].items()},
            crest_height_m=settings['crest_height_m'],crest_half_width_m=settings['crest_half_width_m'],
            uvula_height_m=settings['uvula_height_m'],falloff_into_fluid_m=settings['fluid_falloff_m'])
        amplitude=relief.amplitudes(p[:,:2])
        feature_indices=(fields.shape[1],fields.shape[1]+1)
        # These are finite portions of the existing mucosal tissue under the
        # actual elevated surfaces, not overlapping ridge solids. The front
        # wall at the same X/Y is excluded by an explicit posterior bound.
        feature_fields=np.minimum(amplitude/np.array([settings['crest_height_m'],settings['uvula_height_m']])-.15,
                                  ((-.010-p[:,2])/.01)[:,None])
        fields=np.c_[fields,feature_fields]
    opening_indices={};opening_planes={};opening_windows={}
    if opening_surround_width_m is not None:
        if not feature_plan or not .0005 <= opening_surround_width_m <= .002:
            raise ValueError('Source-defined ostia and an explicit 0.5–2 mm mucosal surround width are required')
        grid=2.**-20
        fixed_root=wall.triangles[list(map(int,info['source_root_faces']))]
        for name in ('left','right','neck'):
            rim=np.asarray(features['actual_ostia'][name],float)
            low=np.floor((rim.min(0)-opening_surround_width_m)/grid)*grid
            high=np.ceil((rim.max(0)+opening_surround_width_m)/grid)*grid
            indices=[]
            for axis in range(3):
                for threshold,sign in ((low[axis],1),(high[axis],-1)):
                    q=fixed_root[:,:,axis]-threshold
                    if np.any((q.min(1)<0)&(q.max(1)>0)):
                        raise ValueError('An opening selection plane would cut an immutable renal-root facet')
                    index=fields.shape[1];indices.append(index);opening_planes[index]=(axis,float(threshold),sign)
                    fields=np.c_[fields,sign*(p[:,axis]-threshold)]
            opening_indices[name]=indices;opening_windows[name]=dict(low_m=low.tolist(),high_m=high.tolist())
    def label(x):
        if x[2] < 0:
            role = 'adventitia' if x[3] < 0 else 'muscularis' if x[4] < 0 else 'mucosa'
            if has_outlet and x[23]>0:
                region='prostatic' if x[21]<0 else 'intermediate' if x[22]<0 else 'spongy'
                return 'urethra-'+region+':'+role
            side = 'left' if x[20] >= 0 else 'right'
            return 'ureter-'+side+':'+role
        is_neck = x[11] <= 0
        trigone = bool(np.all(x[15:20] >= 0))
        region = 'neck' if is_neck else 'fundus' if trigone or x[13] >= 0 and x[14] >= 0 else 'body'
        if x[3] < 0:
            tissue = ('serosa' if x[6] < 0 else 'subserosa') if x[12] >= 0 else 'adventitia'
        elif x[4] < 0:
            if is_neck:
                tissue = 'neck-external-longitudinal' if x[7] < 0 else 'neck-circular' if x[8] < 0 else 'neck-internal-longitudinal'
            elif trigone and x[9] >= 0:
                tissue = 'trigone-deep' if x[10] < 0 else 'trigone-superficial'
            else:
                tissue = 'detrusor-unstratified'
        elif x[5] < 0:
            tissue = 'submucosa'
        else:
            if feature_indices:
                if x[feature_indices[0]]>0:
                    return 'bladder-'+region+':mucosa-interureteric-crest'
                if x[feature_indices[1]]>0:
                    return 'bladder-'+region+':mucosa-uvula'
            for name,indices in opening_indices.items():
                if np.all(x[indices]>0):
                    role='mucosa-internal-urethral-orifice' if name=='neck' else 'mucosa-ureteric-orifice-'+name
                    return 'bladder-'+region+':'+role
            tissue = 'mucosa-trigone' if trigone else 'mucosa'
        return 'bladder-'+region+':'+tissue
    origin = np.round(p.mean(0)*1024)/1024
    for index in info.get('orientation_repair',{}).get('flat_interface_glue_tetrahedra',[]):
        values=fields[tets[index]]
        active=[k for k in np.flatnonzero((values.min(0)<0)&(values.max(0)>0))if k not in (0,1)]
        if active or len({label(v)for v in values})!=1:
            raise ValueError('A flat construction transition crosses tissue ownership; repair the volume complex first')
    local_planes={k:(axis,threshold-origin[axis],sign)for k,(axis,threshold,sign)in opening_planes.items()}
    arrangement = Arrangement(p-origin, fields, label,axis_planes=local_planes)
    arrangement.facet_labels = {tuple(sorted(face)): i for i, face in enumerate(boundary)}
    try:
        positions, surfaces, cells = arrangement.build(tets, skip_fields=[0, 1],
            progress=lambda i,n,v,l: print('partition',i,n,'controls',v,flush=True))
    except ValueError as error:
        if hasattr(error,'details'):
            output.mkdir(parents=True,exist_ok=True)
            atomic_json(output/'partition-failure.json',dict(error.details,
                source_coordinates_sha256=digest_file(coordinates_dir/'coordinates.npz')))
            atomic_npz(output/'partition-failure.npz',positions=np.array(arrangement.points),values=np.array(arrangement.values),
                       source_positions=p,source_tetrahedra=tets)
        raise
    positions += origin
    positions[:len(p)] = p
    output.mkdir(parents=True)
    arrays = {'positions': positions}
    surface_keys = {}; tag_keys = {}
    for i, (sid, fs) in enumerate(surfaces.items()):
        key = 'surface_'+str(i); arrays[key] = fs; surface_keys[sid] = key
        key = 'source_faces_'+str(i); arrays[key] = arrangement.surface_source_faces[sid]; tag_keys[sid] = key
    atomic_npz(output/'arrangement.npz', **arrays)
    roots = {int(k): tuple(v) for k, v in info['source_root_faces'].items()}
    root_actual = defaultdict(list)
    lining_ids = set(map(int, inner_ids)); lined = set()
    for name, uses in cells.items():
        for use in uses:
            sid = use['surface']; tags = arrangement.surface_source_faces[sid]
            for face, tag in zip(surfaces[sid], tags):
                tag = int(tag)
                if tag in roots:
                    side, role = roots[tag]
                    if name != 'ureter-'+side+':'+role or use['neighbor'] is not None:
                        raise ValueError('A renal-root facet changed tissue ownership')
                    root_actual[tag].append(tuple(sorted(tuple(q) for q in positions[face])))
                if tag in lining_ids:
                    if not name.rsplit(':',1)[-1].startswith('mucosa'):
                        raise ValueError('A non-mucosal tissue owns the urinary lining')
                    lined.add(tag)
    if lined != lining_ids:
        raise ValueError('The actual urinary lining is incompletely represented')
    for index in roots:
        expected = tuple(sorted(tuple(q) for q in wall.vertices[wall.faces[index]]))
        if root_actual[index] != [expected]:
            raise ValueError('An existing renal-root facet was subdivided or changed')
    conduit_topology={}
    for name,uses in cells.items():
        if not name.startswith('ureter-'):
            continue
        fs=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
        check=topology(Mesh(positions,fs,wall.frame));conduit_topology[name]=check
        if not check['closed_oriented_manifold'] or len(check['components'])!=1:
            raise ValueError('A ureter coat must be one closed connected tissue: '+name)
    part_ids = {name: ('RENALDETAIL:male:'+name.split(':')[0].split('-')[1]+':ureter-wall:'+name.split(':')[1]
                     if name.startswith('ureter-') else 'URETHRA:male:'+name[8:] if name.startswith('urethra-')
                     else 'BLADDER:male:'+name[8:]) for name in cells}
    source = dict(schema='continuous-urinary-wall-arrangement/v1', frame=wall.frame,
        surfaces=surface_keys, source_face_arrays=tag_keys, cells=cells, part_ids=part_ids,
        root_faces_preserved=len(roots), all_inner_faces_mucosal=True,
        conduit_topology=conduit_topology,explicit_outlet_membership=has_outlet,
        source_wall_sha256=info['wall_sha256'], coordinates_sha256=digest_file(coordinates_dir/'coordinates.npz'),
        feature_plan=str(feature_plan)if feature_plan else None,
        feature_plan_sha256=digest_file(feature_plan)if feature_plan else None,
        opening_surround_width_m=opening_surround_width_m,
        opening_support_method='axis-exact-rim-window'if opening_indices else None,
        opening_windows=opening_windows,
        opening_curves=features['actual_ostia'] if opening_indices else None,
        anatomical_boundaries_provisional=True, installed=False, anatomical_acceptance=False)
    atomic_json(output/'arrangement.json', source)
    if arrangement_only:
        return dict(regions=list(cells), source_root_faces_preserved=len(roots))
    result = pack(positions, surfaces, cells, part_ids, wall.frame, output/'tissues', origin=origin,
                  max_pairs=1000000, construction=dict(coordinates=info,
                    source_renal_end_facets_literal=True, complete_mucosal_lining=True,
                    source_cell_ownership_provisional_at_collars=True))
    report = dict(schema='continuous-urinary-wall-coats/v1', passed=result['passed'],
        regions=list(cells), part_ids=part_ids, geometry_bytes=result['total_bytes'],
        root_faces_preserved=len(roots), all_inner_faces_mucosal=True,
        source_wall_sha256=info['wall_sha256'], coordinates_sha256=digest_file(coordinates_dir/'coordinates.npz'),
        installed=False, anatomical_acceptance=False)
    atomic_json(output/'layers.json', report)
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    for flag in ('joined-dir', 'body-wall', 'pose-dir', 'catalogue', 'output'):
        p.add_argument('--'+flag, required=True)
    p.add_argument('--body-carriers');p.add_argument('--outlet-wall')
    p = sub.add_parser('refine-outlet')
    p.add_argument('coordinates_dir');p.add_argument('--outlet-wall',required=True);p.add_argument('--output',required=True)
    p = sub.add_parser('orient-coordinates')
    p.add_argument('coordinates_dir');p.add_argument('--output',required=True)
    p = sub.add_parser('construct')
    p.add_argument('coordinates_dir'); p.add_argument('--source-carriers', required=True); p.add_argument('--output', required=True)
    p.add_argument('--arrangement-only', action='store_true')
    p.add_argument('--feature-plan')
    p.add_argument('--opening-surround-width-m',type=float)
    args = vars(parser.parse_args()); command = args.pop('command')
    result = prepare(**args) if command == 'prepare' else refine_outlet_membership(**args)if command=='refine-outlet'else repair_coordinate_orientation(**args)if command=='orient-coordinates'else construct(**args)
    print(json.dumps({k:v for k,v in result.items() if k in ('passed', 'regions', 'geometry_bytes', 'mesher','orientations_changed','flat_interface_glue_tetrahedra')}, indent=2))
