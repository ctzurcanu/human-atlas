"""Three bounded, read-only placement stages; no donor surface export."""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import igl
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.reference_assets import ReferenceAssets
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts

SOURCE_SHA = '5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'
BASE = Path('/Users/christiantzurcanu/Documents/dev/3d/ext/meshes')
OUT = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-facial-family-inputs-v1'
CIDS = [20537, 21433, 23808, 22122, 20872, 17991, 21461, 22493, 22534, 24587, 24815, 18490, 18827, 16872, 18832, 16525, 16922, 18755, 16667, 19565]


def sources():
    reader = ReferenceAssets(BASE/'629250ec501a0ed64e00bf2e632c09a4/assets.json', ROOT/'.local-models/models/Male_2024.glb', BASE/'en.json')
    rows = reader.resolve(CIDS)
    arrays = {r['node_id']: reader.node_mesh(r) for r in rows}
    source = ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(source) == SOURCE_SHA
    atlas = json.loads(source.read_bytes())
    native_reader = AtlasReader()
    native = {}
    ids = {'ZA:Mandible:Bone', 'ZA:Occipital bone:Bone-6', 'ZA:Temporal bone.l:Bone-1', 'ZA:Temporal bone.r:Bone-1',
           'ZA:Hyoid bone', 'ZA:Superior pharyngeal constrictor.l', 'ZA:Stylopharyngeus muscle.l',
           'ZA:Medial pterygoid muscle.l', 'TA98TONGUE:male:styloglossus:left', 'ZA:Facial artery.l',
           'ZA:Facial vein.l', 'ZA:Anterior division of retromandibular vein.l', 'ZA:Submandibular gland.l', 'ZA:Palatine tonsil.l'}
    for part in atlas['parts']:
        if part['id'] in ids or (part['id'].startswith('TA98HYOID:') and part.get('system') == 'skeletal' and not part.get('suppressed')):
            native[part['id']] = native_reader.read_part(atlas, part, frame='male-native-za-world-m', catalogue_path=source, catalogue_sha=SOURCE_SHA)
    return reader, rows, arrays, native


def stage1():
    began = time.monotonic()
    reader, rows, arrays, native = sources()
    file = OUT/'reference-node-placement-metadata-v6.json'
    assert not file.exists()
    for row in rows:
        v, f = arrays[row['node_id']]
        used = v[np.unique(f)]
        row.update(vertices=len(v), triangles=len(f), bounds=np.array([used.min(0), used.max(0)]).tolist(),
                   used_vertex_centroid=used.mean(0).tolist(), linear_transform_determinant=float(np.linalg.det(np.asarray(row['node_transform'])[:3,:3])))
    result = dict(schema='supplied-reference-node-placement/v1', source_catalogue_sha256=SOURCE_SHA,
                  source=reader.provenance(), reference_nodes=rows,
                  native_landmark_sources=[dict(id=k, geometry_sha256=m.geometry_sha, bounds=m.bounds.tolist()) for k,m in native.items()],
                  side_evidence=[dict(cid=r['cid'], node_id=r['node_id'], attr=r['attr'], x_interval=[r['bounds'][0][0],r['bounds'][1][0]]) for r in rows if r['cid'] in [20537,21433,23808,21461]],
                  axes_and_units='Not assumed from glTF. Bilateral metadata + homologous native landmarks will resolve side, superior/anterior and metric scale in stage 2.',
                  redistribution_license_verified=False, private_reference_only=True, no_geometry_or_candidate_changed=True,
                  elapsed_seconds=time.monotonic()-began)
    atomic_json(file, result)
    print(json.dumps(dict(file=str(file.relative_to(ROOT)), nodes=len(rows), triangles_decoded=sum(r['triangles'] for r in rows),
                         branches=[{k:r[k] for k in ['cid','node_id','attr','bounds','triangles']} for r in rows if r['cid'] in [20537,21433,23808]],
                         elapsed_seconds=result['elapsed_seconds'])))


def landmark_points(vertices, family):
    """Actual source-vertex witnesses in explicitly homologous bone regions."""
    v = np.asarray(vertices)
    v = np.unique(v, axis=0)
    low, high = v.min(0), v.max(0)
    middle = (low+high)/2
    width = max(abs(low[0]), abs(high[0]))
    rows = []
    def pick(name, mask, axis, maximum):
        ids = np.flatnonzero(mask)
        assert len(ids), (family, name)
        i = ids[np.argmax(v[ids, axis]) if maximum else np.argmin(v[ids, axis])]
        rows.append(dict(name=name, family=family, actual_point=v[i].tolist(), unique_vertex_index=int(i),
                         eligible_source_vertices=len(ids), witness_definition=dict(axis=axis, maximum=maximum)))
    central = abs(v[:,0]) < .08*width
    if family == 'mandible':
        pick('inferior_midline_mandible', central & (v[:,2] > middle[2]), 1, False)
        pick('anterior_low_symphysis', central & (v[:,1] < low[1]+.4*(high[1]-low[1])), 2, True)
        for suffix, sign in [('left',1),('right',-1)]:
            side = sign*v[:,0] > .65*width
            posterior = v[:,2] < np.quantile(v[side,2], .35)
            pick('superior_posterior_condyle_'+suffix, side & posterior, 1, True)
            pick('inferior_posterior_angle_'+suffix, side & posterior, 1, False)
    elif family == 'hyoid':
        pick('anterior_hyoid_body', central, 2, True)
        for suffix, sign in [('left',1),('right',-1)]:
            pick('posterior_greater_horn_'+suffix, sign*v[:,0]>.55*width, 2, False)
    elif family == 'occipital':
        pick('posterior_midline_occiput', central, 2, False)
        pick('superior_midline_occiput', central, 1, True)
    elif family.startswith('temporal-'):
        pick('inferior_styloid_tip_'+family.split('-')[1], np.ones(len(v),bool), 1, False)
    return rows


def proper_similarity(source, target):
    source, target = np.asarray(source), np.asarray(target)
    a, b = source.mean(0), target.mean(0)
    x, y = source-a, target-b
    u, singular, vt = np.linalg.svd(x.T@y)
    diagonal = np.array([1.,1.,np.linalg.det(u@vt)])
    rotation = u@np.diag(diagonal)@vt
    scale = float((singular*diagonal).sum()/(x*x).sum())
    translation = b-scale*a@rotation
    assert scale > 0 and np.linalg.det(rotation) > .999999
    return scale, rotation, translation


def stage2():
    began = time.monotonic()
    reader, rows, arrays, native = sources()
    file = OUT/'reference-native-similarity-v6.json'
    assert not file.exists()
    previous = json.loads((OUT/'reference-node-placement-metadata-v6.json').read_bytes())
    assert previous['source']['assets_sha256'] == reader.assets_sha256 and previous['source']['glb_json_sha256'] == reader.json_sha256
    def reference(cid, attr=None):
        row = next(r for r in rows if r['cid']==cid and r['attr']==attr)
        v,f = arrays[row['node_id']]
        return v[np.unique(f)], row
    ref_specs = [('mandible',22122,None),('hyoid',20872,None),('occipital',17991,None),('temporal-left',21461,4),('temporal-right',21461,8)]
    native_specs = {'mandible':['ZA:Mandible:Bone'], 'hyoid':[k for k in native if k.startswith('TA98HYOID:')],
                    'occipital':['ZA:Occipital bone:Bone-6'], 'temporal-left':['ZA:Temporal bone.l:Bone-1'], 'temporal-right':['ZA:Temporal bone.r:Bone-1']}
    pairs = []
    for family,cid,attr in ref_specs:
        v, reference_row = reference(cid,attr)
        first = landmark_points(v,family)
        target_vertices = np.concatenate([native[k].vertices[np.unique(native[k].faces)] for k in native_specs[family]])
        second = {r['name']:r for r in landmark_points(target_vertices,family)}
        for item in first:
            target = second[item['name']]
            pairs.append(dict(name=item['name'], family=family, reference=item['actual_point'], native_m=target['actual_point'],
                              source_node_id=reference_row['node_id'], source_cid=cid, native_ids=native_specs[family],
                              reference_witness=item, native_witness=target,
                              qualification='Geometric anatomical-region vertex witness, not a manually reviewed clinical landmark'))
    x, y = np.array([r['reference'] for r in pairs]), np.array([r['native_m'] for r in pairs])
    scale, rotation, translation = proper_similarity(x,y)
    predicted = scale*x@rotation+translation
    errors = np.linalg.norm(predicted-y,axis=1)
    for row,point,error in zip(pairs,predicted,errors):
        row.update(transformed_reference_m=point.tolist(), residual_mm=float(error*1000))
    holdout = []
    for family,_,_ in ref_specs:
        mask=np.array([r['family']!=family for r in pairs])
        s,r,t=proper_similarity(x[mask],y[mask])
        error=np.linalg.norm(s*x[~mask]@r+t-y[~mask],axis=1)*1000
        holdout.append(dict(held_out_family=family, fitted_point_count=int(mask.sum()), held_out_point_count=int((~mask).sum()),
                            rms_mm=float(np.sqrt(np.mean(error**2))), maximum_mm=float(error.max())))
    result=dict(schema='proper-reference-native-similarity/v1',source_catalogue_sha256=SOURCE_SHA,source=reader.provenance(),
                transform_convention='native row point in metres = scale * reference row point @ proper_rotation + translation_m',
                scale_metres_per_reference_unit=scale, proper_rotation=rotation.tolist(), translation_m=translation.tolist(),
                rotation_determinant=float(np.linalg.det(rotation)), native_axes_per_reference_positive_basis=rotation.tolist(),
                numerical_handedness_preserved=True, numerical_axis_qualification='Reference attr4/positiveX is left; skull superior positions establish positiveY superior; anterior symphysis/hyoid establish positiveZ anterior. Proper transform is constrained to determinant+1.',
                source_unit_declared=None, source_unit_inference='Metric scale inferred from actual homologous bone landmarks; no undocumented source-inch/metre assumption.',
                rms_residual_mm=float(np.sqrt(np.mean(errors**2))*1000), maximum_residual_mm=float(errors.max()*1000),
                landmark_pairs=pairs, leave_one_bone_family_out=holdout,
                registration_accepted_for_precise_native_routes=False,
                limitation='Cross-model anatomy and geometric extremum witnesses differ; this transform supplies coarse placement guidance only. Native muscular corridors must be rebuilt and checked independently.',
                no_geometry_or_candidate_changed=True, elapsed_seconds=time.monotonic()-began)
    atomic_json(file,result)
    print(json.dumps({k:result[k] for k in ['scale_metres_per_reference_unit','proper_rotation','translation_m','rotation_determinant','rms_residual_mm','maximum_residual_mm','leave_one_bone_family_out','elapsed_seconds']}))


def stage3():
    began=time.monotonic()
    reader,rows,arrays,native=sources()
    file=OUT/'reference-coarse-facial-corridors-v6.json'
    assert not file.exists()
    registration=json.loads((OUT/'reference-native-similarity-v6.json').read_bytes())
    assert registration['source']['assets_sha256']==reader.assets_sha256 and registration['source_catalogue_sha256']==SOURCE_SHA
    scale=registration['scale_metres_per_reference_unit']
    rotation=np.asarray(registration['proper_rotation']);translation=np.asarray(registration['translation_m'])
    def get(cid,attr=4):
        row=next(r for r in rows if r['cid']==cid and r['attr']==attr)
        v,f=arrays[row['node_id']]
        used=np.unique(f)
        unique,inverse=np.unique(v[used],axis=0,return_inverse=True)
        lookup=np.full(len(v),-1,int);lookup[used]=inverse
        return row,Mesh(unique,lookup[f],'supplied-male2024-reference-local-units')
    parent_row,parent=get(18755)
    contexts=[('superior-constrictor',22493,'ZA:Superior pharyngeal constrictor.l'),
              ('styloglossus',24587,'TA98TONGUE:male:styloglossus:left'),
              ('stylopharyngeus',24815,'ZA:Stylopharyngeus muscle.l'),
              ('medial-pterygoid-superficial',18490,'ZA:Medial pterygoid muscle.l'),
              ('medial-pterygoid-deep',18827,'ZA:Medial pterygoid muscle.l'),
              ('facial-vein',16872,'ZA:Facial vein.l'),
              ('retromandibular-anterior-vein',18832,'ZA:Anterior division of retromandibular vein.l')]
    targets=[]
    for name,cid,native_id in contexts:
        row,mesh=get(cid)
        targets.append((name,row,mesh,native[native_id],native_id))
    branches=[]
    for cid,name in [(20537,'ascending-palatine'),(21433,'tonsillar'),(23808,'glandular')]:
        row,mesh=get(cid)
        d2,parent_face,parent_points=igl.point_mesh_squared_distance(mesh.vertices,parent.vertices,parent.faces)
        root_index=int(np.argmin(d2))
        edges=np.unique(np.sort(np.r_[mesh.faces[:,[0,1]],mesh.faces[:,[1,2]],mesh.faces[:,[2,0]]],axis=1),axis=0)
        weights=np.linalg.norm(mesh.vertices[edges[:,0]]-mesh.vertices[edges[:,1]],axis=1)
        graph=coo_matrix((np.r_[weights,weights],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(mesh.vertices),)*2).tocsr()
        component_count,component_labels=connected_components(graph)
        distance=dijkstra(graph,directed=False,indices=root_index)
        finite=np.isfinite(distance)
        maximum=float(distance[finite].max())
        points=[]
        for fraction in [0.,.35,.7,1.]:
            band=finite & (abs(distance-fraction*maximum)<=.045*maximum)
            if not band.any():
                band[np.argmin(abs(distance-fraction*maximum))]=True
            point=mesh.vertices[band].mean(0)
            transformed=scale*point@rotation+translation
            measurements=[]
            for target_name,target_row,target,target_native,native_id in targets:
                sq,face,witness=igl.point_mesh_squared_distance(point[None],target.vertices,target.faces)
                normal=np.cross(target.triangles[face[0],1]-target.triangles[face[0],0],target.triangles[face[0],2]-target.triangles[face[0],0])
                normal/=np.linalg.norm(normal)
                delta=point-witness[0]
                native_sq,native_face,native_witness=igl.point_mesh_squared_distance(transformed[None],target_native.vertices,target_native.faces)
                measurements.append(dict(context=target_name,source_node_id=target_row['node_id'],source_face=int(face[0]),
                                         surface_distance_scaled_mm=float(np.sqrt(sq[0])*scale*1000),
                                         point_minus_nearest_surface_scaled_mm=(delta*scale*1000).tolist(),
                                         local_oriented_surface_normal_offset_scaled_mm=float(delta@normal*scale*1000),
                                         coordinate_side=dict(lateral='lateral' if delta[0]>0 else 'medial',vertical='superior' if delta[1]>0 else 'inferior',anteroposterior='anterior' if delta[2]>0 else 'posterior'),
                                         coordinate_side_qualification='Direction from the nearest actual source triangle witness; not a whole-muscle or closed-space classification',
                                         transformed_reference_to_native_surface_mm=float(np.sqrt(native_sq[0])*1000),native_part_id=native_id,native_face=int(native_face[0])))
            points.append(dict(geodesic_fraction=fraction,source_band_vertices=int(band.sum()),reference_point=point.tolist(),
                               coarse_transformed_reference_m=transformed.tolist(),registration_rms_uncertainty_mm=registration['rms_residual_mm'],
                               qualification='Coarse surface-distance-band centroid; no dense centreline was extracted or copied. Not a candidate waypoint.',relative_surfaces=measurements))
        contacts=[]
        for target_name,target_row,target,_,_ in targets:
            contact,summary=pair_contacts(mesh,target,0.,max_pairs=100000,max_tests=1000000)
            assert summary['complete_within_pair_budget']
            contacts.append(dict(context=target_name,source_node_id=target_row['node_id'],summary=summary,
                                 interpretation='Reference-only actual surface contacts; tonsillar constrictor perforation is a declared tissue interface'))
        branches.append(dict(name=name,cid=cid,source_node_id=row['node_id'],source_attr=row['attr'],source_components=int(component_count),
                             root_connected_vertices=int(finite.sum()),other_component_vertices=int((~finite).sum()),
                             parent_source_node_id=parent_row['node_id'],takeoff_source_vertex_index=root_index,
                             takeoff_reference_point=mesh.vertices[root_index].tolist(),parent_surface_face=int(parent_face[root_index]),
                             parent_surface_reference_point=parent_points[root_index].tolist(),closest_parent_surface_gap_scaled_mm=float(np.sqrt(d2[root_index])*scale*1000),
                             coarse_geodesic_extent_scaled_mm=maximum*scale*1000,coarse_landmarks=points,reference_context_triangle_contacts=contacts,
                             anatomical_acceptance=False))
    result=dict(schema='coarse-reference-relative-facial-corridor/v1',source_catalogue_sha256=SOURCE_SHA,source=reader.provenance(),
                registration_file='reference-native-similarity-v6.json',registration_rms_residual_mm=registration['rms_residual_mm'],
                branches=branches,relative_source_frame='PositiveX left, positiveY superior, positiveZ anterior, inferred from actual named bilateral and bone meshes',
                candidate_waypoints_created=False,dense_reference_route_copied=False,reference_surfaces_exported=False,
                native_route_domain=dict(method='Independently solve between fixed native muscle/vein surfaces, preserving reference-relative topology rather than transporting the reference centreline.',
                                         ascending_palatine=['Use a native facial takeoff beside the SG/SP interval, not the old posterior ECA-crossing shortcut.','Go over an anatomically reviewed superior constrictor border, maintaining native muscle/vein exclusion; broad palate perfusion remains unfinished.'],
                                         tonsillar=['Select the recipient-facing branch departure through the SG/medial-pterygoid interval; the mandibular ramus and retromandibular/maxillary veins remain fixed blockers.','Only superior-constrictor perforation into the native tonsil is explicitly permitted.'],
                                         glandular=['Use short native gland-groove departures; preserve nearby submental artery, facial vein and submandibular duct.'],
                                         direct_reference_registration_rejected=True,full_pharyngeal_lumen_volume_unverified=True),
                no_geometry_or_candidate_changed=True,elapsed_seconds=time.monotonic()-began)
    assert digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')==SOURCE_SHA
    atomic_json(file,result)
    print(json.dumps(dict(file=str(file.relative_to(ROOT)),branches=[dict(name=b['name'],source_node_id=b['source_node_id'],components=b['source_components'],coarse_extent_mm=b['coarse_geodesic_extent_scaled_mm'],takeoff_gap_mm=b['closest_parent_surface_gap_scaled_mm'],points=len(b['coarse_landmarks'])) for b in branches],elapsed_seconds=result['elapsed_seconds'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=int, choices=[1,2,3], required=True)
    args = parser.parse_args()
    {1:stage1,2:stage2,3:stage3}[args.stage]()
