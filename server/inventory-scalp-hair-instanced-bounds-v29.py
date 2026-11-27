"""Conservative bounds of the actual instanced shafts; no contact testing."""
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
REVISION='31fab36e6897adaec9c8ddca201261607902203d75c64f3b214e57d329f460c2'


def main():
    assert digest_file(SOURCE)==EXPECTED
    atlas=json.loads(SOURCE.read_bytes());part=next(p for p in atlas['parts'] if p['id']=='ZA:Scalp hair:metric')
    runtime_path=BASE/'nearby-unresolved-participation-frame-v25.json';runtime=json.loads(runtime_path.read_bytes())
    assert runtime['source_catalogue_sha256']==EXPECTED and runtime['registry_revision']==REVISION
    status=next(r for r in runtime['rows'] if r['id']==part['id'])
    assert not status['normal_runtime_suppressed'] and not status['normal_default_visible_skin0'] and status['normal_visible_skin1']
    assert part.get('geometryScale',1)==1 and not part.get('geometryInstance') and not part.get('sourceOffset') and not part.get('geometryOrigin')
    reader=AtlasReader();field=part['hairField'];path=reader.asset(field['url']);assert digest_file(path)==field['sha256']
    bank=reader.asset(atlas['chunks'][part['chunk']]['url']);assert digest_file(bank)==status['registered_geometry_binding']['bank_sha256']
    assert digest_file(ROOT/'app/metric-hair-field.ts')==runtime['code_bindings']['app/metric-hair-field.ts']
    raw=np.memmap(path,dtype='u1',mode='r');roots=reader.view(raw,field['positions'],field['count'],'<f4',3).astype(float)
    # The normalized signed-short attribute follows WebGL's clamped conversion,
    # then hairVector() normalizes it before constructing an orthogonal basis.
    packed=reader.view(raw,field['normals'],field['count'],'<i2',3)
    directions=np.clip(packed.astype(float)/32767.,-1.,1.);lengths=np.linalg.norm(directions,axis=1)
    assert np.isfinite(roots).all() and np.isfinite(lengths).all() and np.all(lengths>0)
    directions/=lengths[:,None];tips=roots+directions*field['length'];radius=field['diameter']/2
    radial_extent=radius*np.sqrt(np.maximum(0.,1-directions**2))
    # A circular cylinder contains the rendered six-sided cylinder for every
    # hairVector basis. Add a conservative Float32 arithmetic envelope rather
    # than assuming the CPU's exact doubles equal the vertex shader's words.
    eps=np.finfo(np.float32).eps;gamma64=64*eps/(1-64*eps)
    halo=float(4*gamma64*(np.max(abs(roots))+field['length']+radius))
    lows=np.minimum(roots,tips)-radial_extent-halo;highs=np.maximum(roots,tips)+radial_extent+halo
    whole=np.array([lows.min(0),highs.max(0)]);branches=[];candidate_arrays={}
    def overlaps(low,high,other_low,other_high):return np.all(high>=other_low,axis=-1)&np.all(low<=other_high,axis=-1)
    for name in ['ascending-palatine','tonsillar','glandular-0','glandular-1']:
        source=BASE/(name+'-fields.npz')
        with np.load(source,allow_pickle=False) as z:
            values=z['fields'];assert str(z['frame'])=='model:local-ta98:assembled-metres'
            low=values.min(axis=(0,1));high=values.max(axis=(0,1))
        shaft_ids=np.flatnonzero(overlaps(lows,highs,low,high)).astype('<u4')
        root_ids=np.flatnonzero(np.all(roots>=low,axis=1)&np.all(roots<=high,axis=1)).astype('<u4')
        candidate_arrays[name+'_shaft_ids']=shaft_ids;candidate_arrays[name+'_root_ids']=root_ids
        branches.append(dict(branch=name,source_fields=str(source.relative_to(ROOT)),source_fields_sha256=digest_file(source),swept_bounds_m=[low.tolist(),high.tolist()],
            full_occupied_union_aabb_overlaps=bool(overlaps(whole[0],whole[1],low,high)),root_points_in_swept_aabb=int(len(root_ids)),
            conservative_instance_shaft_aabbs_overlapping=int(len(shaft_ids)),all_rendered_shafts_excluded_by_bounds=not len(shaft_ids),
            actual_contacts_tested=False,remaining_contact_test_required=bool(len(shaft_ids))))
    indices=BASE/'scalp-hair-instanced-bounds-v29-indices.npz';atomic_npz(indices,**candidate_arrays)
    result=dict(schema='scalp-hair-instanced-bounds/v1',source_catalogue_sha256=EXPECTED,registry_revision=REVISION,
        part_id=part['id'],world_frame='model:local-ta98:assembled-metres',source_runtime_evidence=str(runtime_path),source_runtime_sha256=digest_file(runtime_path),
        actual_hair_field=dict(path=str(path),sha256=field['sha256'],bytes=field['bytes'],count=field['count'],length_m=field['length'],diameter_m=field['diameter'],spacing_m=field['spacing'],
            root_attribute_dtype='Float32',normal_attribute_dtype='normalized signed Int16; clamp then normalize'),
        proxy_carrier_geometry_sha256=status['registered_geometry_binding']['geometry_sha256'],proxy_carrier_bank_sha256=status['registered_geometry_binding']['bank_sha256'],
        macro_carrier_is_not_occupied_volume_classifier=True,explicit_world_pose=status['explicit_record_pose'],
        runtime_prototype='metricHairGeometry: closed CylinderGeometry(diameter/2,diameter/2,length,6,1,false), translated by length/2; hairVector maps it to each exact root/normal.',
        source_code_bindings={name:runtime['code_bindings'][name] for name in ['app/metric-hair-field.ts','app/scene.tsx','app/depth-layers.ts','app/viewer-state.ts']},
        bounds_method='Per-instance closed circular cylinder envelopes, containing every six-sided rendered shaft, plus conservative Float32 gamma64 arithmetic padding. No carrier hull or synthetic anatomical volume.',
        float32_roundoff_halo_m=halo,normal_zero_count=0,root_bounds_m=[roots.min(0).tolist(),roots.max(0).tolist()],full_occupied_shaft_union_aabb_m=whole.tolist(),
        candidate_indices=str(indices),candidate_indices_sha256=digest_file(indices),branches=branches,
        skin_gating_proof=dict(normal_runtime_suppressed=False,local_default_skin0_visible=False,skin1_visible=True,is_skin_part=True,
            source='Exact ordinary-runtime loader and partVisible evaluation in nearby-unresolved-participation-frame-v25.json; no live visibility state was changed.'),
        no_actual_triangle_contact_or_occupied_intersection_test=True,no_geometry_or_visibility_mutation=True,anatomical_registration_not_inferred=True)
    target=BASE/'scalp-hair-instanced-bounds-v29.json';atomic_json(target,result);assert digest_file(SOURCE)==EXPECTED
    print(json.dumps(dict(output=str(target.relative_to(ROOT)),count=field['count'],halo_m=halo,branches=branches)))


if __name__=='__main__':main()
