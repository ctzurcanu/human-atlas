"""Prepare a female placement of the installed partial male tongue revision.

No registration is inferred and no atlas is installed. One explicitly supplied
proper similarity places all twelve physical carriers together, reusing the
two male geometry banks literally. Actual target-context and render review are
separate required inputs to the guarded install preparation stage.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'server')]
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.geometry_instance import instance_pose,instance_preimage
from geometry_tools.model import AtlasReader,Mesh
from tongue_basal_contract import require

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
MALE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FEMALE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
MALE_SHA='1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39'
FEMALE_SHA='770f47f802c79cf6f3e4336bf2fa7c42cf710998069221f533a1ea3f8322e470'
FRAME='model:local-female-ta98:assembled-metres'
NATIVE_TONGUE='LOCAL:female:tongue:1089'
SUPPRESS=[NATIVE_TONGUE,'LOCAL:female:genioglossus:239','LOCAL:female:l_hyoglossus:242',
          'LOCAL:female:r_hyoglossus:255','LOCAL:female:l_styloglossus:243',
          'LOCAL:female:r_styloglossus:256','HRA:VH_F_dorsal_tongue']
REQUIRED_FIXED=['LOCAL:female:mandible:960','LOCAL:female:hyoid:1002',
                'LOCAL:female:skull_cranium:999','LOCAL:female:skull_cranium:1000',
                'LOCAL:female:skull_cranium:1001','LOCAL:female:mylohyoid:238',
                'LOCAL:female:l_geniohyoid:241','LOCAL:female:r_geniohyoid:254',
                'LOCAL:female:uvula:1088','LOCAL:female:pharynx:1087','LOCAL:female:nasopharynx:1047']
LIMIT=('Transferred partial male tongue revision with 17 controls and 12 physical meshes. '
       'Intrinsic muscle extents and typed regional references remain approximate, non-owning representations. '
       'Target female origins/insertions are placement proposals; no literal female cortical attachment, '
       'histological enthesis, exact interdigitation or complete TA98 coverage is claimed.')
CORE={'TA98TONGUE:male:mucosa','TA98TONGUE:male:aponeurosis',
      'TA98TONGUE:male:intrinsic:superior-longitudinal','TA98TONGUE:male:intrinsic:transverse-vertical-shared'}
OVERRIDE_ORIGIN=np.array([0.,1.5,0.])


def read(path):return json.loads(Path(path).read_text())
def bound(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=digest_file(p))


def load_sources():
    require(digest_file(MALE)==MALE_SHA and digest_file(FEMALE)==FEMALE_SHA,'An explicitly selected installed catalogue changed')
    male,female=read(MALE),read(FEMALE);marker=male.get('finiteTongueIntegration',{})
    require(marker.get('schema')=='tongue-local/v1' and marker.get('stage')=='partial-local-review'
            and marker.get('anatomicalAcceptance') is False and marker.get('tissueOwnershipClaim') is False
            and len(marker.get('termGroups',[]))==17 and len(set(marker.get('partIds',[])))==12,
            'Use the explicit installed partial male revision')
    require(not female.get('finiteTongueIntegration'),'Existing female tongue integration requires an explicit revision')
    ids={p['id'] for p in female['parts']};require(set(SUPPRESS+REQUIRED_FIXED)<=ids,'A female native source/context ID changed')
    return male,female


def mapping(male,female):
    marker=male['finiteTongueIntegration']
    aliases={}
    for group in marker['termGroups']:
        aliases[group['conceptId']]=group['conceptId'].replace('TA98TONGUE:male:','TA98TONGUE:female:')
    aliases['TA98TONGUE:male:carrier:transverse-vertical-shared']='TA98TONGUE:female:carrier:transverse-vertical-shared'
    parts={}
    for pid in marker['partIds']:
        if pid.startswith('TA98TONGUE:male:'):parts[pid]=pid.replace('TA98TONGUE:male:','TA98TONGUE:female:')
        else:
            family='genioglossus' if 'Genioglossus' in pid else 'hyoglossus'
            side='left' if pid.endswith('.l') else 'right'
            parts[pid]='TA98TONGUE:female:'+family+':'+side
            aliases[pid]=parts[pid]
    return parts,aliases


def rewrite(value,ids):
    if isinstance(value,str):return ids.get(value,value)
    if isinstance(value,list):return [rewrite(v,ids) for v in value]
    if isinstance(value,dict):return {ids.get(k,k):rewrite(v,ids) for k,v in value.items()}
    return value


def alias_targets(marker,part_map):
    groups={g['key']:[part_map[p] for p in g['partIds']] for g in marker['termGroups']}
    return {
        'LOCAL:female:tongue':groups['tongue'],'TA98FEMALE:A05.1.04.001':groups['tongue'],
        'LOCAL:female:genioglossus':groups['genioglossus'],'TA98FEMALE:A05.1.04.101':groups['genioglossus'],
        'LOCAL:female:l_hyoglossus':[part_map['ZA:Hyoglossus muscle.l']],
        'LOCAL:female:r_hyoglossus':[part_map['ZA:Hyoglossus muscle.r']],
        'TA98FEMALE:A05.1.04.102':groups['hyoglossus'],
        'LOCAL:female:l_styloglossus':[part_map['TA98TONGUE:male:styloglossus:left']],
        'LOCAL:female:r_styloglossus':[part_map['TA98TONGUE:male:styloglossus:right']],
        'TA98FEMALE:A05.1.04.105':groups['styloglossus'],
        'HRA:VH_F_dorsal_tongue':groups['dorsum'],'TA98FEMALE:A05.1.04.004':groups['dorsum']}


def proper_similarity(matrix):
    m=np.asarray(matrix,dtype=float)
    require(m.shape==(4,4) and np.isfinite(m).all() and np.array_equal(m[3],[0.,0.,0.,1.]),'Use one finite affine 4x4 matrix')
    a=m[:3,:3];scale=float(np.linalg.norm(a,axis=0).mean());require(scale>0,'Positive uniform scale required')
    rotation=a/scale
    require(np.allclose(rotation.T@rotation,np.eye(3),rtol=0,atol=1e-10)
            and abs(np.linalg.det(rotation)-1)<1e-10,'Only a proper rotation and uniform scale are allowed; no reflection/shear')
    quat=Rotation.from_matrix(rotation).as_quat();quat/=np.linalg.norm(quat)
    # Use exactly the same reconstructed rotation as runtime instance_pose.
    runtime_rotation,_=instance_pose(dict(geometryInstance=dict(rotation=quat.tolist(),scale=scale)))
    require(np.max(np.abs(runtime_rotation-rotation))<1e-12,'Quaternion conversion changed the declared proper pose')
    return runtime_rotation,scale,m[:3,3].copy(),quat


def protected_overrides(spec,male):
    """Replay the one explicitly authorized source-bound protected basal field."""
    if not spec.get('protectedField'):return {},None
    binding=spec['protectedField'];require(digest_file(binding['path'])==binding['sha256'],'Protected field proof changed')
    fit=read(binding['path']);contract_path=Path(fit['basal_contract']);contract=read(contract_path)
    require(fit.get('schema')=='installed-tongue-female-protected-basal-fit/v1'
            and digest_file(contract_path)==fit['basal_contract_sha256']
            and fit.get('source_catalogue_sha256')==MALE_SHA and fit.get('target_catalogue_sha256')==FEMALE_SHA
            and fit.get('skeleton_modified') is False and fit.get('anatomical_acceptance') is False,
            'Protected overrides lack exact donor/target/field proof')
    require(contract.get('schema')=='female-tongue-protected-basal-contract/v1'
            and contract.get('fixed_skeleton') is True
            and set(contract.get('protected_core_carriers',[]))==CORE
            and digest_file(contract['source_fit'])==contract['source_fit_sha256'],
            'Protected core or exact source pose changed')
    prior=read(contract['source_fit'])
    require(np.array_equal(np.asarray(spec['matrix']),np.asarray(prior['matrix']))
            and np.array_equal(np.asarray(fit['matrix']),np.asarray(prior['matrix'])),
            'Protected buffers must follow the explicitly selected common source pose')
    for path,sha in fit.get('source_bindings',{}).items():require(digest_file(path)==sha,'Protected field source changed: '+path)
    low,high=float(contract['low_Y_m']),float(contract['protected_Y_m']);delta=np.asarray(contract['delta_m'],float)
    lower=1-1.5*delta[1]/(high-low)
    require(high>low and delta.shape==(3,) and np.isfinite(delta).all() and delta[1]>0
            and lower>0 and abs(lower-contract['analytic_determinant_lower_bound'])<1e-12
            and contract.get('globally_invertible') is True,
            'Protected field lacks its actual positive analytic Jacobian bound')
    rows={r['partId']:r for r in fit['parts']};ids=set(male['finiteTongueIntegration']['partIds'])
    require(set(rows)==ids and len(rows)==12,'The source-bound field must account for the whole twelve-part closure')
    override={};proofs=[]
    for pid,row in rows.items():
        require(digest_file(row['mesh'])==row['sha256'] and digest_file(row['pre_basal_mesh'])==row['pre_basal_sha256']
                and digest_file(row['source_mesh'])==row['source_sha256'],'Protected mesh/index ancestry changed')
        before,after,donor=(Mesh.load(row[k]) for k in ('pre_basal_mesh','mesh','source_mesh'))
        require(after.frame==before.frame==FRAME and after.geometry_sha==row['geometry_sha256']
                and donor.geometry_sha==row['source_geometry_sha256']
                and np.array_equal(after.faces,before.faces) and np.array_equal(after.faces,donor.faces),
                'Protected field must preserve literal source face indices')
        u=np.clip((before.vertices[:,1]-low)/(high-low),0.,1.);w=1-u*u*(3-2*u);changed=w>0
        expected=before.vertices+w[:,None]*delta
        jac=np.broadcast_to(np.eye(3),(len(expected),3,3)).copy()
        jac[:,:,1]+=(-6*u*(1-u)/(high-low))[:,None]*delta[None,:]
        require(np.array_equal(after.vertices,expected) and float(np.linalg.det(jac).min())>0,
                'Protected coordinates do not replay the declared finite source field')
        expected_attrs={k:np.asarray(v).copy() for k,v in before.attributes.items()}
        if np.any(changed):
            raw=expected_attrs['normals'];n=np.linalg.solve(jac[changed].transpose(0,2,1),raw[changed].astype(float)[...,None])[...,0]
            n/=np.linalg.norm(n,axis=1)[:,None];raw[changed]=np.rint(n*32767).astype('<i2')
        require(set(after.attributes)==set(expected_attrs) and all(np.array_equal(after.attributes[k],v) for k,v in expected_attrs.items()),
                'Protected normal field, signs, UVs or other source attributes changed without declared provenance')
        if pid in CORE:
            require(not np.any(changed) and np.array_equal(after.vertices,before.vertices)
                    and after.attributes_sha==before.attributes_sha,'One of four protected core carriers changed')
        else:
            override[pid]=dict(mesh=after,row=row)
        proofs.append(dict(part_id=pid,source=bound(row['source_mesh']),pre_field=bound(row['pre_basal_mesh']),
                           field_output=bound(row['mesh']),face_indices_literal=True,UVs_and_other_fields_literal=True,
                           inverse_transpose_normal_words_replayed=True,minimum_sampled_jacobian=float(np.linalg.det(jac).min()),
                           protected_core_literal=pid in CORE))
    require(set(override)==ids-CORE and len(override)==8,'Only two IL and six extrinsic buffers may be overridden')
    return override,dict(schema='tongue-female-protected-override-source-proof/v1',field_fit=binding,
        field_contract=bound(contract_path),analytic_jacobian_lower_bound=lower,meshes=proofs,
        actual_decoded_geometry_audits_pending=True,anatomical_acceptance=False,qualifying_candidate=False)


def pack_override(bank,mesh,part):
    """One actual Float32 bank; the decoded result is a new audit candidate."""
    def put(values):
        bank.extend(b'\0'*((-len(bank))%4));offset=len(bank);bank.extend(np.asarray(values).tobytes());return offset
    require(mesh.attributes['normals'].dtype==np.dtype('int16') and 'normalSigns' in mesh.attributes
            and np.isin(mesh.attributes['normalSigns'],[-1,1]).all(),'Overrides require actual encoded words/signs')
    relative=(mesh.vertices-OVERRIDE_ORIGIN).astype('<f4');decoded=relative.astype(float)+OVERRIDE_ORIGIN
    error=float(np.max(np.linalg.norm(decoded-mesh.vertices,axis=1)))
    require(error<=2e-8,'Override Float32 encoding exceeds the declared20nm bound')
    for key in ('geometryInstance','geometryScale','sourceOffset','geometryOrigin','positions','normals','normalSigns','indices','uvs','colors'):
        part.pop(key,None)
    part.update(geometryOrigin=OVERRIDE_ORIGIN.tolist(),vertexCount=len(relative),indexCount=mesh.faces.size,
                positions=put(relative),normals=put(mesh.attributes['normals'].astype('<i2')),
                normalSigns=put(mesh.attributes['normalSigns'].astype('i1')),indices=put(mesh.faces.astype('<u4').reshape(-1)))
    for key,dtype in [('uvs','<f4'),('colors','u1')]:
        if key in mesh.attributes:
            require(mesh.attributes[key].dtype==np.dtype(dtype),'Override source attribute must retain its storage dtype')
            part[key]=put(mesh.attributes[key])
    return error


def readiness(output):
    male,female=load_sources();part_map,concept_map=mapping(male,female);aliases=alias_targets(male['finiteTongueIntegration'],part_map)
    concepts={c['id'] for c in female['concepts']};require(set(aliases)<=concepts,'Legacy female concept aliases changed')
    report=dict(schema='tongue-female-transfer-control-readiness/v1',source=bound(MALE),target=bound(FEMALE),
        target_frame=FRAME,part_mapping=part_map,concept_mapping=concept_map,
        existing_alias_elements=aliases,suppressed_source_parts=SUPPRESS,required_fixed_context=REQUIRED_FIXED,
        controls=17,physical_meshes=12,typed_references=6,additional_geometry_bytes=0,
        sourcePartId=NATIVE_TONGUE,sourceConceptId='LOCAL:female:tongue',
        geometry_generated=False,fit_inferred=False,installed=False,anatomical_acceptance=False,
        placement_input_schema='tongue-female-proper-similarity/v1',
        target_audit_schema='tongue-female-target-context-audit/v1')
    if output:
        require(not output.exists(),'Use a new readiness stage');output.mkdir(parents=True);atomic_json(output/'readiness.json',report)
    print(json.dumps(report))


def prepare(spec_path,output):
    began=time.monotonic();male,female=load_sources();spec=read(spec_path)
    require(spec.get('schema')=='tongue-female-proper-similarity/v1' and spec.get('source')==bound(MALE)
            and spec.get('target')==bound(FEMALE),'Placement must bind the exact installed male and original female catalogues')
    require(spec.get('qualification') and spec.get('fixed_skeleton') is True,'Declare the placement rationale and fixed skeleton')
    rotation,scale,translation,quat=proper_similarity(spec['matrix'])
    overrides,override_proof=protected_overrides(spec,male)
    require(output.resolve().is_relative_to(WORK) and not output.exists(),'Use a new local female placement stage')
    part_map,concept_map=mapping(male,female);all_ids={**part_map,**concept_map}
    before_parts={p['id']:p for p in female['parts']};male_parts={p['id']:p for p in male['parts']}
    candidate=copy.deepcopy(female);marker=male['finiteTongueIntegration'];reader=AtlasReader()
    material_map={};chunk_map={};bank_bindings={};actual={};decoded=[];new_parts=[]
    output.mkdir(parents=True);meshdir=output/'meshes';meshdir.mkdir()
    override_bank=bytearray();override_errors={};override_chunk=None;override_bank_path=output/'override-geometry.bin'
    for pid in marker['partIds']:
        source=male_parts[pid];part=rewrite(copy.deepcopy(source),all_ids);old_chunk=source['chunk']
        if old_chunk not in chunk_map:
            chunk=copy.deepcopy(male['chunks'][old_chunk]);path=reader.asset(chunk['url']);sha=digest_file(path)
            require(path.is_file() and path.stat().st_size==chunk['bytes'],'A literal source bank changed size')
            chunk_map[old_chunk]=len(candidate['chunks']);candidate['chunks'].append(chunk);bank_bindings[str(path)]=sha
        part['chunk']=chunk_map[old_chunk]
        mid=source['material']
        if mid not in material_map:
            new_mid=mid.replace('TA98TONGUE:male:','TA98TONGUE:female:')
            require(new_mid not in candidate['materials'],'Female transfer material identity collision')
            material_map[mid]=new_mid;candidate['materials'][new_mid]=copy.deepcopy(male['materials'][mid])
        part['material']=material_map[mid]
        previous_pose=instance_pose(source);r0,s0=previous_pose if previous_pose is not None else (np.eye(3),1.)
        total_rotation=rotation@r0;total_quat=Rotation.from_matrix(total_rotation).as_quat();total_quat/=np.linalg.norm(total_quat)
        part['geometryInstance']=dict(rotation=total_quat.tolist(),scale=scale*s0)
        part['geometryOrigin']=((np.asarray(source.get('geometryOrigin',[0.,0.,0.]))@rotation.T)*scale+translation).tolist()
        if pid in overrides:
            if override_chunk is None:
                override_chunk=len(candidate['chunks'])
                candidate['chunks'].append(dict(url='/local-models/'+str(override_bank_path.relative_to(ROOT/'.local-models')),bytes=0))
            override_errors[pid]=pack_override(override_bank,overrides[pid]['mesh'],part)
            part['chunk']=override_chunk
        donor_attachment={k:copy.deepcopy(v) for k,v in source.items()
            if k!='geometryOrigin' and any(word in k.lower() for word in ('attachment','styloid','contact','storedsourcereplay','histological'))}
        for key in donor_attachment:part.pop(key,None)
        part.update(geometricAttachment='similarity-transferred-approximate-origin-insertion',
                    histologicalEnthesisComplete=False,anatomicalAcceptance=False,
                    coverageLimitation=(source.get('coverageLimitation','')+' '+LIMIT).strip(),
                    transferProvenance=dict(sourceCatalogue=bound(MALE),sourcePartId=pid,
                        sourceAttachmentEvidence=donor_attachment,
                        donorCorticalChartsAreNotTargetSharedInterfaces=True,
                        targetAttachmentStatus='actual target context audit and review required'))
        if pid in overrides:
            part['geometricAttachment']='source-bound-protected-field-origin-insertion-proposal'
            part['transferProvenance'].update(protectedField=copy.deepcopy(spec['protectedField']),
                sourceIndexAndUVFieldsLiteral=True,sourceNormalWordsTransportedByInverseTranspose=True,
                nonAffineOverrideHasNoGeometryInstancePreimage=True)
        new_parts.append(part)
    if overrides:
        require(len(override_bank)<=2000000 and sum(len(r['mesh'].faces) for r in overrides.values())<=30000,
                'Protected override bank exceeds the declared2MB/30000-triangle limit')
        override_bank_path.write_bytes(override_bank);candidate['chunks'][override_chunk]['bytes']=len(override_bank)
    for part in new_parts:
        pid=next(k for k,v in part_map.items() if v==part['id']);source=male_parts[pid]
        world=reader.read_part(dict(chunks=candidate['chunks']),part,frame=FRAME)
        donor=reader.read_part(male,source,frame='model:local-ta98:assembled-metres')
        if pid in overrides:
            expected=overrides[pid]['mesh'];wanted=(expected.vertices-OVERRIDE_ORIGIN).astype('<f4').astype(float)+OVERRIDE_ORIGIN
            require(np.array_equal(world.triangles,wanted[expected.faces]) and instance_preimage(world) is None,
                    'Actual non-affine decoded triangles differ or incorrectly claim a similarity preimage')
            require(np.array_equal(world.faces,expected.faces),'Protected override changed literal source face indices')
            used=np.asarray(world.source['source_vertex_ids'])
            require(all(np.array_equal(world.attributes[k],v[used]) for k,v in expected.attributes.items()),
                    'Actual override normal words/signs, UVs or other fields changed during packing')
            raw=np.cross(world.triangles[:,1]-world.triangles[:,0],world.triangles[:,2]-world.triangles[:,0])
            effective=world.attributes['normals'][world.faces].astype(float)*world.attributes['normalSigns'][world.faces,None]
            bad=int(np.count_nonzero(np.einsum('fci,fi->fc',effective,raw)<=0))
            override_errors[pid]=dict(position_error_m=override_errors[pid],nonpositive_decoded_corner_dots=bad,
                source_mesh=bound(overrides[pid]['row']['mesh']),source_face_indices_literal=True,
                actual_decoded_attributes_literal=True,geometry_instance_preimage_claimed=False)
        else:
            wanted=(donor.vertices@rotation.T)*scale+translation
            require(np.max(np.linalg.norm(world.vertices-wanted,axis=1))<1e-12,'World coordinates do not replay the explicit similarity')
            require(instance_preimage(world) is not None,'Literal geometryInstance preimage did not replay')
            require(np.array_equal(world.faces,donor.faces),'A shared-bank face chart changed')
            sign=donor.attributes.get('normalSigns',np.ones(len(donor.vertices),dtype='i1'))
            require(np.array_equal(world.attributes.get('normalSigns',np.ones(len(world.vertices),dtype='i1')),sign)
                    and np.allclose(world.attributes['normals'],donor.attributes['normals'].astype(float)@rotation.T,rtol=0,atol=1e-9),
                    'Raw signs or actual transformed normal field changed')
        part['bounds']=[world.vertices.min(0).tolist(),world.vertices.max(0).tolist()]
        mesh_path=meshdir/(part['id'].replace(':','_')+'.npz');world.save(mesh_path)
        actual[part['id']]=world
        decoded.append(dict(part_id=part['id'],source_part_id=pid,mesh=bound(mesh_path),
                            geometry_sha256=world.geometry_sha,source_bank_sha256=world.source['bank_sha256'],
                            proper_similarity_world_replay=pid not in overrides,literal_shared_bank=pid not in overrides,
                            protected_field_override=pid in overrides,geometry_instance_preimage_claimed=pid not in overrides))
    override_evidence_path=None
    if overrides:
        override_proof.update(actual_decoded_buffers=override_errors,emitted_bank=bound(override_bank_path),
            encoded_origin=OVERRIDE_ORIGIN.tolist(),all_decoded_corner_normals_positive=all(r['nonpositive_decoded_corner_dots']==0 for r in override_errors.values()),
            actual_self_bone_covering_checks_pending=True,qualifying_candidate=False)
        override_evidence_path=output/'protected-override-proof.json';atomic_json(override_evidence_path,override_proof)
    require(len(new_parts)==len(set(p['id'] for p in new_parts))==12,'Transfer physical identity duplication')
    for part in candidate['parts']:
        if part['id'] in SUPPRESS:part['suppressed']=True
    candidate['parts'].extend(new_parts);candidate['defaultHidden']=list(dict.fromkeys(candidate.get('defaultHidden',[])+SUPPRESS))
    groups=rewrite(copy.deepcopy(marker['termGroups']),all_ids)
    for group in groups:
        if group.get('anchor'):
            anchor=group['anchor'];mesh=actual[anchor['partId']]
            transformed=(np.asarray(anchor['point'])@rotation.T)*scale+translation
            point=np.asarray(anchor['barycentric'])@mesh.triangles[anchor['face']]
            require(np.linalg.norm(point-transformed)<1e-12,'Transferred anchor no longer agrees with its literal face/barycentric coordinates')
            anchor['point']=point.tolist()
        for curve in group.get('curves',[]):curve['points']=((np.asarray(curve['points'])@rotation.T)*scale+translation).tolist()
        group['anatomicalAcceptance']=False
    newconcepts=[]
    for source in male['concepts']:
        if source['id'] not in concept_map:continue
        concept=rewrite(copy.deepcopy(source),all_ids)
        concept.update(anatomicalAcceptance=False,coverageLimitation=(source.get('coverageLimitation','')+' '+LIMIT).strip())
        newconcepts.append(concept)
    require(not ({c['id'] for c in candidate['concepts']} & {c['id'] for c in newconcepts}),'Female concept collision')
    candidate['concepts'].extend(newconcepts)
    aliases=alias_targets(marker,part_map)
    for concept in candidate['concepts']:
        if concept['id'] in aliases:concept.update(elements=aliases[concept['id']],coverageLimitation=LIMIT,anatomicalAcceptance=False)
    new_marker=rewrite(copy.deepcopy(marker),all_ids)
    new_marker.update(sex='female',sourcePartId=NATIVE_TONGUE,sourceConceptId='LOCAL:female:tongue',
        sourceCatalogueSha256=FEMALE_SHA,sourceMaleCatalogue=bound(MALE),termGroups=groups,
        stage='partial-local-review',anatomicalAcceptance=False,tissueOwnershipClaim=False,
        coverageLimitation=LIMIT,transfer=dict(schema='tongue-female-similarity-transfer/v1',placement=bound(spec_path),
            scale=scale,rotation=quat.tolist(),translation_m=translation.tolist(),
            fixed_skeleton=True,literal_source_banks=not bool(overrides),additional_geometry_bytes=len(override_bank),
            target_context_checks_pending=True,target_render_review_required=True,
            donor_cortical_contact_not_target_contact=True))
    new_marker['carrierBankHashes']={row['part_id']:row['source_bank_sha256'] for row in decoded}
    if overrides:
        new_marker.update(bankSha256=digest_file(override_bank_path),bankUrl=candidate['chunks'][override_chunk]['url'])
        new_marker['transfer'].update(protectedField=copy.deepcopy(spec['protectedField']),
            protectedOverrideProof=bound(override_evidence_path),literal_core_instance_count=4,
            non_affine_override_count=8,decoded_geometry_checks_pending=True,qualifying_candidate=False)
    candidate['finiteTongueIntegration']=new_marker
    candidate['partsCount']=len(candidate['parts'])
    candidate['triangles']=sum(p['indexCount']//3 for p in candidate['parts'] if not p.get('suppressed'))
    after_parts={p['id']:p for p in candidate['parts']}
    require(all(after_parts[pid]==(dict(part,suppressed=True) if pid in SUPPRESS else part) for pid,part in before_parts.items()),
            'A native female physical record changed beyond explicit old-source suppression')
    require(all(candidate['materials'][k]==v for k,v in female['materials'].items()),'An unrelated female material changed')
    require(candidate['chunks'][:len(female['chunks'])]==female['chunks'],'A native female bank changed')
    fixed=[]
    for pid in REQUIRED_FIXED:
        mesh=reader.read_part(female,before_parts[pid],frame=FRAME);mp=meshdir/(pid.replace(':','_')+'.npz');mesh.save(mp)
        fixed.append(dict(part_id=pid,mesh=bound(mp),geometry_sha256=mesh.geometry_sha))
    evidence=dict(schema='tongue-female-transfer-geometry/v1',source=bound(MALE),target=bound(FEMALE),placement=bound(spec_path),
        active_part_ids=new_marker['partIds'],part_mapping=part_map,concept_mapping=concept_map,
        native_alias_elements=aliases,suppressed_source_parts=SUPPRESS,required_fixed_context=REQUIRED_FIXED,
        actual_world_meshes=decoded,fixed_context_meshes=fixed,source_bank_bindings=bank_bindings,
        controls=17,physical_meshes=12,typed_references=6,geometry_bytes_added=len(override_bank),
        fixed_skeleton=True,shared_charts_and_raw_bank_normals_literal=not bool(overrides),
        shared_core_charts_and_raw_bank_normals_literal=True,
        literal_core_instance_count=4 if overrides else 12,protected_override_count=len(overrides),
        protected_override_proof=bound(override_evidence_path) if override_evidence_path else None,
        emitted_bank_bindings={str(override_bank_path):digest_file(override_bank_path)} if overrides else {},
        target_context_checks_pending=True,actual_render_review_required=True,
        anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atlas_path=output/'candidate-atlas.json';atomic_json(atlas_path,candidate)
    evidence['candidate']=bound(atlas_path);atomic_json(output/'transfer-evidence.json',evidence)
    atomic_json(output/'placement.json',spec)
    print(json.dumps(dict(candidate=str(atlas_path.resolve()),candidate_sha256=digest_file(atlas_path),
                         evidence=str((output/'transfer-evidence.json').resolve()),controls=17,meshes=12,
                         geometry_bytes_added=len(override_bank),qualifying_candidate=False,
                         decoded_geometry_checks_pending=True,installed=False,elapsed_seconds=time.monotonic()-began)))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stage',choices=['readiness','prepare'],required=True)
    ap.add_argument('--output',type=Path);ap.add_argument('--spec',type=Path);args=ap.parse_args()
    if args.stage=='readiness':readiness(args.output)
    else:
        require(args.output is not None and args.spec is not None,'Prepare requires explicit --spec and new --output')
        prepare(args.spec.resolve(),args.output.resolve())
