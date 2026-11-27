"""Prepare an explicitly selected coupled tongue candidate; never install it.

--inputs is a tongue-coupled-inputs/v1 JSON with bound {path,sha256} values:
source_audit, inventory, whole, whole_fit, intrinsics, coverings,
extrinsics, styloglossus, styloglossus_checks, features, feature_ancestry,
covering_muscle_audit, retained_context_audit, selected_review.
covering_meshes separately binds mucosa/aponeurosis exports from the SAME bank.

No input is auto-selected.  Failed or absent coupled-interface checks prevent
even creating the output directory.  The current full mucosal shell is known
to fail these checks and must not be packed.  No catalogue or source is edited.
"""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'server'))
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.constraints import check_constraint
from geometry_tools.model import AtlasReader, Mesh
from tongue_packing_stages import (begin_geometry, seal_geometry, verify_stage,
                                  finalize_stage)
from tongue_whole_input import whole_input
from tongue_basal_contract import (APO, FRAME, INTRINSIC_KEYS, MODE, MUCO,
    bind_features, checked_binding, corner_normal_mesh, covering_muscle_gate,
    require, retained_envelope_gate, selected_extrinsics_gate)

_spec = importlib.util.spec_from_file_location('tongue_original_packer', ROOT/'server/build-tongue-local-candidate.py')
old = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(old)
WORK = old.WORK
WHOLE, MUSCLES = old.WHOLE, old.MUSCLES
REQUIRED = frozenset(('source_audit', 'inventory', 'whole', 'whole_fit', 'intrinsics',
    'coverings', 'extrinsics', 'styloglossus', 'styloglossus_checks', 'features',
    'feature_ancestry', 'covering_muscle_audit', 'retained_context_audit', 'selected_review'))
LIMIT = ('Partial local tongue proposal; exact term identity is not anatomical acceptance. '
    'Intrinsic ROIs are approximate, non-owning macroscopic muscle-extent envelopes. '
    'Transverse/vertical controls share one physical carrier. Native hyoglossus attachments, '
    'histological entheses, exact interdigitation and all absent terms remain pending. '
    'Selected GG/HG origin and insertion candidates are source-bound numerical proposals; '
    'their exact anatomical attachment extents and histological entheses remain unaccepted.')


def read(path):
    return json.loads(Path(path).read_text())


def source_record(path, sha, geometry=None):
    record = dict(path=str(Path(path).resolve()), sha256=sha)
    checked_binding(record, ROOT/'.local-models')
    if geometry is not None:
        record['geometry_sha256'] = geometry
    return record


def qualified_intrinsics(construction,path):
    """Adapt explicit shared-TV semantics without editing producer evidence."""
    result=copy.deepcopy(construction)
    for row in result.get('meshes',[]):
        terms=row.get('terms',[row['term']]if 'term'in row else [])
        if len(terms)<=1:
            continue
        shared=construction.get('shared_tv_representation',{})
        require(set(terms)=={'A05.1.04.108','A05.1.04.109'}
                and set(shared.get('terms',[]))==set(terms)
                and shared.get('one_mesh')is True and shared.get('duplicate_occupied_volumes')is False
                and shared.get('orientation_fields_are_illustrative')is True
                and row.get('representation')==MODE and row.get('exclusive_tissue_volume')is False
                and row.get('anatomical_acceptance')is False
                and construction.get('exclusive_tissue_partition_claimed')is False
                and construction.get('primary_reference')=='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
                'Shared-TV producer row lacks explicit qualified one-carrier anatomical evidence')
        if not row.get('coverageLimitation'):
            row['coverageLimitation']=('One approximate macroscopic transverse/vertical muscle-extent '
                'reference envelope, qualified by the selected primary anatomical study. Two term controls '
                'share this single physical carrier; directional fibre fields are illustrative. This is '
                'not two exclusive histological volumes or accepted donor-specific anatomy.')
            row['coverageQualificationSource']=dict(path=str(Path(path).resolve()),sha256=digest_file(path),
                primary_reference_url=construction['primary_reference'],representation_mode=MODE,
                explicit_shared_tv_representation=copy.deepcopy(shared))
    return result


def stored_extrinsic_replay(mesh,row):
    """Read the actual saved field names; do not coerce an unchecked bank."""
    path=checked_binding(dict(path=row['storedBank'],sha256=row['storedBankSha256']),ROOT/'.local-models')
    with np.load(path,allow_pickle=False)as stored:
        require('positions'in stored and 'geometryOrigin'in stored,'Stored extrinsic lacks explicit local coordinates')
        if 'faces'in stored and 'attribute_normals'in stored:
            faces,words=stored['faces'],stored['attribute_normals']
            schema='local-float32-faces-attribute-normals'
        else:
            require('indices'in stored and 'normals'in stored,'Unknown stored extrinsic array schema')
            faces,words=stored['indices'],stored['normals']
            schema='local-float32-indices-normals'
        require(stored['positions'].dtype==np.dtype('float32')and words.dtype==np.dtype('int16')
                and np.array_equal(stored['geometryOrigin'],np.array([0.,1.625,0.])),
                'Stored extrinsic must retain actual Float32 positions, Int16 words and reviewed dyadic origin')
        require(np.array_equal(stored['positions'].astype(float)+stored['geometryOrigin'],mesh.vertices)
                and np.array_equal(np.asarray(faces).reshape(-1,3),mesh.faces)
                and np.array_equal(words,mesh.attributes['normals']),
                'Selected extrinsic differs from actual reviewed stored coordinates/indices/normal words')
        if 'frame'in stored:
            require(str(stored['frame'])==mesh.frame,'Stored extrinsic frame changed')
        for attribute,field in [('uvs','uvs'if 'uvs'in stored else 'attribute_uvs'),
                                ('normalSigns','normalSigns'if 'normalSigns'in stored else 'attribute_normalSigns')]:
            if attribute in mesh.attributes:
                require(field in stored and np.array_equal(stored[field],mesh.attributes[attribute]),
                        'Selected extrinsic attribute differs from reviewed storage: '+attribute)
    return dict(path=str(path),sha256=digest_file(path),schema=schema,actual_stored_replay_passed=True)


def preflight(inputs_path):
    """Validate saved physical evidence only; no fabricated successful rows."""
    plan = read(inputs_path)
    require(plan.get('schema') == 'tongue-coupled-inputs/v1', 'Unknown coupled input schema')
    require(set(plan.get('bindings', {})) == REQUIRED, 'Explicitly bind every required coupled input')
    paths = {key:checked_binding(value, ROOT/'.local-models' if key != 'inventory' else ROOT)
             for key, value in plan['bindings'].items()}
    docs = {key:read(path) for key, path in paths.items() if path.suffix == '.json'}
    source, inventory = docs['source_audit'], docs['inventory']
    catalogue = Path(source['sourceCatalogues']['male']['path']).resolve()
    require(digest_file(catalogue) == source['sourceCatalogues']['male']['sha256'],
            'Current source catalogue has changed')
    native = read(catalogue)
    require(not native.get('finiteTongueIntegration'), 'A previously integrated tongue needs an explicit revision workflow')
    whole, _, whole_binding = whole_input(paths['whole'], paths['whole_fit'])
    intrinsic=qualified_intrinsics(docs['intrinsics'],paths['intrinsics'])
    docs['intrinsics']=intrinsic
    covers,stylo=docs['coverings'],docs['styloglossus']
    require(intrinsic.get('whole_input') == covers.get('whole_input') == whole_binding,
            'Regenerate all derived tissue from the SAME actual reviewed corrected whole')
    require(intrinsic.get('schema') == 'tongue-macroscopic-muscle-reference-envelopes/v1'
            and intrinsic.get('numerical_and_resource_passed') is True
            and intrinsic.get('intrinsic_pairs_crossing_free') is True
            and intrinsic.get('fixed_bone_crossings_zero') is True
            and intrinsic.get('fixed_bone_pairs_complete') is True
            and intrinsic.get('normal_chart_variant_passed') is True
            and intrinsic.get('exclusive_tissue_partition_claimed') is False,
            'Selected intrinsic construction lacks complete actual numerical/corner-chart evidence')
    require(covers.get('schema') in ('tongue-finite-coverings/v1','tongue-shared-exposed-coverings/v1')
            and covers.get('actual_emitted_checks', {}).get('passed') is True,
            'Selected coverings must have actual emitted finite-cell checks')
    require(inventory.get('term_count') == 36 and len(inventory.get('terms', [])) == 36,
            'Complete classified tongue inventory is required')
    source_sha = digest_file(catalogue)
    features = docs['features']
    require(features.get('wholeInput') == whole_binding
            and features.get('sourceCatalogueSha256') == source_sha
            and features.get('inventorySha256') == digest_file(paths['inventory'])
            and features.get('sourceWholeMeshSha256', features.get('carrierMeshSha256'))
                == digest_file(paths['whole']),
            'Typed feature definitions do not bind the actual selected corrected whole/inventory')
    require(stylo.get('source_catalogue_sha256') == source_sha
            and stylo.get('numerical_candidate_pass') is True
            and docs['styloglossus_checks'].get('passed') is True
            and docs['styloglossus_checks'].get('candidateSha256') == digest_file(paths['styloglossus']),
            'Selected styloglossus evidence is stale or incomplete')
    extrinsics = selected_extrinsics_gate(docs['extrinsics'],source_sha,ROOT/'.local-models')
    context = {r['partId']:r for r in source['contextParts']}
    muscle_sources, intrinsic_sources = {}, {}
    actual_cover_manifest_path = Path(covers['packed_manifest'])
    actual_cover_manifest = read(actual_cover_manifest_path)
    require({r['key'] for r in intrinsic['meshes']} == INTRINSIC_KEYS
            and len(intrinsic['meshes']) == 4, 'Exactly four unique physical intrinsic carriers are required')
    for row in intrinsic['meshes']:
        if row['key']=='superior-longitudinal' and actual_cover_manifest.get('part_ids',{}).get(row['key']):
            actual_sl_path=actual_cover_manifest_path.with_name(actual_cover_manifest['part_ids'][row['key']].replace(':','_')+'.npz')
            actual_sl=Mesh.load(actual_sl_path)
            require(actual_sl.geometry_sha==row['geometry_sha256'],
                    'Actual shared-bank SL does not retain the selected physical face order')
            intrinsic_sources[row['key']]=source_record(actual_sl_path,digest_file(actual_sl_path),actual_sl.geometry_sha)
        else:
            intrinsic_sources[row['key']] = source_record(row['mesh'], row['sha256'], row['geometry_sha256'])
        muscle_sources[row['key']] = intrinsic_sources[row['key']]
    retained = {}
    for pid,row in extrinsics.items():
        path = checked_binding(dict(path=row['mesh'], sha256=row['sha256']), ROOT/'.local-models')
        retained[pid] = source_record(path, row['sha256'], Mesh.load(path).geometry_sha)
        muscle_sources[pid] = retained[pid]
    for side, row in stylo['sides'].items():
        key = 'styloglossus-v2-'+side
        retained[key] = source_record(row['mesh'], row['sha256'], row['geometry_sha256'])
        muscle_sources[key] = retained[key]
    cover_sources = {key:dict(value) for key, value in plan.get('covering_meshes', {}).items()}
    require(set(cover_sources) == {'mucosa', 'aponeurosis'}, 'Bind both actual covering exports explicitly')
    for value in cover_sources.values():
        checked_binding(value, ROOT/'.local-models')
    # This gate intentionally precedes bank writes and rejects the current shell.
    shared_cover_interfaces = {}
    for row in plan.get('covering_muscle_shared_interfaces', []):
        cover, muscle = row['cover'], row['muscle']
        require(cover in cover_sources and muscle in muscle_sources,
                'Shared interface must name selected actual covering/muscle identities')
        a = Mesh.load(cover_sources[cover]['path'])
        muscle_source = muscle_sources[muscle]
        b = Mesh.load(muscle_source['path'])
        if 'normals'not in b.attributes:
            b=corner_normal_mesh(muscle_source['path'],muscle_source['geometry_sha256'])
        parameters = dict(first_faces=row['cover_faces'], second_faces=row['muscle_faces'])
        actual = check_constraint(dict(kind='shared-interface', parameters=parameters), a, b)
        require(actual['passed'], 'Actual selected covering/muscle seam lacks opposed facets/normals')
        key = (cover, muscle)
        require(key not in shared_cover_interfaces, 'Duplicate covering/muscle seam')
        shared_cover_interfaces[key] = dict(row, actual_check=actual)
    cover_gate = covering_muscle_gate(docs['covering_muscle_audit'], cover_sources, muscle_sources,
                                     shared_interfaces=shared_cover_interfaces)
    retained_gate = retained_envelope_gate(docs['retained_context_audit'], intrinsic_sources, retained)
    review = docs['selected_review']
    require(review.get('schema') == 'tongue-coupled-prepack-selection/v1'
            and review.get('selected_for_local_candidate') is True
            and review.get('exposed_mucosa_entry_boundaries_reviewed') is True
            and review.get('aponeurosis_muscle_interface_reviewed') is True
            and review.get('anatomical_acceptance') is False,
            'Actual surface/section selection must include exposed-entry and aponeurosis interfaces')
    expected = {key:value['sha256'] for key, value in plan['bindings'].items() if key != 'selected_review'}
    require(review.get('selected_input_sha256') == expected, 'Visual selection binds another physical input set')
    images = review.get('inspected_images', [])
    require(images and {'surface', 'section'} <= {r.get('kind') for r in images},
            'Actual inspected emitted surface and section images are required')
    for image in images:
        checked_binding(image, ROOT/'.local-models')
    return dict(plan=plan, paths=paths, docs=docs, original=native, catalogue=catalogue,
                source_sha=source_sha, source=source, whole=whole,
                cover_gate=cover_gate, retained_gate=retained_gate,
                shared_cover_interfaces=shared_cover_interfaces,extrinsics=extrinsics)


def intrinsic_storage_origin(row, mesh):
    """Honor only a bound dyadic intrinsic origin with literal decoded replay."""
    if 'geometryOrigin' not in row:
        return old.ORIGIN, False
    origin = np.asarray(row['geometryOrigin'], dtype=float)
    require(origin.shape == (3,) and np.array_equal(origin, [0., 1.625, 0.]),
            'Explicit intrinsic origin must be the reviewed dyadic origin')
    proof_path = checked_binding(row.get('storageOriginProof'), ROOT/'.local-models')
    proof = read(proof_path)
    require(proof.get('schema') == 'tongue-literal-intrinsic-origin-replay/v1'
            and proof.get('mesh') == dict(path=str(Path(row['mesh']).resolve()),
                                          sha256=row['sha256'],
                                          geometry_sha256=row['geometry_sha256'])
            and proof.get('geometryOrigin') == origin.tolist()
            and proof.get('literal_decoded_replay') is True
            and mesh.geometry_sha == row['geometry_sha256']
            and digest_file(row['mesh']) == row['sha256'],
            'Explicit intrinsic storage origin lacks exact source/hash evidence')
    decoded = (mesh.vertices-origin).astype('<f4').astype(float)+origin
    require(np.array_equal(decoded, mesh.vertices),
            'Explicit intrinsic origin changed literal physical coordinates')
    return origin, True


def pack_record(buffer, mesh, record, origin, *, literal=False):
    """Local Float32 storage; GG uses its reviewed dyadic origin, UVs literal."""
    result = copy.deepcopy(record)
    for key in ('positions', 'normals', 'normalSigns', 'indices', 'uvs', 'colors',
                'sourceOffset', 'geometryScale', 'geometryInstance'):
        result.pop(key, None)
    words, signs = old.normal_words(mesh)
    relative = (mesh.vertices-origin).astype('<f4')
    decoded = relative.astype(float)+origin
    if literal:
        require(np.array_equal(decoded, mesh.vertices),
                'Reviewed intrinsic storage must replay literal physical coordinates')
    require(np.linalg.norm(decoded-mesh.vertices, axis=1).max() <= 2e-8,
            'Actual position encoding exceeds the declared 20nm bound')
    result.update(vertexCount=len(relative), indexCount=mesh.faces.size,
                  geometryOrigin=origin.tolist(), positions=old.add_array(buffer, relative),
                  normals=old.add_array(buffer, words), normalSigns=old.add_array(buffer, signs),
                  indices=old.add_array(buffer, mesh.faces.astype('<u4').reshape(-1)),
                  bounds=[decoded[mesh.faces].reshape(-1, 3).min(0).tolist(),
                          decoded[mesh.faces].reshape(-1, 3).max(0).tolist()])
    if 'uvs' in mesh.attributes:
        uvs = np.asarray(mesh.attributes['uvs'])
        require(uvs.dtype == np.dtype('float32'), 'Copied native GG UVs must remain literal Float32')
        result['uvs'] = old.add_array(buffer, uvs.astype('<f4'))
    return result


def selected_extrinsic_material(part, mesh, materials):
    """The selected four producers lack UVs; never bind a UV image to them."""
    require(part['id'] in {'ZA:Genioglossus muscle.l','ZA:Genioglossus muscle.r',
                           'ZA:Hyoglossus muscle.l','ZA:Hyoglossus muscle.r'},
            'A UV-free correction is restricted to the four selected GG/HG replacements')
    require('uvs' not in mesh.attributes and 'uvs' not in part,
            'This correction only applies to the explicitly selected UV-free producer')
    family = 'genioglossus' if 'Genioglossus' in part['id'] else 'hyoglossus'
    mid = 'TA98TONGUE:male:material:selected-'+family
    materials[mid] = old.material([.74,.32,.34])
    previous = part['material']
    part.update(material=mid, materialCorrection=dict(
        reason='Selected corrected GG/HG producer has no UV buffer; use existing UV-free tissue shader',
        previousMaterial=previous, procedural=True, sampledTextures=False,
        shader='app/tissue-wall-material.ts', pattern='fibrous'))
    return part


def geometry_stage(inputs, output):
    start = time.monotonic()
    def remaining():
        if time.monotonic()-start >= 45:
            raise TimeoutError('Coupled tongue geometry exceeded its 45 second worker limit')
    output = output.resolve()
    require(output.is_relative_to(WORK) and output != WORK and not output.exists(),
            'Use a new local tongue stage; no overwrite is allowed')
    snapshot = begin_geometry(inputs)
    state = preflight(inputs)
    remaining()
    plan, paths, docs = state['plan'], state['paths'], state['docs']
    candidate = copy.deepcopy(state['original'])
    parts = {p['id']:p for p in candidate['parts']}
    source_rows = {r['partId']:r for r in state['source']['contextParts']}
    bank, records, meshes, origins, materials = bytearray(), [], {}, {}, {}
    chunk = len(candidate['chunks'])
    cover_manifest = read(docs['coverings']['packed_manifest'])
    shared_sl = cover_manifest.get('part_ids',{}).get('superior-longitudinal')
    require(shared_sl in (None,'TA98TONGUE:male:intrinsic:superior-longitudinal'),
            'Unknown shared-bank SL identity')
    intrinsic_ids = []
    for row in docs['intrinsics']['meshes']:
        remaining()
        pid = 'TA98TONGUE:male:intrinsic:'+row['key']
        intrinsic_ids.append(pid)
        if shared_sl and row['key']=='superior-longitudinal':
            continue
        mesh = corner_normal_mesh(row['mesh'], row['geometry_sha256'])
        old.bone_crossing_gate(mesh, source_rows)
        material_id = 'TA98TONGUE:male:material:'+row['key']
        materials[material_id] = old.material(row['color'])
        terms = old.intrinsic_terms(row)
        if len(terms) == 2:
            materials[material_id]['tissueAppearance']['fibres'] = dict(
                directions=[[1.,0.,0.],[0.,1.,0.]], weights=[.5,.5], control='tongue-tv')
        single = next((key for key, term in old.INTRINSIC_TERMS.items() if terms == [term]), None)
        physical_concept = ('TA98TONGUE:male:muscle:'+single if single
                            else 'TA98TONGUE:male:carrier:'+row['key'])
        template = old.append_mesh(bytearray(), mesh, part_id=pid,
            concept_id=physical_concept, name=row.get('name', row['key']),
            material=material_id, terms=terms, chunk=chunk,
            limitation=(row.get('coverageLimitation','')+' '+LIMIT).strip())
        if row.get('coverageQualificationSource'):
            template['coverageQualificationSource']=copy.deepcopy(row['coverageQualificationSource'])
        template.update(representationMode=MODE, tissueOwnershipClaim=False, roiBoundaries='approximate')
        if row.get('tissueDrawGroups'):
            template.update(tissueDrawGroups=copy.deepcopy(row['tissueDrawGroups']), pairedTissue=True)
        origin, literal_origin = intrinsic_storage_origin(row, mesh)
        part = pack_record(bank, mesh, template, origin, literal=literal_origin)
        records.append((part, mesh, Path(row['mesh']), row['sha256']))
    reader = AtlasReader()
    for side, row in docs['styloglossus']['sides'].items():
        remaining()
        source = Mesh.load(row['mesh'])
        bone_id = row['origin']['host_part_id']
        bone = reader.read_part(state['original'], parts[bone_id], frame=FRAME)
        mesh, constraint, chart = old.opposed_origin_chart(source, bone, row['origin'])
        pid = 'TA98TONGUE:male:styloglossus:'+side
        mid = 'TA98TONGUE:male:material:styloglossus'
        materials[mid] = old.material([.74,.32,.34])
        template = old.append_mesh(bytearray(), mesh, part_id=pid,
            concept_id='TA98TONGUE:male:muscle:styloglossus', name='Styloglossus ('+side+')',
            material=mid, terms=['A05.1.04.105'], chunk=chunk, limitation=LIMIT)
        records.append((pack_record(bank, mesh, template, old.ORIGIN), mesh, Path(row['mesh']), row['sha256']))
        origins[pid] = dict(host_part_id=bone_id, constraint=constraint, chart=chart)
    for pid,row in state['extrinsics'].items():
        remaining()
        with np.load(row['mesh'],allow_pickle=False)as raw:
            has_corners='corner_normals'in raw
        mesh = (corner_normal_mesh(row['mesh'],row['geometrySha256']) if has_corners
                else Mesh.load(row['mesh']))
        require('normals'in mesh.attributes and np.issubdtype(mesh.attributes['normals'].dtype,np.integer),
                'Selected extrinsic producer must provide actual stored normal words')
        old.bone_crossing_gate(mesh,source_rows)
        part = pack_record(bank, mesh, dict(parts[pid], chunk=chunk), np.array([0.,1.625,0.]))
        selected_extrinsic_material(part,mesh,materials)
        part.update(geometricAttachment=row.get('geometricAttachment','source-bound-origin-and-insertion-proposal'),
                    histologicalEnthesisComplete=False,reviewStatus='local corrected origin/insertion; anatomy pending',
                    attachmentEvidence=dict(path=row['attachmentReportPath'],sha256=row['attachmentReportSha256']))
        if row.get('storedBank'):
            part['storedSourceReplay']=stored_extrinsic_replay(mesh,row)
        records.append((part, mesh, Path(row['mesh']), row['sha256']))
    require(len(records) == (9 if shared_sl else 10),
            'Expected selected standalone muscles plus four GG/HG replacements; SL is counted once')
    require(len(bank) <= 2_000_000 and sum(len(mesh.faces) for _,mesh,_,_ in records) <= 30_000,
            'Coupled muscle bank exceeds the declared 2MB / 30,000-triangle resource limit')
    require(cover_manifest.get('schema') == 'shared-finite-tissues/v1'
            and cover_manifest.get('passed') is True
            and cover_manifest.get('reciprocal_normals_use_same_words_with_opposite_sign') is True,
            'Finite covering bank must retain its actual reciprocal normal fields')
    cover_bank = checked_binding(dict(path=cover_manifest['bank'], sha256=cover_manifest['bank_sha256']), ROOT/'.local-models')
    cover_records = []
    covering_items=[('mucosa','A05.1.04.011',[.76,.39,.40]),
                    ('aponeurosis','A05.1.04.025',[.76,.66,.61])]
    if shared_sl:
        covering_items.append(('superior-longitudinal','A05.1.04.106',[.72,.30,.29]))
    for key, term, color in covering_items:
        part = copy.deepcopy(cover_manifest['parts'][key])
        material_id = 'TA98TONGUE:male:material:'+key
        materials[material_id] = old.material(color)
        concept_id=('TA98TONGUE:male:muscle:'+key if key=='superior-longitudinal'
                    else 'TA98TONGUE:male:'+key)
        part.update(conceptId=concept_id, name=key.capitalize(), system='muscular' if key=='superior-longitudinal' else 'digestive',
            chunk=chunk+1, material=material_id, ta98Term=term, sectionSolid=True, pairedTissue=True,
            license='CC BY-SA 4.0', provenance=copy.deepcopy(records[0][0]['provenance']),
            coverageLimitation=LIMIT, anatomicalAcceptance=False)
        if key=='superior-longitudinal':
            part.update(representationMode=MODE,tissueOwnershipClaim=False,roiBoundaries='approximate')
        cover_records.append(part)
    # All unsafe interface gates and the actual visual selection already passed.
    output.mkdir(parents=True)
    bank_path = output/'geometry.bin'
    bank_path.write_bytes(bank)
    url = '/local-models/'+str(bank_path.relative_to(ROOT/'.local-models'))
    cover_url = '/local-models/'+str(cover_bank.relative_to(ROOT/'.local-models'))
    candidate['chunks'].extend([dict(url=url, bytes=len(bank)),
                                dict(url=cover_url, bytes=cover_bank.stat().st_size)])
    candidate['materials'].update(materials)
    for part, mesh, _, _ in records:
        if part['id'] in parts:
            parts[part['id']].clear()
            parts[part['id']].update(part)
        else:
            candidate['parts'].append(part)
    candidate['parts'].extend(cover_records)
    parts = {p['id']:p for p in candidate['parts']}
    parts['ZA:Tongue']['suppressed'] = True
    candidate['defaultHidden'] = list(dict.fromkeys(candidate.get('defaultHidden', [])+['ZA:Tongue']))
    decoded_rows = []
    for part, expected, path, sha in records:
        remaining()
        actual = reader.read_part(candidate, part, frame=FRAME)
        ids = np.asarray(actual.source['source_vertex_ids'])
        require(np.array_equal(actual.triangles, expected.triangles),
                'Actual emitted storage changed a reviewed physical triangle')
        words, signs = old.normal_words(expected)
        require(np.array_equal(actual.attributes['normals'], words[ids])
                and np.array_equal(actual.attributes['normalSigns'], signs[ids]),
                'Actual emitted selected normal words or signs changed')
        if 'uvs' in expected.attributes:
            require(np.array_equal(actual.attributes['uvs'], expected.attributes['uvs'][ids]),
                    'Actual emitted native GG UVs changed')
        topo = topology(actual)
        require(topo['closed_oriented_manifold'], 'Actual emitted muscle is not closed/oriented')
        meshes[part['id']] = actual
        item = dict(part_id=part['id'], input_npz=str(path), input_sha256=sha,
                    decoded_geometry_sha256=actual.geometry_sha, normal_words_and_signs_literal=True,
                    topology=topo)
        if part['id'] in origins:
            origin = origins[part['id']]
            bone = reader.read_part(candidate, parts[origin['host_part_id']], frame=FRAME)
            proof = check_constraint(origin['constraint'], bone, actual)
            require(proof['passed'] and proof['literal_opposed_normal_words'],
                    'Actual emitted SG lost its literal opposed bone seam')
            item['styloid_interface'] = dict(origin, actual_emitted_check=proof)
        decoded_rows.append(item)
    for part in cover_records:
        remaining()
        actual = reader.read_part(candidate, part, frame=FRAME)
        key = ('mucosa' if part['id']==MUCO else 'aponeurosis' if part['id']==APO else 'superior-longitudinal')
        if key=='superior-longitudinal':
            selected=next(r for r in docs['intrinsics']['meshes']if r['key']==key)
            expected_path,expected_sha=selected['mesh'],selected['sha256']
        else:
            expected_path,expected_sha=plan['covering_meshes'][key]['path'],plan['covering_meshes'][key]['sha256']
        expected = Mesh.load(expected_path)
        require(actual.geometry_sha == expected.geometry_sha
                and np.array_equal(actual.triangles, expected.triangles),
                'Appended shared covering bank differs from its actually audited export')
        meshes[part['id']] = actual
        decoded_rows.append(dict(part_id=part['id'], input_npz=expected_path,
            input_sha256=expected_sha, decoded_geometry_sha256=actual.geometry_sha,
            topology=topology(actual), literal_shared_bank_retained=True))
    for pid in intrinsic_ids:
        require(check_constraint(dict(kind='nested', parameters={}), meshes[pid], state['whole'])['passed'],
                'Actual emitted ROI escaped the corrected reviewed whole')
    from itertools import combinations
    for a, b in combinations(intrinsic_ids, 2):
        remaining()
        require(check_constraint(dict(kind='disjoint', parameters=dict(occupancy='solid')),
                                 meshes[a], meshes[b])['passed'],
                'Actual emitted intrinsic tissue carriers are not mutually separated')
    # Every reciprocal runtime draw group needs its literal emitted seam proof.
    active = set(meshes)
    for pid in active:
        for group in parts[pid].get('tissueDrawGroups', []):
            peer = group.get('neighbor')
            if not peer:
                continue
            require(peer in active, 'Reciprocal tissue neighbor is absent from active assembly')
            first_faces = [f for g in parts[pid]['tissueDrawGroups'] if g.get('neighbor') == peer
                           for f in range(g['start']//3, (g['start']+g['count'])//3)]
            second_faces = [f for g in parts[peer].get('tissueDrawGroups', []) if g.get('neighbor') == pid
                            for f in range(g['start']//3, (g['start']+g['count'])//3)]
            require(second_faces, 'Reciprocal tissue draw group is absent')
            proof = check_constraint(dict(kind='shared-interface', parameters=dict(
                first_faces=first_faces, second_faces=second_faces)), meshes[pid], meshes[peer])
            require(proof['passed'], 'Actual shared-bank charts lost opposed coordinates/words')
    feature_groups = bind_features(docs['features'], docs['feature_ancestry'], state['whole'],
                                   meshes[MUCO], digest_file(cover_bank))
    part_ids = [MUCO, APO]+intrinsic_ids+[r[0]['id'] for r in records if r[0]['id'].startswith('TA98TONGUE:male:styloglossus:')]
    part_ids += old.NATIVE_MUSCLES['genioglossus'][1]+old.NATIVE_MUSCLES['hyoglossus'][1]
    require(len(part_ids) == len(set(part_ids)) == 12, 'Twelve actual unique active meshes are required')
    groups = [dict(key='tongue',name='Tongue',term='A05.1.04.001',conceptId=WHOLE,partIds=part_ids,kind='exact'),
        dict(key='muscles',parent='tongue',name='Muscles of tongue',term='A04.1.05.001',conceptId=MUSCLES,
             partIds=part_ids[2:],kind='exact')]
    for key, (term, ids) in old.NATIVE_MUSCLES.items():
        groups.append(dict(key=key,parent='muscles',name=old.NAMES[key],term=term,
            conceptId='TA98TONGUE:male:muscle:'+key,partIds=ids,kind='exact'))
    for key, term in {**old.INTRINSIC_TERMS, 'styloglossus':'A05.1.04.105'}.items():
        ids = [pid for pid in part_ids if term in parts[pid].get('ta98Terms', [parts[pid].get('ta98Term')])]
        groups.append(dict(key=key,parent='muscles',name=old.NAMES[key],term=term,
            conceptId='TA98TONGUE:male:muscle:'+key,partIds=ids,kind='exact'))
    for key, term in [('mucosa','A05.1.04.011'),('aponeurosis','A05.1.04.025')]:
        groups.append(dict(key=key,parent='tongue',name=key.capitalize(),term=term,
            conceptId='TA98TONGUE:male:'+key,partIds=['TA98TONGUE:male:'+key],kind='exact'))
    groups.extend(feature_groups)
    require(len(groups) == len({g['term'] for g in groups}) == 17, 'Expected seventeen identity-distinct controls')
    existing = {c['id']:c for c in candidate['concepts']}
    require(not any(g['conceptId'] in existing for g in groups), 'New tongue concepts collide with current catalogue')
    existing['ZA:Tongue']['elements'] = part_ids.copy()
    existing['ZA:Tongue']['coverageLimitation'] = LIMIT
    by_key = {g['key']:g for g in groups}
    for pid in intrinsic_ids:
        part = parts[pid]
        if len(part.get('ta98Terms', [])) > 1:
            candidate['concepts'].append(dict(id=part['conceptId'],name=part['name'],elements=[pid],
                sharedFibreTerritory=True,coverageLimitation=LIMIT,anatomicalAcceptance=False))
    for group in groups:
        concept = dict(id=group['conceptId'],name=group['name'],elements=group['partIds'],
            ta98Term=group['term'],ta98Kind='exact',coverageLimitation=LIMIT,anatomicalAcceptance=False)
        if group.get('parent'):
            concept['partOf'] = by_key[group['parent']]['conceptId']
        candidate['concepts'].append(concept)
    marker = dict(schema='tongue-local/v1',sex='male',wholeConceptId=WHOLE,partIds=part_ids,
        termGroups=groups,bankSha256=digest_file(bank_path),bankUrl=url,
        carrierBankHashes={MUCO:digest_file(cover_bank)},sourcePartId='ZA:Tongue',
        sourceCatalogueSha256=state['source_sha'],representationMode=MODE,tissueOwnershipClaim=False,
        roiBoundaries='approximate',stage='partial-local-review',anatomicalAcceptance=False,
        coverageLimitation=LIMIT,legacyWholeAlias=True,activePhysicalMeshes=12,
        sharedBankPartIds=[part['id']for part in cover_records])
    candidate['finiteTongueIntegration'] = marker
    candidate['triangles'] = sum(p['indexCount']//3 for p in candidate['parts'] if not p.get('suppressed'))
    inventory = copy.deepcopy(docs['inventory'])
    inventory.update(schema='ta98-tongue-local-control-inventory/v1', controls=17,physicalMeshes=12,
        newPhysicalMeshes=8,installed=False,anatomicalAcceptance=False)
    by_term = {g['term']:g for g in groups}
    for row in inventory['terms']:
        group = by_term.get(row['authority']['id'])
        row['candidate_control'] = copy.deepcopy(group) if group else None
        row['candidate_control_absent'] = group is None
        row['distinct_representation_verified'] = False
    evidence = dict(schema='tongue-coupled-pending-geometry/v1',source_catalogue=str(state['catalogue']),
        source_catalogue_sha256=state['source_sha'],records=decoded_rows,
        intrinsic_part_ids=intrinsic_ids,styloid_interfaces=origins,
        output_bank_sha256=digest_file(bank_path),shared_covering_bank_sha256=digest_file(cover_bank),
        new_meshes=8,assembly_meshes=12,represented_controls=17,unrepresented_terms=19,
        original_source_tongue_geometry_literal_and_suppressed=True,legacy_whole_alias=True,
        extrinsic_attachment_reports={pid:row['attachmentReport']for pid,row in state['extrinsics'].items()},
        extrinsic_histological_attachments='unaccepted',representation_mode=MODE,
        tissue_ownership_claim=False,anatomical_acceptance=False,all_input_gates_passed=True,
        emitted_render_review_required=True,installed=False,elapsed_seconds=time.monotonic()-start)
    remaining()
    atomic_json(output/'pending-inventory.json', inventory)
    atomic_json(output/'decoded-packing-proof.json', evidence)
    atomic_json(output/'pending-atlas.json', candidate)
    checkpoint = seal_geometry(output,snapshot,remaining)
    remaining()
    return dict(output=str(output),new_meshes=8,assembly_meshes=12,controls=17,
        installed=False,candidate_published=False,emitted_render_review_required=True,
        verification_batches=checkpoint['verification_batches'],elapsed_seconds=time.monotonic()-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('geometry','verify','finalize'), required=True,
                        help='One saved bounded stage; never continue automatically after a timeout')
    parser.add_argument('--inputs', type=Path, help='Explicit hash-bound input JSON for the geometry stage')
    parser.add_argument('--output', type=Path, required=True, help='New geometry stage or its exact saved checkpoint')
    parser.add_argument('--batch', type=int, help='Explicit verification batch from checkpoint.json')
    args = parser.parse_args()
    require(args.output.resolve().is_relative_to(WORK) and args.output.resolve()!=WORK,
            'Use a local tongue workbench stage')
    if args.stage=='geometry':
        require(args.inputs is not None and args.batch is None, 'Geometry requires --inputs, without --batch')
        result=geometry_stage(args.inputs,args.output)
    elif args.stage=='verify':
        require(args.batch is not None and args.inputs is None, 'Verify requires --batch, without --inputs')
        result=verify_stage(args.output,args.batch)
    else:
        require(args.inputs is None and args.batch is None, 'Finalize reads only the complete saved checkpoint')
        result=finalize_stage(args.output)
    print(json.dumps(result))
