"""Female partial tongue: six new core/IL carriers, five literal native muscles.

No installation. The selected core pose, source-bound basal IL field and its
explicit normal-only charts are the only new geometry inputs. Native female
extrinsics, their identities and their baseline defects remain unchanged.
"""
from __future__ import annotations
import argparse,copy,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'server')]
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.geometry_instance import instance_preimage
from geometry_tools.model import AtlasReader,Mesh
from tongue_basal_contract import require
_spec=importlib.util.spec_from_file_location('female_full_transfer',ROOT/'server/build-tongue-female-transfer-candidate.py')
base=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(base)
WORK=base.WORK
NATIVE=['LOCAL:female:genioglossus:239','LOCAL:female:l_hyoglossus:242','LOCAL:female:r_hyoglossus:255',
        'LOCAL:female:l_styloglossus:243','LOCAL:female:r_styloglossus:256']
SUPPRESS=[base.NATIVE_TONGUE,'HRA:VH_F_dorsal_tongue']
IL=['TA98TONGUE:male:intrinsic:inferior-longitudinal-left','TA98TONGUE:male:intrinsic:inferior-longitudinal-right']
LIMIT=('Partial female tongue installation: male-derived mucosa, aponeurosis and four intrinsic carriers; '
       'the five native female genioglossus/hyoglossus/styloglossus meshes are retained literally. '
       'Their existing attachment/crossing defects are not repaired. Intrinsic regions are approximate '
       'non-owning macroscopic envelopes, not exclusive histological tissues. Seventeen controls do not '
       'establish complete TA98 anatomy; missing terms, exact interdigitation and attachments remain unresolved.')


def read(p):return json.loads(Path(p).read_text())
def bound(p):return base.bound(p)


def selected_IL(male):
    fit_path=WORK/'female-basal-tongue-fit-v1/fit.json';fit=read(fit_path)
    spec=dict(matrix=fit['matrix'],protectedField=bound(fit_path))
    field,fieldproof=base.protected_overrides(spec,male)
    chart_stage=WORK/'female-basal-IL-normal-charts-v1';charted_fit=read(chart_stage/'fit.json')
    proof_path=chart_stage/'evidence.json';proof=read(proof_path)
    require(proof.get('only_normal_attribute_charts_changed') is True,'Explicit normal-only IL provenance required')
    rows={r['partId']:r for r in charted_fit['parts']};proofs={r['partId']:r for r in proof['proof']}
    result={}
    for pid in IL:
        row,pr=rows[pid],proofs[pid];source=field[pid]['mesh'];mesh=Mesh.load(row['mesh'])
        require(digest_file(row['mesh'])==row['sha256'] and pr['source_sha256']==digest_file(field[pid]['row']['mesh'])
                and pr['physical_triangles_literal'] is True and pr['all_changed_corners_below_protected_plane'] is True
                and pr['actual_geometry_audit']['passed'] is True and pr['backward_corner_normals']==0
                and np.array_equal(mesh.triangles,source.triangles),
                'IL normal chart changed actual physical triangles or lost its completed source checks')
        n=len(source.vertices);require(np.array_equal(mesh.vertices[:n],source.vertices),'Original IL controls must remain literal')
        changes={(r['face'],r['corner']):r for r in pr['changes']}
        require(len(changes)==pr['actual_normal_corners_repaired'],'IL chart provenance repeats corners')
        for f in range(len(source.faces)):
            for c in range(3):
                key=(f,c);v=int(mesh.faces[f,c]);old=int(source.faces[f,c])
                if key not in changes:
                    require(v==old and all(np.array_equal(mesh.attributes[k][v],a[old]) for k,a in source.attributes.items()),
                            'Unlisted IL corner or normal/UV field changed')
                else:
                    change=changes[key]
                    require(v==change['new_chart_vertex'] and old==change['source_vertex']
                            and np.array_equal(mesh.vertices[v],source.vertices[old])
                            and np.array_equal(mesh.attributes['normals'][v],change['new_word'])
                            and all(np.array_equal(mesh.attributes[k][v],a[old]) for k,a in source.attributes.items() if k!='normals'),
                            'IL corrected-index/normal provenance does not replay its literal source corner')
        result[pid]=dict(row=row,mesh=mesh)
    return result,dict(field_fit=bound(fit_path),field_source_proof=fieldproof,IL_normal_chart_proof=bound(proof_path))


def main(output):
    started=time.monotonic();male,female=base.load_sources();marker=male['finiteTongueIntegration'];reader=AtlasReader()
    axis_path=WORK/'female-axis-tongue-fit-v1/fit.json';axis=read(axis_path)
    rotation,scale,translation,quat=base.proper_similarity(axis['matrix'])
    require(axis['source_catalogue_sha256']==base.MALE_SHA and axis['target_catalogue_sha256']==base.FEMALE_SHA,
            'Selected body pose binds a different donor/target')
    il,provenance=selected_IL(male)
    require(output.resolve().is_relative_to(WORK) and not output.exists(),'Use a new partial female output stage')
    mp={p['id']:p for p in male['parts']};bp={p['id']:p for p in female['parts']}
    six=base.CORE|set(IL)
    mapping={pid:pid.replace('TA98TONGUE:male:','TA98TONGUE:female:') for pid in six}
    mapping.update({'ZA:Genioglossus muscle.l':NATIVE[0],'ZA:Genioglossus muscle.r':NATIVE[0],
        'ZA:Hyoglossus muscle.l':NATIVE[1],'ZA:Hyoglossus muscle.r':NATIVE[2],
        'TA98TONGUE:male:styloglossus:left':NATIVE[3],'TA98TONGUE:male:styloglossus:right':NATIVE[4]})
    concepts={c['id']:c['id'].replace('TA98TONGUE:male:','TA98TONGUE:female:') for c in male['concepts'] if c['id'].startswith('TA98TONGUE:male:')}
    ids={**mapping,**concepts};candidate=copy.deepcopy(female);parts=[];chunks={};materials={};bindings={}
    output.mkdir(parents=True);meshdir=output/'meshes';meshdir.mkdir();bank=bytearray();il_chunk=None;il_path=output/'IL-geometry.bin';encoding={}
    for pid in marker['partIds']:
        if pid not in six:continue
        source=mp[pid];part=base.rewrite(copy.deepcopy(source),ids)
        material=source['material']
        if material not in materials:
            mid=material.replace('TA98TONGUE:male:','TA98TONGUE:female:');require(mid not in candidate['materials'],'Partial material collision')
            candidate['materials'][mid]=copy.deepcopy(male['materials'][material]);materials[material]=mid
        part['material']=materials[material]
        if pid in IL:
            if il_chunk is None:
                il_chunk=len(candidate['chunks']);candidate['chunks'].append(dict(url='/local-models/'+str(il_path.relative_to(ROOT/'.local-models')),bytes=0))
            encoding[pid]=base.pack_override(bank,il[pid]['mesh'],part);part['chunk']=il_chunk
            part['transferMode']='protected-basal-field-with-source-bound-normal-only-charts'
        else:
            old_chunk=source['chunk']
            if old_chunk not in chunks:
                row=copy.deepcopy(male['chunks'][old_chunk]);path=reader.asset(row['url'])
                chunks[old_chunk]=len(candidate['chunks']);candidate['chunks'].append(row);bindings[str(path)]=digest_file(path)
            part.update(chunk=chunks[old_chunk],geometryInstance=dict(rotation=quat.tolist(),scale=scale),
                        geometryOrigin=((np.asarray(source.get('geometryOrigin',[0,0,0]))@rotation.T)*scale+translation).tolist(),
                        transferMode='literal-core-bank-instance')
        part.update(anatomicalAcceptance=False,tissueOwnershipClaim=False,coverageLimitation=LIMIT,
                    transferSource=dict(catalogue=bound(base.MALE),part_id=pid,axis_pose=bound(axis_path)))
        parts.append(part)
    require(len(parts)==6 and len(bank)<=200000,'Expected six new economical core/IL carriers')
    il_path.write_bytes(bank);candidate['chunks'][il_chunk]['bytes']=len(bank);bindings[str(il_path)]=digest_file(il_path)
    actual={};records=[]
    for part in parts:
        pid=next(k for k,v in mapping.items() if v==part['id']);mesh=reader.read_part(candidate,part,frame=base.FRAME)
        if pid in IL:
            source=il[pid]['mesh'];decoded=(source.vertices-base.OVERRIDE_ORIGIN).astype('<f4').astype(float)+base.OVERRIDE_ORIGIN
            require(np.array_equal(mesh.triangles,decoded[source.faces]) and instance_preimage(mesh) is None,
                    'Actual IL buffer differs or falsely claims a similarity preimage')
            used=np.asarray(mesh.source['source_vertex_ids'])
            require(all(np.array_equal(mesh.attributes[k],v[used]) for k,v in source.attributes.items()),'Actual IL buffer changed encoded normal/UV fields')
        else:
            donor=reader.read_part(male,mp[pid],frame='model:local-ta98:assembled-metres')
            require(instance_preimage(mesh) is not None and np.max(np.abs(mesh.triangles-((donor.triangles@rotation.T)*scale+translation)))<1e-12,
                    'Core is not the literal common proper-pose instance')
        raw=np.cross(mesh.triangles[:,1]-mesh.triangles[:,0],mesh.triangles[:,2]-mesh.triangles[:,0])
        word=mesh.attributes['normals'][mesh.faces].astype(float)*mesh.attributes['normalSigns'][mesh.faces,None]
        bad=int(np.count_nonzero(np.einsum('fci,fi->fc',word,raw)<=0))
        part['bounds']=[mesh.vertices.min(0).tolist(),mesh.vertices.max(0).tolist()]
        path=meshdir/(part['id'].replace(':','_')+'.npz');mesh.save(path);actual[part['id']]=mesh
        records.append(dict(part_id=part['id'],source_part_id=pid,mesh=bound(path),geometry_sha256=mesh.geometry_sha,
                            bank_sha256=mesh.source['bank_sha256'],nonpositive_corner_dots=bad,
                            geometry_instance_preimage_claimed=pid in base.CORE))
    for part in candidate['parts']:
        if part['id'] in SUPPRESS:part['suppressed']=True
    candidate['parts'].extend(parts);candidate['defaultHidden']=list(dict.fromkeys(candidate.get('defaultHidden',[])+SUPPRESS))
    groups=base.rewrite(copy.deepcopy(marker['termGroups']),ids)
    for g in groups:
        g['partIds']=list(dict.fromkeys(g['partIds']));g['anatomicalAcceptance']=False
        if g.get('anchor'):
            a=g['anchor'];a['point']=(np.asarray(a['barycentric'])@actual[a['partId']].triangles[a['face']]).tolist()
        for curve in g.get('curves',[]):curve['points']=((np.asarray(curve['points'])@rotation.T)*scale+translation).tolist()
    active=list(dict.fromkeys(mapping[p] for p in marker['partIds']));require(len(active)==11,'Expected six new plus five unique retained native meshes')
    newconcepts=[]
    for source in male['concepts']:
        if source['id'] not in concepts:continue
        c=base.rewrite(copy.deepcopy(source),ids);c['elements']=list(dict.fromkeys(c['elements']))
        c.update(coverageLimitation=LIMIT,anatomicalAcceptance=False);newconcepts.append(c)
    require(not(set(c['id'] for c in female['concepts'])&set(c['id'] for c in newconcepts)),'Partial concept collision')
    candidate['concepts'].extend(newconcepts)
    aliases={'LOCAL:female:tongue':active,'TA98FEMALE:A05.1.04.001':active,
             'HRA:VH_F_dorsal_tongue':[mapping['TA98TONGUE:male:mucosa']],
             'TA98FEMALE:A05.1.04.004':[mapping['TA98TONGUE:male:mucosa']]}
    for c in candidate['concepts']:
        if c['id'] in aliases:c.update(elements=aliases[c['id']],coverageLimitation=LIMIT,anatomicalAcceptance=False)
    carrier_hashes={r['part_id']:r['bank_sha256'] for r in records}
    for pid in NATIVE:
        mesh=reader.read_part(female,bp[pid],frame=base.FRAME);path=meshdir/(pid.replace(':','_')+'.npz');mesh.save(path)
        carrier_hashes[pid]=mesh.source['bank_sha256'];records.append(dict(part_id=pid,source_part_id=pid,mesh=bound(path),
            geometry_sha256=mesh.geometry_sha,bank_sha256=mesh.source['bank_sha256'],retained_native_literal=True))
    newmarker=base.rewrite(copy.deepcopy(marker),ids)
    newmarker.update(sex='female',partIds=active,termGroups=groups,activePhysicalMeshes=11,
        sourcePartId=base.NATIVE_TONGUE,sourceConceptId='LOCAL:female:tongue',sourceCatalogueSha256=base.FEMALE_SHA,
        bankSha256=digest_file(il_path),bankUrl=candidate['chunks'][il_chunk]['url'],carrierBankHashes=carrier_hashes,
        coverageLimitation=LIMIT,partialNativeExtrinsicsRetained=True,retainedNativePartIds=NATIVE,
        anatomicalAcceptance=False,tissueOwnershipClaim=False,stage='partial-local-review',
        transfer=dict(source=bound(base.MALE),axis_pose=bound(axis_path),new_physical_meshes=6,retained_native_meshes=5,
            additional_geometry_bytes=len(bank),unchanged_core_instances=4,protected_IL_buffers=2,
            target_context_checks_pending=True,actual_render_review_required=True,qualifying_candidate=False))
    candidate['finiteTongueIntegration']=newmarker;candidate['partsCount']=len(candidate['parts'])
    candidate['triangles']=sum(p['indexCount']//3 for p in candidate['parts'] if not p.get('suppressed'))
    after={p['id']:p for p in candidate['parts']}
    require(all(after[pid]==(dict(p,suppressed=True) if pid in SUPPRESS else p) for pid,p in bp.items()),'A retained native/source physical record changed')
    before_concepts={c['id']:c for c in female['concepts']};after_concepts={c['id']:c for c in candidate['concepts']}
    require(all(after_concepts[cid]==c for cid,c in before_concepts.items() if cid not in aliases),'An unrelated/native extrinsic concept changed')
    fixed=[]
    for pid in base.REQUIRED_FIXED:
        mesh=reader.read_part(female,bp[pid],frame=base.FRAME);path=meshdir/(pid.replace(':','_')+'.npz');mesh.save(path)
        fixed.append(dict(part_id=pid,mesh=bound(path),geometry_sha256=mesh.geometry_sha))
    out=output/'candidate-atlas.json';atomic_json(out,candidate)
    evidence=dict(schema='tongue-female-core-native-partial/v1',source=bound(base.MALE),target=bound(base.FEMALE),
        candidate=bound(out),axis_pose=bound(axis_path),IL_source_proof=provenance,new_part_ids=[p['id'] for p in parts],
        retained_native_part_ids=NATIVE,suppressed_source_parts=SUPPRESS,legacy_alias_elements=aliases,
        actual_world_meshes=records,fixed_context_meshes=fixed,part_mapping=mapping,source_bank_bindings=bindings,
        controls=17,physical_meshes=11,new_meshes=6,retained_meshes=5,typed_references=6,
        actual_new_corner_normal_gate_passed=all(r.get('nonpositive_corner_dots',0)==0 for r in records),
        encoded_IL_origin=base.OVERRIDE_ORIGIN.tolist(),IL_encoding_error_m=encoding,
        fixed_skeleton=True,native_extrinsics_literal=True,native_defects_not_solved=True,
        target_context_checks_pending=True,qualifying_candidate=False,anatomical_acceptance=False,installed=False,
        elapsed_seconds=time.monotonic()-started)
    atomic_json(output/'transfer-evidence.json',evidence)
    print(json.dumps(dict(candidate=str(out),sha256=digest_file(out),new_meshes=6,retained_meshes=5,physical_meshes=11,controls=17,
        added_bank_bytes=len(bank),new_corner_normals_passed=evidence['actual_new_corner_normal_gate_passed'],
        qualifying_candidate=False,installed=False,elapsed_seconds=time.monotonic()-started)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
