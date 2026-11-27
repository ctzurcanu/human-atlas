"""Actual 20 covering pairs +24 non-owning muscle-envelope diagnostics.

The only allowed covering coplanarity is the literal declared APO/SL chart.
Writes bounded witnesses and a deliberately review-pending input manifest;
never packs, installs, changes normals or implies anatomical acceptance.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.constraints import check_constraint
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from tongue_basal_contract import MODE,covering_muscle_gate,retained_envelope_gate,selected_extrinsics_gate,require

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
REFERENCE='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/'


def read(path):
    return json.loads(Path(path).read_text())


def bind(path):
    path=Path(path).resolve();return dict(path=str(path),sha256=digest_file(path))


def source_record(key,path,mesh):
    path=Path(path).resolve()
    with np.load(path,allow_pickle=False)as raw:
        if 'corner_normals'in raw:
            words=raw['corner_normals'].copy();signs=np.ones(mesh.faces.shape,dtype='i1')
        else:
            words=mesh.attributes['normals'][mesh.faces]
            signs=mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),dtype='i1'))[mesh.faces]
    require(words.dtype==np.dtype('int16')and words.shape==(*mesh.faces.shape,3)
            and np.all(np.linalg.norm(words.astype(float),axis=2)>0),
            'Actual raw Int16 corner words are required: '+key)
    require(np.isin(signs,[-1,1]).all(),'Actual normal signs must remain literal ±1')
    return dict(key=key,path=str(path),mesh=str(path),sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,
        oriented_surface_sha256=mesh.surface_sha,frame=mesh.frame,triangles=len(mesh.faces),
        raw_int16_corner_words_sha256=hashlib.sha256(words.tobytes()).hexdigest(),
        raw_corner_signs_sha256=hashlib.sha256(signs.tobytes()).hexdigest(),
        raw_words_normalized_or_reencoded=False)


def main(coverings,intrinsics,extrinsics,styloglossus,features,ancestry,output):
    began=time.monotonic()
    def remaining():
        if time.monotonic()-began>50:
            raise TimeoutError('Selected coupling audit exceeded its 50 second worker budget')
    coverings,intrinsics,extrinsics,styloglossus,features,ancestry,output=(p.resolve()for p in
        (coverings,intrinsics,extrinsics,styloglossus,features,ancestry,output))
    require(all(p.is_relative_to(WORK)and p.is_file()for p in
        (coverings,intrinsics,extrinsics,styloglossus,features,ancestry)), 'Use explicit existing local inputs')
    require(output.is_relative_to(WORK)and output!=WORK and not output.exists(),'Use a new local audit stage')
    cover,intrinsic,selected,stylo=(read(p)for p in (coverings,intrinsics,extrinsics,styloglossus))
    source_sha=cover['whole_input']['source_catalogue_sha256']
    selected_rows=selected_extrinsics_gate(selected,source_sha,ROOT/'.local-models')
    require(cover['passed']is True and intrinsic['numerical_and_resource_passed']is True
            and stylo['numerical_candidate_pass']is True,'Selected source numerical evidence is incomplete')
    manifest_path=Path(cover['packed_manifest']);manifest=read(manifest_path)
    require(digest_file(manifest['bank'])==manifest['bank_sha256'],'Actual shared bank changed')
    mesh_paths={key:manifest_path.with_name(pid.replace(':','_')+'.npz')for key,pid in manifest['part_ids'].items()}
    require(set(mesh_paths)=={'mucosa','aponeurosis','superior-longitudinal'},'Use the actual three-cell bank')
    for row in intrinsic['meshes']:
        if row['key']!='superior-longitudinal':
            require(digest_file(row['mesh'])==row['sha256'],'Actual intrinsic NPZ changed')
            mesh_paths[row['key']]=Path(row['mesh'])
    for pid,row in selected_rows.items():
        mesh_paths[pid]=Path(row['mesh'])
    for side,row in stylo['sides'].items():
        require(digest_file(row['mesh'])==row['sha256'],'Actual SG source changed')
        mesh_paths['styloglossus-v2-'+side]=Path(row['mesh'])
    meshes={key:Mesh.load(path)for key,path in mesh_paths.items()}
    sources={key:source_record(key,mesh_paths[key],mesh)for key,mesh in meshes.items()}
    trees={key:BVH(mesh)for key,mesh in meshes.items()}
    intrinsic_keys=[r['key']for r in intrinsic['meshes']]
    retained_keys=list(selected_rows)+['styloglossus-v2-left','styloglossus-v2-right']
    muscle_keys=intrinsic_keys+retained_keys
    require(len(intrinsic_keys)==4 and len(retained_keys)==6 and len(muscle_keys)==10,
            'Require four current intrinsic and six actual retained muscles')
    def shared_faces(first,second):
        return [i for g in manifest['parts'][first]['tissueDrawGroups']
                if g.get('neighbor')==manifest['part_ids'][second]
                for i in range(g['start']//3,(g['start']+g['count'])//3)]
    apo_faces,sl_faces=shared_faces('aponeurosis','superior-longitudinal'),shared_faces('superior-longitudinal','aponeurosis')
    shared_check=check_constraint(dict(kind='shared-interface',parameters=dict(first_faces=apo_faces,second_faces=sl_faces)),
                                  meshes['aponeurosis'],meshes['superior-longitudinal'])
    require(len(apo_faces)==len(sl_faces)==122 and shared_check['passed']
            and shared_check['literal_opposed_normal_words'],'Only the actual literal122-face APO/SL chart may be coplanar')
    shared={('aponeurosis','superior-longitudinal'):dict(cover_faces=apo_faces,muscle_faces=sl_faces,actual_check=shared_check)}
    output.mkdir(parents=True);(output/'covering-muscle').mkdir();(output/'retained-context').mkdir()
    measured=[]
    def measure(first,second,folder):
        remaining()
        arrays,summary=pair_contacts(meshes[first],meshes[second],0.,first_bvh=trees[first],second_bvh=trees[second],
                                     max_pairs=100000,max_tests=1000000)
        require(summary['complete_within_pair_budget']and not summary['truncated'],
                'Actual source pair exceeded its finite measurement budget')
        filename=(first+'--'+second).replace(':','_').replace('/','_')+'.npz'
        path=folder/filename;atomic_npz(path,**arrays)
        measured.append(dict(first=first,second=second,summary=summary,witnesses=str(path),witnesses_sha256=digest_file(path)))
        atomic_json(output/'progress.json',dict(schema='tongue-selected-coupling-progress/v1',
            completed_pairs=len(measured),expected_pairs=44,complete=False,measurements=measured,
            elapsed_seconds=time.monotonic()-began,installed=False),overwrite=True)
        return arrays,summary,path
    cover_pairs=[]
    for first in ('mucosa','aponeurosis'):
        for second in muscle_keys:
            arrays,summary,path=measure(first,second,output/'covering-muscle')
            overlap=arrays['pairs'][arrays['relation']=='coplanar-overlap'].tolist()
            cover_pairs.append(dict(cover=first,muscle=second,cover_path=str(mesh_paths[first]),
                cover_sha256=sources[first]['sha256'],muscle_path=str(mesh_paths[second]),muscle_sha256=sources[second]['sha256'],
                summary=summary,coplanar_face_pairs=overlap,witnesses=str(path),witnesses_sha256=digest_file(path)))
    covering_audit=dict(schema='tongue-finite-covering-muscle-audit/v1',passed=False,pairs=cover_pairs,
        source_catalogue_sha256=source_sha,sources=list(sources.values()),bank=manifest['bank'],bank_sha256=manifest['bank_sha256'],
        source_construction=str(coverings),source_construction_sha256=digest_file(coverings),
        selected_extrinsics=str(extrinsics),selected_extrinsics_sha256=digest_file(extrinsics),
        literal_shared_interface=shared_check,anatomical_acceptance=False)
    covering_audit['passed']=True
    try:
        covering_muscle_gate(covering_audit,{k:sources[k]for k in ('mucosa','aponeurosis')},
            {k:sources[k]for k in muscle_keys},shared_interfaces=shared)
    except ValueError as exc:
        covering_audit.update(passed=False,failure=str(exc))
    cover_path=output/'covering-muscle/summary.json';atomic_json(cover_path,covering_audit)
    retained_pairs=[]
    for first in intrinsic_keys:
        for second in retained_keys:
            arrays,summary,path=measure(first,second,output/'retained-context')
            overlap=summary['crossing_face_pairs']+summary['coplanar_overlap_face_pairs']
            retained_pairs.append(dict(intrinsic=first,retained=second,
                intrinsic_sha256=sources[first]['sha256'],retained_sha256=sources[second]['sha256'],
                intrinsic_geometry_sha256=sources[first]['geometry_sha256'],retained_geometry_sha256=sources[second]['geometry_sha256'],
                summary=summary,witnesses=str(path),witnesses_sha256=digest_file(path),
                classification='allowed-macroscopic-envelope-overlap'if overlap else 'no-measured-overlap',
                primary_reference_url=REFERENCE,tissue_ownership_claim=False))
    retained_audit=dict(schema='tongue-intrinsic-retained-context-audit/v1',passed=True,
        passed_scope='Complete source-bound measurements only; not exclusive tissue ownership or anatomical acceptance',
        representation_mode=MODE,tissue_ownership_claim=False,roi_boundaries='approximate',
        primary_reference_url=REFERENCE,complete_within_pair_budget=True,pair_count=24,expected_pair_count=24,
        intrinsic_construction=str(intrinsics),construction_sha256=digest_file(intrinsics),
        source_catalogue_sha256=source_sha,selected_extrinsics=str(extrinsics),selected_extrinsics_sha256=digest_file(extrinsics),
        intrinsic=[sources[k]for k in intrinsic_keys],context=[sources[k]for k in retained_keys],pairs=retained_pairs,
        crossing_face_pairs=sum(p['summary']['crossing_face_pairs']for p in retained_pairs),
        coplanar_overlap_face_pairs=sum(p['summary']['coplanar_overlap_face_pairs']for p in retained_pairs),
        anatomical_acceptance=False)
    retained_envelope_gate(retained_audit,{k:sources[k]for k in intrinsic_keys},{k:sources[k]for k in retained_keys})
    retained_path=output/'retained-context/summary.json';atomic_json(retained_path,retained_audit)
    features_doc=read(features);whole=features_doc['wholeInput']
    bindings={key:bind(path)for key,path in dict(source_audit=whole['source_audit_path'],
        inventory=ROOT/'reports/ta98-tongue-inventory.json',whole=whole['whole_mesh'],whole_fit=whole['fit_result_path'],
        intrinsics=intrinsics,coverings=coverings,extrinsics=extrinsics,styloglossus=styloglossus,
        styloglossus_checks=styloglossus.with_name('context-checks.json'),features=features,feature_ancestry=ancestry,
        covering_muscle_audit=cover_path,retained_context_audit=retained_path).items()}
    bindings['selected_review']=dict(path=str(output/'selected-review.json'),sha256=None,pending=True)
    pending=dict(schema='tongue-coupled-inputs/v1',bindings=bindings,
        covering_meshes={key:sources[key]for key in ('mucosa','aponeurosis')},
        covering_muscle_shared_interfaces=[dict(cover='aponeurosis',muscle='superior-longitudinal',
            cover_faces=apo_faces,muscle_faces=sl_faces)],reviewPending=True,packingReady=False,
        numericalCoveringGatePassed=covering_audit['passed'],rootHGSurfaceSectionSelectionRequired=True)
    atomic_json(output/'coupled-inputs.pending.json',pending)
    report=dict(schema='tongue-selected-coupling-audit/v1',complete=True,expected_pairs=44,measured_pairs=44,
        covering_gate_passed=covering_audit['passed'],covering_crossing_pairs=sum(p['summary']['crossing_face_pairs']for p in cover_pairs),
        covering_coplanar_pairs=sum(p['summary']['coplanar_overlap_face_pairs']for p in cover_pairs),
        only_allowed_covering_contact='literal122-face APO/SL reciprocal chart',
        retained_measurements_passed=True,retained_crossing_pairs=retained_audit['crossing_face_pairs'],
        retained_tissue_ownership_claim=False,pending_manifest=str(output/'coupled-inputs.pending.json'),
        root_review_pending=True,installed=False,geometry_modified=False,anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-began)
    remaining();atomic_json(output/'audit.json',report)
    atomic_json(output/'progress.json',dict(report,measurements=measured),overwrite=True)
    print(json.dumps(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ('coverings','intrinsics','extrinsics','styloglossus','features','ancestry','output'):
        parser.add_argument('--'+arg,type=Path,required=True)
    a=parser.parse_args();main(a.coverings,a.intrinsics,a.extrinsics,a.styloglossus,a.features,a.ancestry,a.output)
