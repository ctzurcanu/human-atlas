"""One three-cell mucosa/aponeurosis/SL bank with two literal interfaces.

The exposed outer mucosa and the passed APO/SL arrangement are immutable.
Only source-identified inner mucosal controls are reconciled to the already
emitted APO upper chart; this records a nanometre encoding correction, not
distance welding. No GG, SG, skeleton or current catalogue is modified.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.constraints import check_constraint
from geometry_tools.model import Mesh
from geometry_tools.tissue import pack
from tongue_whole_input import fixed_bone_sources
from tongue_basal_contract import FRAME,complete_pair,require

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
ORIGIN=np.array([0.,1.625,0.])
IDS={'mucosa':'TA98TONGUE:male:mucosa','aponeurosis':'TA98TONGUE:male:aponeurosis',
     'superior-longitudinal':'TA98TONGUE:male:intrinsic:superior-longitudinal'}


def read(path):
    return json.loads(Path(path).read_text())


def effective_words(mesh):
    signs=mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),dtype='i1'))[mesh.faces]
    return (mesh.attributes['normals'][mesh.faces].astype(np.int64)*signs[:,:,None]).astype('<i2')


def triangle_word_multiset(triangles,words):
    result=[]
    for tri,field in zip(triangles,words):
        phases=[tuple(tuple(xyz)+tuple(map(int,n))for xyz,n in zip(np.roll(tri,-r,axis=0),np.roll(field,-r,axis=0)))for r in range(3)]
        result.append(min(phases))
    return sorted(result)


def triangle_word_face_map(source_triangles,source_words,actual_triangles,actual_words):
    def key(triangle,field):
        return min(tuple(tuple(xyz)+tuple(map(int,n))for xyz,n in zip(np.roll(triangle,-r,axis=0),
                   np.roll(field,-r,axis=0)))for r in range(3))
    actual={key(triangle,field):i for i,(triangle,field)in enumerate(zip(actual_triangles,actual_words))}
    require(len(actual)==len(actual_triangles),'Duplicate actual face/word chart is not allowed')
    mapped=[actual.get(key(triangle,field))for triangle,field in zip(source_triangles,source_words)]
    require(None not in mapped and len(set(mapped))==len(source_triangles)==len(actual_triangles),
            'Physical source-face/corner-word bijection is incomplete')
    return mapped


def main(exposed,aponeurosis,source_coverings,output):
    began=time.monotonic()
    def remaining():
        if time.monotonic()-began>50:
            raise TimeoutError('Shared tongue coverings exceeded the 50 second worker budget')
    exposed,aponeurosis,source_coverings,output=(p.resolve()for p in (exposed,aponeurosis,source_coverings,output))
    require(all(p.is_relative_to(WORK)and p.is_file()for p in (exposed,aponeurosis,source_coverings)),
            'Use explicit existing local construction inputs')
    require(output.is_relative_to(WORK)and output!=WORK and not output.exists(),'Use a new local stage')
    e,a,old=read(exposed),read(aponeurosis),read(source_coverings)
    role_json=aponeurosis.with_name('surface-role-map.json')
    role=read(role_json)
    require(a.get('passed')is True and role['interface_result_sha256']==digest_file(aponeurosis),
            'Selected actual APO/SL interface has not passed or changed')
    require(a['whole_input']==old['whole_input']and a['covering_construction_sha256']==digest_file(source_coverings),
            'APO/SL and covering source do not share the reviewed source')
    source_path=Path(e['source']);source=Mesh.load(source_path)
    require(digest_file(source_path)==e['sourceSha256']==old['source_sha256'],
            'Exposed and aponeurotic sources changed')
    m_manifest_path=Path(e['packedManifest']);m_manifest=read(m_manifest_path)
    m_path=m_manifest_path.with_name('TA98TONGUE_male_mucosa.npz');muco=Mesh.load(m_path)
    require(m_manifest.get('passed')is True and digest_file(m_manifest['bank'])==m_manifest['bank_sha256'],
            'Actual exposed finite cell bank changed')
    arrangement_path=Path(a['arrangement'])
    require(digest_file(arrangement_path)==a['arrangement_sha256']==role['arrangement_sha256'],
            'Actual APO/SL arrangement changed')
    role_path=Path(role['role_map'])
    require(digest_file(role_path)==role['role_map_sha256'],'Surface/source identity map changed')
    with np.load(arrangement_path,allow_pickle=False)as z:
        av=z['positions'].copy()
        surfaces={key:z[value].copy()for key,value in a['surfaces'].items()}
        words={key:z[value].copy()for key,value in a['normal_words'].items()}
    with np.load(role_path,allow_pickle=False)as z:
        ordinals=z['APO_mucosa_free_face_ordinals'].copy()
        source_parents=z['APO_free_source_parent_faces'].copy()[ordinals]
    require(len(ordinals)==len(set(ordinals.tolist()))==62 and set(source_parents.tolist())==set(range(62)),
            'APO top needs its exact sixty-two original source-face identities')
    canonical_tri=av[surfaces['vascular-free:aponeurosis'][ordinals]]
    canonical_words=words['vascular-free:aponeurosis'][ordinals]
    provenance_path=exposed.with_name('source-face-provenance.json')
    provenance=read(provenance_path)
    literal={r['sourceFace']:r for r in provenance['faces']}
    outer_count=e['carrierTriangles']
    require(np.array_equal(muco.triangles[:outer_count],Mesh.load(e['carrier']).triangles)
            or np.linalg.norm(muco.triangles[:outer_count]-Mesh.load(e['carrier']).triangles,axis=2).max()<2e-8,
            'Actual exposed outer carrier order changed')
    mt=muco.triangles.copy();mw=effective_words(muco)
    normals=source.attributes['normals'].astype(float);normals/=np.linalg.norm(normals,axis=1)[:,None]
    thickness=.00035+.00007*np.clip(normals[:,1],0,1)-.00005*np.clip(-normals[:,1],0,1)
    inner=source.vertices-thickness[:,None]*normals
    actual_inner_1625=(inner-ORIGIN).astype('<f4').astype(float)+ORIGIN
    actual_inner_1500=(inner-np.array([0.,1.5,0.])).astype('<f4').astype(float)+np.array([0.,1.5,0.])
    replacements={};shared_muco_faces=[];corner_mappings=[]
    for i,parent in enumerate(source_parents):
        source_face=old['aponeurosis_source_faces'][int(parent)]
        require(source_face in literal,'The aponeurosis top must remain a full protected exposed source facet')
        outer_face=literal[source_face]['carrierFace']
        inner_face=outer_count+outer_face
        source_ids=source.faces[source_face]
        expected_old=actual_inner_1500[source_ids]
        selected_order=next((order for order in ([0,1,2],[1,2,0],[2,0,1])
            if np.array_equal(canonical_tri[i],expected_old[order])),None)
        require(selected_order is not None,'APO top differs from its literal original source-control identities')
        source_ids=source_ids[selected_order]
        current=actual_inner_1625[source_ids]
        require({tuple(x)for x in mt[inner_face]}=={tuple(x)for x in current},
                'Mucosal inner face is not the same three actual source controls')
        for sid,old_point,new_point in zip(source_ids,current,canonical_tri[i]):
            key=tuple(old_point)
            if key in replacements:
                require(np.array_equal(replacements[key]['point'],new_point),'A shared source control has inconsistent owners')
            replacements[key]=dict(point=new_point.copy(),source_vertex=int(sid))
        shared_muco_faces.append(inner_face)
        corner_mappings.append(dict(sourceFace=int(source_face),aponeurosisFreeFace=int(ordinals[i]),
                                    mucosaOriginalFace=int(inner_face),canonicalChartFace=i))
    old_tri=mt.copy()
    displacement=[]
    for key,row in replacements.items():
        d=float(np.linalg.norm(np.asarray(key)-row['point']))
        require(d<=1e-8,'Source-owned encoding correction exceeds ten nanometres')
        displacement.append(dict(sourceVertex=row['source_vertex'],before=list(key),after=row['point'].tolist(),distanceM=d))
        where=np.all(mt==np.asarray(key),axis=2)
        require(not np.any(where[:outer_count]),'Inner correction may not touch the exposed outer patch')
        mt[where]=row['point']
    shared_set=set(shared_muco_faces)
    muco_free=np.array([i for i in range(outer_count,len(mt))if i not in shared_set])
    apo_rest=np.array([i for i in range(len(surfaces['vascular-free:aponeurosis']))if i not in set(ordinals.tolist())])
    surface_triangles={'mucosa:outer':mt[:outer_count],'mucosa:free':mt[muco_free],
        'mucosa:aponeurosis':canonical_tri,
        'aponeurosis:free':av[surfaces['vascular-free:aponeurosis'][apo_rest]],
        'aponeurosis:SL':av[surfaces['vascular-contact:superior-longitudinal:aponeurosis']],
        'SL:free':av[surfaces['SL:free']]}
    fields={'mucosa:outer':mw[:outer_count],'mucosa:free':mw[muco_free],
        'mucosa:aponeurosis':canonical_words,
        'aponeurosis:free':words['vascular-free:aponeurosis'][apo_rest],
        'aponeurosis:SL':words['vascular-contact:superior-longitudinal:aponeurosis'],
        'SL:free':words['SL:free']}
    positions=[];lookup={};packed_surfaces={}
    for key,triangles in surface_triangles.items():
        faces=[]
        for triangle in triangles:
            face=[]
            for xyz in triangle:
                literal_key=tuple(xyz)
                if literal_key not in lookup:
                    lookup[literal_key]=len(positions);positions.append(xyz)
                face.append(lookup[literal_key])
            faces.append(face)
        packed_surfaces[key]=np.asarray(faces,int)
    positions=np.asarray(positions)
    require(np.array_equal((positions-ORIGIN).astype('<f4').astype(float)+ORIGIN,positions),
            'The chosen dyadic origin changed an already emitted source position')
    cells={'mucosa':[dict(surface='mucosa:outer',side=1,neighbor=None),
                    dict(surface='mucosa:free',side=1,neighbor=None),
                    dict(surface='mucosa:aponeurosis',side=-1,neighbor='aponeurosis')],
        'aponeurosis':[dict(surface='aponeurosis:free',side=1,neighbor=None),
                      dict(surface='mucosa:aponeurosis',side=1,neighbor='mucosa'),
                      dict(surface='aponeurosis:SL',side=-1,neighbor='superior-longitudinal')],
        'superior-longitudinal':[dict(surface='SL:free',side=1,neighbor=None),
                                dict(surface='aponeurosis:SL',side=1,neighbor='aponeurosis')]}
    remaining()
    output.mkdir(parents=True)
    result=pack(positions,packed_surfaces,cells,IDS,FRAME,output/'packed',origin=ORIGIN,
                corner_words=fields,max_pairs=250000,max_tests=2000000)
    require(result['passed'],'Actual emitted cells failed closure/self/occupied-volume gates')
    decoded={key:Mesh.load(output/'packed'/(pid.replace(':','_')+'.npz'))for key,pid in IDS.items()}
    require(np.array_equal(decoded['mucosa'].triangles[:outer_count],old_tri[:outer_count])
            and np.array_equal(effective_words(decoded['mucosa'])[:outer_count],mw[:outer_count]),
            'Actual emitted exposed outer coordinates or normal words changed')
    literal_proofs={}
    for key in ('aponeurosis','superior-longitudinal'):
        source_mesh=next(r for r in role['meshes']if r['key']==key)
        path=Path(source_mesh['mesh']);require(digest_file(path)==source_mesh['sha256'],'APO/SL saved source changed')
        original=Mesh.load(path)
        with np.load(path,allow_pickle=False)as z:
            original_words=z['corner_normals'].copy()
        require(original.surface_sha==decoded[key].surface_sha
                and triangle_word_multiset(original.triangles,original_words)
                ==triangle_word_multiset(decoded[key].triangles,effective_words(decoded[key])),
                'The actual shared bank changed a passed APO/SL facet or raw corner word')
        face_map=triangle_word_face_map(original.triangles,original_words,decoded[key].triangles,
                                       effective_words(decoded[key]))
        literal_proofs[key]=dict(sourceMesh=str(path),sourceSha256=digest_file(path),
            sourceGeometrySha256=original.geometry_sha,actualGeometrySha256=decoded[key].geometry_sha,
            sourceOrientedSurfaceSha256=original.surface_sha,actualOrientedSurfaceSha256=decoded[key].surface_sha,
            sourceFaceToActualFace=face_map,faceOrderChanged=face_map!=list(range(len(face_map))),
            orientedFacetAndRawCornerWordMultisetLiteral=True)
    interface_proofs=[]
    for first,second in [('mucosa','aponeurosis'),('aponeurosis','superior-longitudinal')]:
        def owned(cell,peer):
            return [i for g in result['parts'][cell]['tissueDrawGroups']if g.get('neighbor')==IDS[peer]
                    for i in range(g['start']//3,(g['start']+g['count'])//3)]
        parameters=dict(first_faces=owned(first,second),second_faces=owned(second,first))
        check=check_constraint(dict(kind='shared-interface',parameters=parameters),decoded[first],decoded[second])
        require(check['passed']and check['literal_opposed_normal_words'],'An emitted reciprocal interface lost literal coordinates/words')
        interface_proofs.append(dict(first=IDS[first],second=IDS[second],parameters=parameters,actual_check=check))
    bone_proofs=[]
    for key,mesh in decoded.items():
        for row in fixed_bone_sources():
            remaining()
            bone=Mesh.load(row['meshSavePath'])
            arrays,summary=pair_contacts(mesh,bone,0.,max_pairs=250000,max_tests=2000000)
            complete_pair(summary)
            witness=output/(key+'--'+row['partId'].replace(':','_').replace('/','_')+'.npz')
            atomic_npz(witness,**arrays)
            bone_proofs.append(dict(first=IDS[key],second=row['partId'],boneGeometrySha256=bone.geometry_sha,
                summary=summary,witnesses=str(witness),witnessesSha256=digest_file(witness)))
    changes=[]
    for row in displacement:
        sid=row['sourceVertex'];outer=(source.vertices[sid]-ORIGIN).astype('<f4').astype(float)+ORIGIN
        before=float(np.linalg.norm(outer-row['before']));after=float(np.linalg.norm(outer-row['after']))
        changes.append(dict(sourceVertex=sid,beforeM=before,afterM=after,changeM=after-before))
    require(max(abs(r['changeM'])for r in changes)<=1e-8,'Actual wall-thickness norm changed beyond the declared encoding correction')
    report=dict(schema='tongue-shared-exposed-coverings/v1',passed=True,installed=False,anatomical_acceptance=False,
        exposedSource=str(exposed),exposedSourceSha256=digest_file(exposed),
        aponeurosisSource=str(aponeurosis),aponeurosisSourceSha256=digest_file(aponeurosis),
        originalCoverings=str(source_coverings),originalCoveringsSha256=digest_file(source_coverings),
        whole_input=old['whole_input'],packed_manifest=str(output/'packed/tissue.json'),
        actual_emitted_checks=result['checks'],source_geometry_modified=False,
        sourceMeshesModified=False,derivedMucosaInnerGeometryChanged=True,outerMucosaCoordinatesAndWordsLiteral=True,
        APO_SLSourceSurfacesAndRawWordsLiteral=literal_proofs,canonicalSourceFaceMapping=corner_mappings,
        sourceOwnedInnerEncodingCorrections=displacement,wallThicknessNormChanges=changes,
        maximumInnerControlDisplacementM=max(r['distanceM']for r in displacement),
        maximumWallThicknessNormChangeM=max(abs(r['changeM'])for r in changes),
        actualPositionReencodingErrorM=result['maximum_position_encoding_error_m'],storageOrigin=ORIGIN.tolist(),
        interfaces=interface_proofs,fixedBoneChecks=bone_proofs,threeCellsOneSharedBank=True,
        actual_meshes=3,additionalAnatomicalMeshes=0,SG_GGReadOrModified=False,
        sourceNormalFieldsIndependentlyRegenerated=False,referenceRenderingStillRequired=True,
        remaining=['Exposed-muscle fits and complete tongue assembly still need source-bound cross checks and actual rendered review.',
                   'The floor-of-mouth reflection and regional anatomical boundaries remain approximate.'],
        elapsed_seconds=time.monotonic()-began)
    remaining()
    atomic_json(output/'construction.json',report)
    print(json.dumps(dict(output=str(output),passed=True,cells=3,bytes=result['total_bytes'],
        maximumInnerControlDisplacementM=report['maximumInnerControlDisplacementM'],
        sharedFaces=[len(p['parameters']['first_faces'])for p in interface_proofs],
        bonePairs=len(bone_proofs),elapsedSeconds=report['elapsed_seconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exposed',type=Path,required=True)
    parser.add_argument('--aponeurosis',type=Path,required=True)
    parser.add_argument('--source-coverings',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    main(args.exposed,args.aponeurosis,args.source_coverings,args.output)
