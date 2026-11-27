"""One source-bound SL/aponeurosis complement with a reciprocal finite chart.

Only the portion of the approximate superior-longitudinal extent penetrating
the actual finite aponeurosis is removed. The aponeurosis is not thinned or
deformed. Its boundary is subdivided to reuse one literal contact chart with
opposed raw corner-normal words. This script does not pack or install a bank.
"""
from pathlib import Path
import argparse
import itertools
import json
import sys
import time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.tissue import validate_cells
from geometry_tools.algorithms.finite_tissue_boolean import subtract_occupied,make_solid,same_oriented_triangles
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature
from geometry_tools.algorithms.float32_diagonal_repair import repair_diagonals
from tongue_whole_input import WORK,whole_input,fixed_bone_sources,new_output
from tongue_interface_triangulation import repair as repair_common_storage_triangulation


def corner_words(path,mesh):
    with np.load(path,allow_pickle=False)as z:
        if 'corner_normals'in z:
            words=z['corner_normals'].copy()
        else:
            words=mesh.attributes['normals'][mesh.faces].astype(np.int32)
            if 'normalSigns'in mesh.attributes:
                words*=mesh.attributes['normalSigns'][mesh.faces][:,:,None]
            words=words.astype('<i2')
    if words.shape!=(len(mesh.faces),3,3):raise ValueError('Actual source corner-normal field is missing')
    return words


def emit(path,positions,faces,words,frame,source):
    used,local=np.unique(faces,return_inverse=True);local=local.reshape(-1,3)
    mesh=Mesh(positions[used],local,frame,source=source)
    atomic_npz(path,vertices=mesh.vertices,faces=mesh.faces,frame=np.array(frame),
        source_json=np.array(json.dumps(source,sort_keys=True)),corner_normals=words)
    return Mesh.load(path)


def normal_proof(mesh,words):
    t=mesh.triangles;g=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);length=np.linalg.norm(g,axis=1)
    if np.any(length==0):raise ValueError('Emitted zero-area face')
    g/=length[:,None];n=words.astype(float);nl=np.linalg.norm(n,axis=2)
    if np.any(nl==0):raise ValueError('Emitted zero normal word')
    n/=nl[:,:,None];dots=np.einsum('fki,fi->fk',n,g)
    return dict(complete_face_positive_plane_dot=bool(np.all(dots>0)),minimum_corner_plane_dot=float(dots.min()),
        complete_face_proof='Positive barycentric combinations of corners retain positive face-plane dot and cannot be zero')


def diagnose_storage(output,stage,checkpoint):
    """One saved-array test of the precise Float32 triangulation producer fix."""
    began=time.monotonic();output=new_output(output);output.mkdir(parents=True)
    record=json.loads((stage/'interface-result.json').read_text())
    with np.load(checkpoint/'boolean-checkpoint.npz',allow_pickle=False)as raw:
        v=raw['positions'].copy()+np.array([0.,1.5,0.]);free=raw['faces_2'].copy();parents=raw['parents_2'].copy()
        shared=raw['faces_1'].copy();sl_faces=np.r_[raw['faces_0'],shared]
    with np.load(stage/'reciprocal-arrangement.npz',allow_pickle=False)as saved:
        q=saved['positions'].copy();cn=saved['normal_words_2'].copy();contact_words=saved['normal_words_1'].copy()
    tri=q[free];g=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);g/=np.linalg.norm(g,axis=1)[:,None]
    bad=np.flatnonzero(np.any(np.einsum('fki,fi->fk',cn.astype(float),g)<=0,axis=1))
    witnesses=[]
    for i in bad:
        raw_tri=v[free[i]];rn=np.cross(raw_tri[1]-raw_tri[0],raw_tri[2]-raw_tri[0]);rn/=np.linalg.norm(rn)
        widths=np.linalg.norm(tri[i]-np.roll(tri[i],-1,axis=0),axis=1)
        witnesses.append(dict(free_face=int(i),source_parent_face=int(parents[i]),
            original_to_encoded_geometric_normal_angle_degrees=float(np.degrees(np.arccos(np.clip(rn@g[i],-1,1)))),
            minimum_encoded_edge_m=float(widths.min())))
    try:
        repaired,proof=repair_diagonals(v,free,np.array([0.,1.5,0.]),stored_positions=q,
            corner_normals=cn,plane_tolerance=2e-12)
        fields=cn.copy()
        for change in proof['changes']:
            ids=change['faceIds'];old=fields[ids].copy()
            fields[ids]=np.array([[old[s[0],s[1]]for s in row]for row in change['newCornerSources']])
            change['source_parent_faces']=[int(parents[i])for i in ids]
        if boundary_signature(repaired)!=boundary_signature(free):raise ValueError('Free-chart boundary changed')
        a=Mesh(q,sl_faces,'model:local-ta98:assembled-metres')
        b=Mesh(q,np.r_[repaired,shared[:,[0,2,1]]],a.frame)
        bwords=np.r_[fields,-contact_words[:,[0,2,1]]]
        _,contacts=pair_contacts(a,b,0.,max_pairs=100000,max_tests=1000000)
        validation=validate_cells({'SL':a,'aponeurosis':b},max_pairs=250000,max_tests=2000000)
        npf=normal_proof(b,bwords)
        passed=bool(validation['passed']and npf['complete_face_positive_plane_dot']and
            contacts['complete_within_pair_budget']and contacts['crossing_face_pairs']==0 and
            contacts['coplanar_overlap_face_pairs']==len(shared))
        atomic_npz(output/'diagnostic-arrangement.npz',positions=q,SL_faces=sl_faces,apo_free_faces=repaired,
            shared_faces=shared,apo_free_corner_words=fields,shared_corner_words=contact_words)
        report=dict(passed=passed,diagnostic_only=True,geometry_coordinates_unchanged=True,
            reciprocal_chart_unchanged=True,free_chart_boundary_unchanged=True,
            normal_witnesses=witnesses,diagonal_repair=proof,contacts=contacts,cell_checks=validation,
            normal_proof=npf,anatomical_parameters_changed=False,anatomical_acceptance=False,
            source_interface_result_sha256=digest_file(stage/'interface-result.json'),
            source_checkpoint_sha256=digest_file(checkpoint/'boolean-checkpoint.npz'),
            producer_correction='Apply exact-boundary-preserving coplanar interior-diagonal repair to the aponeurosis free chart after Float32 embedding, using actual inherited corner normals as the inversion guard; reuse shared interface chart unchanged',
            elapsed_seconds=time.monotonic()-began)
    except ValueError as error:
        report=dict(passed=False,diagnostic_only=True,normal_witnesses=witnesses,error=str(error),
            source_interface_result_sha256=digest_file(stage/'interface-result.json'),
            source_checkpoint_sha256=digest_file(checkpoint/'boolean-checkpoint.npz'),
            no_thresholds_relaxed=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'diagnosis.json',report)
    print((output/'diagnosis.json').read_text()[:13000],flush=True)


def main(output,base,coverings,checkpoint_only=False,resume_checkpoint=None,storage_aware_retriangulation=False):
    began=time.monotonic();output=new_output(output)
    construction=base/'construction.json';record=json.loads(construction.read_text())
    if not record['numerical_and_resource_passed']or not record['fixed_bone_crossings_zero']:
        raise ValueError('Use the passed corrected-whole intrinsic construction')
    binding=record['whole_input'];whole,whole_path,actual_binding=whole_input(binding['whole_mesh'],binding['fit_result_path'])
    if actual_binding!=binding:raise ValueError('Whole input evidence changed')
    covering_report=json.loads((coverings/'construction.json').read_text())
    if not covering_report['passed']or covering_report['whole_input']!=binding:
        raise ValueError('Coverings must derive from the exact same corrected whole')
    tissue=json.loads((coverings/'packed/tissue.json').read_text())
    if not tissue['passed']or digest_file(tissue['bank'])!=tissue['bank_sha256']:
        raise ValueError('Actual covering bank is stale or failed')
    rows={r['key']:r for r in record['meshes']}
    for row in rows.values():
        if digest_file(row['mesh'])!=row['sha256']:raise ValueError('Intrinsic source changed')
    sl_row=rows['superior-longitudinal'];sl=Mesh.load(sl_row['mesh'])
    apo_path=coverings/'packed/TA98TONGUE_male_aponeurosis.npz';apo=Mesh.load(apo_path)
    if apo.geometry_sha!=tissue['checks']['cells']['aponeurosis']['geometry_sha256']:
        raise ValueError('Decoded actual aponeurosis cell changed')
    muco_path=coverings/'packed/TA98TONGUE_male_mucosa.npz';muco=Mesh.load(muco_path)
    output.mkdir(parents=True)
    origin=np.array([0.,1.5,0.]);sp=sl.vertices-origin;ap=apo.vertices-origin
    uses={'superior-longitudinal':[dict(surface='SL:free',side=1,neighbor=None)]}
    old_surfaces={'SL:free':sl.faces}
    if resume_checkpoint is None:
        result=subtract_occupied(sp,old_surfaces,{'SL:free':dict(id='SL:free',kind='outer')},
            {'superior-longitudinal':sl.faces},
            {'superior-longitudinal':[('SL:free',i,1)for i in range(len(sl.faces))]},uses,
            ap,apo.faces,occupied_face_owners=np.full(len(apo.faces),'aponeurosis',object),
            partition_free_outer=True,exact_source_replay=True,reconcile_shared_subdivisions=True,
            diagnostic_directory=output/'boolean-replay')
    else:
        prior=json.loads((resume_checkpoint/'boolean-checkpoint.json').read_text())
        if prior['base_intrinsic_construction_sha256']!=digest_file(construction)or prior['aponeurosis_mesh_sha256']!=digest_file(apo_path):
            raise ValueError('Resume checkpoint belongs to different literal inputs')
        with np.load(resume_checkpoint/'boolean-checkpoint.npz',allow_pickle=False)as z:
            saved_surfaces={r['surface']:z[f'faces_{r["index"]}'].copy()for r in prior['charts']}
            saved_charts={r['surface']:dict(parentFaces=z[f'parents_{r["index"]}'].copy(),
                barycentrics=z[f'barycentrics_{r["index"]}'].copy())for r in prior['charts']}
            contact_name=next(s for s in saved_surfaces if s.startswith('vascular-contact:'))
            free_name=next(s for s in saved_surfaces if s.startswith('vascular-free:'))
            result=dict(positions=z['positions'].copy(),surfaces=saved_surfaces,charts=saved_charts,evidence=prior['evidence'],
                cells={'superior-longitudinal':[dict(surface='SL:free',side=1,neighbor=None),
                    dict(surface=contact_name,side=1,neighbor='aponeurosis')]},
                vascularContacts=[dict(surface=contact_name)],vascularFreeSurfaces=[dict(surface=free_name)])
    surfaces=result['surfaces'];contact=[u['surface']for u in result['vascularContacts']]
    free=[u['surface']for u in result['vascularFreeSurfaces']]
    if len(contact)!=1 or len(free)!=1:raise ValueError('Expected one reciprocal contact chart and one free aponeurosis chart')
    contact=contact[0];free=free[0]
    checkpoint_arrays={'positions':result['positions']}
    checkpoint_rows=[]
    for i,(sid,faces)in enumerate(surfaces.items()):
        checkpoint_arrays[f'faces_{i}']=faces
        checkpoint_arrays[f'parents_{i}']=result['charts'][sid]['parentFaces']
        checkpoint_arrays[f'barycentrics_{i}']=result['charts'][sid]['barycentrics']
        bary=result['charts'][sid]['barycentrics']
        bad=np.argwhere((bary< -2e-10)|(bary>1+2e-10))
        checkpoint_rows.append(dict(surface=sid,index=i,minimum_barycentric=float(bary.min()),
            maximum_barycentric=float(bary.max()),out_of_guard_entries=len(bad),
            witnesses=[dict(entry=ids.tolist(),weights=bary[ids[0],ids[1]].tolist(),
                parent_face=int(result['charts'][sid]['parentFaces'][ids[0]]),
                point=result['positions'][faces[ids[0],ids[1]]].tolist())for ids in bad[:6]]))
    atomic_npz(output/'boolean-checkpoint.npz',**checkpoint_arrays)
    atomic_json(output/'boolean-checkpoint.json',dict(geometry_parameters_unchanged=True,
        base_intrinsic_construction_sha256=digest_file(construction),
        aponeurosis_mesh_sha256=digest_file(apo_path),charts=checkpoint_rows,evidence=result['evidence'],
        checkpoint_only=checkpoint_only,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began))
    if checkpoint_only:
        print(json.dumps(dict(checkpoint_saved=True,geometry_parameters_unchanged=True,charts=checkpoint_rows,
            elapsed_seconds=time.monotonic()-began)),flush=True)
        return
    cells={**result['cells'],'aponeurosis':[dict(surface=free,side=1,neighbor=None),
        dict(surface=contact,side=-1,neighbor='superior-longitudinal')]}
    local=result['positions'];positions=local.astype('<f4').astype(float)+origin
    encoding=float(np.linalg.norm(positions-(local+origin),axis=1).max())
    # Literal source prefixes and the three untouched intrinsic NPZs remain
    # unchanged. Only newly computed cut controls are rounded for actual storage.
    if not np.array_equal(positions[:len(sl.vertices)],sl.vertices):raise ValueError('An old SL control moved')
    if not np.array_equal(positions[len(sl.vertices):len(sl.vertices)+len(apo.vertices)],apo.vertices):
        raise ValueError('An old aponeurosis control moved')
    original_words={'SL:free':corner_words(sl_row['mesh'],sl),'aponeurosis':corner_words(apo_path,apo)}
    words={};source_normal_proof=[]
    for sid,faces in surfaces.items():
        p=positions[faces];g=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);g/=np.linalg.norm(g,axis=1)[:,None]
        if sid==contact:
            out=np.repeat(np.rint(g*32767).astype('<i2')[:,None,:],3,axis=1)
            family='literal reciprocal geometric normal chart'
        else:
            source_mesh=sl if sid=='SL:free'else apo;source_key='SL:free'if sid=='SL:free'else'aponeurosis'
            parent=result['charts'][sid]['parentFaces'];bary=result['charts'][sid]['barycentrics'].copy()
            # A literal corner has exact one-hot coordinates. Do this before
            # any least-squares weight guard; picometre source edges can make
            # an otherwise exact original corner ill-conditioned numerically.
            literal_weight_replays=0
            source_local=source_mesh.vertices-origin
            for i,face in enumerate(faces):
                t=source_local[source_mesh.faces[parent[i]]]
                for k,index in enumerate(face):
                    same=np.flatnonzero(np.all(t==local[index],axis=1))
                    if len(same):bary[i,k]=np.eye(3)[int(same[0])];literal_weight_replays+=1
            if float(bary.min())< -2e-10 or float(bary.max())>1+2e-10:raise ValueError('Normal chart left an actual parent face')
            interpolated=np.einsum('fki,fij->fkj',bary,original_words[source_key][parent].astype(float))
            interpolated/=np.linalg.norm(interpolated,axis=2)[:,:,None]
            out=np.rint(interpolated*32767).astype('<i2')
            literal=0
            for i,face in enumerate(faces):
                source_face=source_mesh.faces[parent[i]]
                for k,index in enumerate(face):
                    same=np.flatnonzero(np.all(source_mesh.vertices[source_face]==positions[index],axis=1))
                    if len(same):out[i,k]=original_words[source_key][parent[i],int(same[0])];literal+=1
            family='Actual source-face/barycentric artist chart, literal original words retained at source corners'
            source_normal_proof.append(dict(surface=sid,literal_source_corners=literal,parent_family=source_key,
                exact_one_hot_parent_corner_weights=literal_weight_replays))
        words[sid]=out
    storage_repair=None
    if storage_aware_retriangulation:
        repaired,repaired_words,storage_repair=repair_common_storage_triangulation(local+origin,positions,
            surfaces[free],result['charts'][free]['parentFaces'],words[free],surfaces[contact],
            np.r_[surfaces['SL:free'],surfaces[contact]],sl.frame,parent_faces=(75,78,80),time_limit=45.)
        if not storage_repair['passed']:
            witnesses=[];sl_parents=result['charts']['SL:free']['parentFaces']
            for row in storage_repair.get('candidate_results',[]):
                for _,face in row.get('SL_defect_face_pairs',[]):
                    witnesses.append(dict(apo_parent=row['parent'],candidate=row['ordinal'],SL_face=int(face),
                        SL_original_parent_face=int(sl_parents[face])if face<len(sl_parents)else None,
                        hits_immutable_reciprocal_chart=face>=len(sl_parents)))
            storage_repair['coupled_SL_witnesses']=witnesses
        atomic_json(output/'storage-repair-proof.json',storage_repair)
        if not storage_repair['passed']:
            atomic_json(output/'interface-result.json',dict(passed=False,anatomical_acceptance=False,installed=False,
                finite_constraint_failure=True,storage_repair_proof_sha256=digest_file(output/'storage-repair-proof.json'),
                whole_input=binding,base_intrinsic_construction_sha256=digest_file(construction),
                actual_aponeurosis_sha256=digest_file(apo_path),geometry_coordinates_changed=False,
                reciprocal_chart_changed=False,elapsed_seconds=time.monotonic()-began))
            print(json.dumps(dict(passed=False,finite_constraint_failure=True,
                failed_parent_faces=storage_repair.get('failed_parent_faces'),
                elapsed_seconds=time.monotonic()-began)),flush=True)
            return
        surfaces[free]=repaired;words[free]=repaired_words
    atomic_npz(output/'reciprocal-arrangement.npz',positions=positions,storage_origin=origin,
        **{f'surface_{i}':f for i,f in enumerate(surfaces.values())},
        **{f'normal_words_{i}':words[s]for i,s in enumerate(surfaces)})
    shared_meshes={};emitted_rows=[];normals=[]
    for name,cell_uses in cells.items():
        faces=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in cell_uses])
        nw=np.concatenate([words[u['surface']]if u['side']==1 else -words[u['surface']][:,[0,2,1]]for u in cell_uses])
        path=output/(name+'.npz')
        mesh=emit(path,positions,faces,nw,sl.frame,dict(kind='source-bound-reciprocal-aponeurosis-SL-interface',
            whole_input=binding,parent_intrinsic_construction_sha256=digest_file(construction),
            actual_aponeurosis_mesh_sha256=digest_file(apo_path),aponeurosis_thickness_changed=False,
            representation='approximate-macroscopic-muscle-extent-reference-envelope'if name=='superior-longitudinal'else'finite-aponeurosis-cell',
            anatomical_acceptance=False))
        shared_meshes[name]=mesh;proof=normal_proof(mesh,nw);normals.append(dict(cell=name,**proof))
        emitted_rows.append(dict(key=name,mesh=str(path),sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,
            vertices=len(mesh.vertices),triangles=len(mesh.faces),corner_normal_source=dict(path=str(path),sha256=digest_file(path),array='corner_normals')))
    # The old aponeurosis/mucosa interface consists of the first62 decoded
    # aponeurosis facets. Any unexpected subdivision there would require a
    # coupled mucosal reconstruction and is rejected in this focused stage.
    outer_count=next(u['count']//3 for u in tissue['parts']['aponeurosis']['tissueDrawGroups']if 'neighbor'in u)
    parents=result['charts'][free]['parentFaces'];outer=surfaces[free][parents<outer_count]
    def coordinate_faces(v,f):
        points=v[f];keys=[]
        for t in points:
            rows=[tuple(p)for p in t];keys.append(min(tuple(rows),tuple(rows[1:]+rows[:1]),tuple(rows[2:]+rows[:2])))
        return sorted(keys)
    outer_literal=coordinate_faces(positions,outer)==coordinate_faces(apo.vertices,apo.faces[:outer_count])
    if not outer_literal:raise ValueError('SL complement subdivided the mucosa/aponeurosis outer interface; coupled propagation required')
    validation=validate_cells(shared_meshes,max_pairs=250000,max_tests=2000000)
    source_solid=make_solid(sp,sl.faces)[0];apo_solid=make_solid(ap,apo.faces)[0]
    removed=float((source_solid^apo_solid).volume());before=float(source_solid.volume())
    after=validation['cell_volumes_m3'].get('superior-longitudinal',0.)
    actual_apo=validation['cell_volumes_m3'].get('aponeurosis',0.)
    apo_before=float(apo_solid.volume());extent=float(np.ptp(np.r_[sp,ap],axis=0).max())
    volume_bound=2048*np.finfo(float).eps*extent**3*(len(sl.faces)+len(apo.faces))
    preservation=dict(before_SL_volume_m3=before,removed_overlap_volume_m3=removed,after_SL_volume_m3=after,
        SL_retained_volume_fraction=after/before,SL_complement_volume_residual_m3=after-(before-removed),
        apo_before_volume_m3=apo_before,apo_after_volume_m3=actual_apo,apo_volume_difference_m3=actual_apo-apo_before,
        actual_volume_roundoff_bound_m3=volume_bound,original_source_controls_unchanged=True,
        original_mucosa_aponeurosis_chart_literal=outer_literal,aponeurosis_thickness_parameters_modified=False)
    all_intrinsics={key:(shared_meshes[key]if key=='superior-longitudinal'else Mesh.load(row['mesh']))for key,row in rows.items()}
    pairs=[]
    for(a,ma),(b,mb)in itertools.combinations(all_intrinsics.items(),2):
        _,summary=pair_contacts(ma,mb,0.,max_pairs=100000,max_tests=1000000)
        pairs.append(dict(first=a,second=b,summary=summary))
    bone_pairs=[]
    bone_sources=[(row,Mesh.load(row['meshSavePath']))for row in fixed_bone_sources()]
    for name,mesh in {**all_intrinsics,'aponeurosis':shared_meshes['aponeurosis']}.items():
        for row,bone in bone_sources:
            _,summary=pair_contacts(mesh,bone,0.,max_pairs=100000,max_tests=1000000)
            bone_pairs.append(dict(first=name,second=row['partId'],bone_geometry_sha256=row['geometrySha256'],summary=summary))
    sl_mesh=shared_meshes['superior-longitudinal'];q=np.r_[sl_mesh.vertices,sl_mesh.triangles.mean(1)]
    minimum_winding=float(np.abs(igl.winding_number(whole.vertices,whole.faces,q)).min())
    _,muco_contact=pair_contacts(sl_mesh,muco,0.,max_pairs=100000,max_tests=1000000)
    _,shared_contact=pair_contacts(sl_mesh,shared_meshes['aponeurosis'],0.,max_pairs=100000,max_tests=1000000)
    clear=lambda r:r['complete_within_pair_budget']and r['crossing_face_pairs']==0 and r['coplanar_overlap_face_pairs']==0
    passed=bool(validation['passed']and all(clear(r['summary'])for r in pairs+bone_pairs)and minimum_winding>.99 and clear(muco_contact)
        and shared_contact['complete_within_pair_budget']and shared_contact['crossing_face_pairs']==0
        and all(r['complete_face_positive_plane_dot']for r in normals)
        and abs(actual_apo-apo_before)<=volume_bound and abs(after-(before-removed))<=volume_bound and after/before>.8)
    manifest=dict(schema='tongue-aponeurosis-SL-reciprocal-interface/v1',whole_input=binding,
        base_intrinsic_construction=str(construction),base_intrinsic_construction_sha256=digest_file(construction),
        covering_construction_sha256=digest_file(coverings/'construction.json'),aponeurosis_source_sha256=digest_file(apo_path),
        source_aponeurosis_geometry_sha256=apo.geometry_sha,meshes=emitted_rows,
        arrangement=str(output/'reciprocal-arrangement.npz'),arrangement_sha256=digest_file(output/'reciprocal-arrangement.npz'),
        surfaces={s:f'surface_{i}'for i,s in enumerate(surfaces)},normal_words={s:f'normal_words_{i}'for i,s in enumerate(surfaces)},
        cells=cells,reciprocal_interface=dict(surface=contact,triangles=len(surfaces[contact]),
            positions_stored_once=True,exact_opposite_raw_normal_words=True,
            first='superior-longitudinal',first_side=1,second='aponeurosis',second_side=-1),
        boolean_evidence=result['evidence'],source_normal_proof=source_normal_proof,normal_checks=normals,
        storage_aware_triangulation=storage_repair,
        actual_position_encoding_error_m=encoding,occupied_volume_preservation=preservation,actual_emitted_cell_checks=validation,
        intrinsic_pairs=pairs,fixed_bone_pairs=bone_pairs,source_minimum_absolute_winding=minimum_winding,
        mucosa_contact=muco_contact,shared_interface_contact=shared_contact,
        unchanged_intrinsic_NPZs=[dict(key=k,path=r['mesh'],sha256=r['sha256'])for k,r in rows.items()if k!='superior-longitudinal'],
        resumed_boolean_checkpoint=None if resume_checkpoint is None else dict(path=str(resume_checkpoint),
            sha256=digest_file(resume_checkpoint/'boolean-checkpoint.npz'),geometry_recomputed=False),
        passed=passed,anatomical_acceptance=False,installed=False,bank_packed=False,
        construction_script_sha256=digest_file(__file__),elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'interface-result.json',manifest)
    new_rows=[]
    sl_emitted=next(r for r in emitted_rows if r['key']=='superior-longitudinal')
    for old in record['meshes']:
        new_rows.append({**old,**sl_emitted}if old['key']=='superior-longitudinal'else old)
    resource={**record['resource'],'mesh_records':len(new_rows),
        'total_triangles':sum(r['triangles']for r in new_rows),
        'estimated_packed_raw_bytes':sum(r['vertices']*18+r['triangles']*12 for r in new_rows)}
    resource['passed']=(resource['mesh_records']<=resource['maximum_mesh_records']and
        resource['total_triangles']<=resource['maximum_triangles']and
        resource['estimated_packed_raw_bytes']<=resource['maximum_raw_bytes'])
    atomic_json(output/'construction.json',{**record,'meshes':new_rows,
        'parent_construction_sha256':digest_file(construction),'interface_result_sha256':digest_file(output/'interface-result.json'),
        'interface_aponeurosis':next(r for r in emitted_rows if r['key']=='aponeurosis'),
        'numerical_and_resource_passed':passed and resource['passed'],'resource':resource,
        'visual_review':'pending','installed':False,
        'remaining_limits':['Macro muscle extents remain approximate; retained extrinsic intralingual interfaces are not exclusive tissue segmentation']})
    show=[dict(mesh=r['mesh'],name=r['key'],color=rows[r['key']]['color'],corner_normal_source=r.get('corner_normal_source'))for r in new_rows]
    apo_row=next(r for r in emitted_rows if r['key']=='aponeurosis')
    show.append(dict(mesh=apo_row['mesh'],name='Actual aponeurosis, reciprocal shared boundary',color=[.75,.8,.85,.7],corner_normal_source=apo_row['corner_normal_source']))
    atomic_json(output/'review-spec.json',dict(frame=whole.frame,resolution=640,wall_seconds=60,meshes=show,
        views=[dict(name='right',direction=[1,0,0]),dict(name='top',direction=[0,1,0]),
            dict(name='sagittal-interface',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.625,.055],clear_outer=True))]))
    print(json.dumps(dict(passed=passed,SL_triangles=len(sl_mesh.faces),aponeurosis_triangles=len(shared_meshes['aponeurosis'].faces),
        reciprocal_faces=len(surfaces[contact]),SL_retained_fraction=after/before,
        bone_pairs=len(bone_pairs),bone_crossings=sum(r['summary']['crossing_face_pairs']for r in bone_pairs),
        apo_volume_change_m3=actual_apo-apo_before,elapsed_seconds=manifest['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--base',type=Path,default=WORK/'basal-macro-envelopes-v1-corner-charts')
    ap.add_argument('--coverings',type=Path,default=WORK/'coverings-corrected-whole-v1')
    ap.add_argument('--checkpoint-only',action='store_true')
    ap.add_argument('--resume-checkpoint',type=Path)
    ap.add_argument('--diagnose-saved-stage',type=Path)
    ap.add_argument('--storage-aware-retriangulation',action='store_true')
    args=ap.parse_args()
    if args.diagnose_saved_stage:
        if args.resume_checkpoint is None:raise ValueError('A saved Boolean checkpoint is required for diagnosis')
        diagnose_storage(args.output,args.diagnose_saved_stage.resolve(),args.resume_checkpoint.resolve())
    else:
        main(args.output,args.base.resolve(),args.coverings.resolve(),args.checkpoint_only,
            args.resume_checkpoint.resolve()if args.resume_checkpoint else None,args.storage_aware_retriangulation)
