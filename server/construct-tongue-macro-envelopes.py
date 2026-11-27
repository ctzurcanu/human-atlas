"""One original, explicitly nonexclusive macroscopic muscle-extent proposal.

Four render records: continuous SL extent, two full IL bellies, one shared T/V
carrier. No native external-muscle envelope is treated as an exclusive cutter.
"""
from pathlib import Path
import argparse
import importlib.util
import itertools
import json
import sys
import time
import numpy as np
import igl
import manifold3d as md

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature
from geometry_tools.algorithms.float32_diagonal_repair import repair_diagonals
from geometry_tools.tissue import layer_arrangement
from tongue_whole_input import WORK,whole_input,checked_inner_envelope,fixed_bone_sources,new_output


def regenerated_il(native):
    """Re-query the selected whole; never transport a previously emitted cell."""
    path=ROOT/'server/construct-tongue-intrinsic-ribbons.py'
    spec=importlib.util.spec_from_file_location('tongue_il_whole_regeneration',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    env=module.UpperBodyEnvelope(native);result={}
    for side,sign in [('left',1.),('right',-1.)]:
        v,f=module.longitudinal(env,sign*.62,sign,False)
        result[side]=Mesh(v,f,native.frame,{'normals':igl.per_vertex_normals(v,f)},
            dict(kind='inferior-longitudinal-macro-extent-regenerated-from-selected-whole',
                 whole_tongue_geometry_sha256=native.geometry_sha,
                 method='Original longitudinal sweep re-evaluated against actual corrected whole ray intervals',
                 prior_cells_deformed=False,approximate_reference_extent=True))
    return result


def main(out,storage_origin=None,whole=None,fit_result=None,inner_evidence=None,validate_inputs_only=False):
    began=time.monotonic();out=new_output(out);work=WORK
    native,whole_path,binding=whole_input(whole,fit_result)
    source=Path(binding['source_catalogue_path']);source_sha=binding['source_catalogue_sha256']
    origin=np.array([0.,1.625,.055]);tolerance=.00001;maximum_measured_error=.00002;separation=.00002;depth=.0035
    storage_origin=origin if storage_origin is None else np.asarray(storage_origin,float)
    def solid(m):return make_solid(m.vertices-origin,m.faces)[0]
    if inner_evidence is None:
        if whole is not None:
            raise ValueError('Corrected-whole construction requires explicit --inner-envelope-evidence from its new erosion stage')
        inner_evidence=work/'intrinsic-territories-v2/erosion-check.json'
    inner_mesh,inner_evidence,inner_check=checked_inner_envelope(inner_evidence,native,binding)
    mask=json.loads((work/'intrinsic-territories-v2/dorsum-mask-regularization.json').read_text())
    controls=work/'surface-controls-v1/surface-controls.json';assert digest_file(controls)==mask['source_controls_sha256']
    ids=np.array(sorted(set(mask['original_source_faces'])|set(mask['added_source_faces'])),int)
    bone_rows=fixed_bone_sources()
    if validate_inputs_only:
        print(json.dumps(dict(validated=True,generation_started=False,whole_input=binding,
            inner_envelope_evidence=str(inner_evidence),inner_envelope_evidence_sha256=digest_file(inner_evidence),
            fixed_bone_ids=[r['partId']for r in bone_rows],
            corrected_intrinsic_policy='Regenerate SL, both IL and shared T/V from selected whole; no prior cells deformed')),
            flush=True)
        return
    out.mkdir(parents=True)
    native_solid=solid(native)
    inner=solid(inner_mesh)
    used,lookup=np.unique(native.faces[ids],return_inverse=True);faces=lookup.reshape(-1,3)
    top=native.vertices[used];bottom=top-np.array([0.,depth,0.])
    positions,surfaces,cells,proof=layer_arrangement(np.array([top,bottom]),faces,['dorsal-depth-classifier'])
    closed=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]
        for u in cells['dorsal-depth-classifier']])
    classifier=make_solid(positions-origin,closed)[0]
    kernels=[];erosion={}
    for d in [.0007,.0025]:
        kernel=md.Manifold.sphere(d,16);km=kernel.to_mesh64();v=np.asarray(km.vert_properties)[:,:3];f=np.asarray(km.tri_verts);t=v[f]
        n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);rmin=float(np.min(abs(np.einsum('ij,ij->i',n,t[:,0]))/np.linalg.norm(n,axis=1)))
        rmax=float(np.linalg.norm(v,axis=1).max());assert rmax-rmin<=.00015
        kernels.append(dict(requested_radius_m=d,inscribed_radius_m=rmin,circumscribed_radius_m=rmax,radial_error_m=rmax-rmin))
        erosion[d]=native_solid.minkowski_difference(kernel)
    sl_raw=(erosion[.0007]-erosion[.0025])^classifier
    assert sl_raw.status()==md.Error.NoError and not sl_raw.is_empty()
    sl=sl_raw.simplify(tolerance)
    il=(regenerated_il(native)if whole is not None else
        {side:Mesh.load(work/('intrinsics-v5/inferior-longitudinal-'+side+'.npz'))for side in ['left','right']})
    # The tiny construction separation is visual/numerical, not a histological interstitium.
    cut=md.Manifold.batch_boolean([sl]+[solid(m)for m in il.values()],md.OpType.Add)
    cut=cut.minkowski_sum(md.Manifold.cube([2*separation]*3,center=True))
    tv_raw=inner-cut;assert tv_raw.status()==md.Error.NoError and not tv_raw.is_empty()
    tv=tv_raw.simplify(tolerance)
    rows=[];meshes=[];checks=[];simplification=[]
    definitions=[('superior-longitudinal',sl,sl_raw,['A05.1.04.106']),
        ('inferior-longitudinal-left',il['left'],None,['A05.1.04.107']),
        ('inferior-longitudinal-right',il['right'],None,['A05.1.04.107']),
        ('transverse-vertical-shared',tv,tv_raw,['A05.1.04.108','A05.1.04.109'])]
    colors={'superior-longitudinal':[.72,.3,.29,1.],'inferior-longitudinal-left':[.68,.34,.31,1.],
        'inferior-longitudinal-right':[.68,.34,.31,1.],'transverse-vertical-shared':[.78,.43,.36,1.]}
    for key,value,raw,terms in definitions:
        if isinstance(value,Mesh):
            vertices=value.vertices;f=value.faces;n=value.attributes['normals'];encoding=0.
        else:
            assert value.status()==md.Error.NoError
            vm=value.to_mesh64();vertices=np.asarray(vm.vert_properties)[:,:3]+origin;f=np.asarray(vm.tri_verts,dtype=np.int64)
            rawm=raw.to_mesh64();rv=np.asarray(rawm.vert_properties)[:,:3]+origin;rf=np.asarray(rawm.tri_verts,dtype=np.int64)
            def samples(v,f):
                t=v[f];return np.r_[v,t.mean(1),(t[:,0]+t[:,1])/2,(t[:,1]+t[:,2])/2,(t[:,2]+t[:,0])/2]
            qa=samples(rv,rf);qb=samples(vertices,f)
            measured=max(float(np.sqrt(igl.point_mesh_squared_distance(qa,vertices,f)[0]).max()),
                         float(np.sqrt(igl.point_mesh_squared_distance(qb,rv,rf)[0]).max()))
            simplification_row=dict(key=key,raw_triangles=len(rf),simplified_triangles=len(f),
                documented_surface_movement_bound_m=tolerance,bidirectional_sampled_error_m=measured,
                actual_error_acceptance_bound_m=maximum_measured_error,
                documented_tolerance_met_by_independent_samples=measured<=tolerance+1e-10,
                official_installed_api_documentation=md.Manifold.simplify.__doc__,passed=measured<=maximum_measured_error+1e-10)
            simplification.append(simplification_row)
            atomic_json(out/(key+'-simplification.json'),simplification_row)
            if measured>maximum_measured_error+1e-10:
                raw_path=out/(key+'-raw-float64-diagnostic.npz')
                Mesh(rv,rf,native.frame,source=dict(kind='rejected-before-encoding-boolean')).save(raw_path)
                Mesh(vertices,f,native.frame,source=dict(kind='rejected-simplified-float64')).save(out/(key+'-simplified-float64-diagnostic.npz'))
                raise ValueError(f'{key}: actual sampled simplification error {measured} m exceeds explicit {maximum_measured_error} m bound')
            n=igl.per_vertex_normals(vertices,f)
        pre_storage=Mesh(vertices,f,native.frame,{'normals':n},dict(kind='identical-macro-roi-before-storage'))
        pre_folder=out/'pre-storage';pre_folder.mkdir(exist_ok=True)
        pre_path=pre_folder/(key+'.npz');pre_storage.save(pre_path)
        stored=(vertices-storage_origin).astype(np.float32).astype(float)+storage_origin
        encoding=float(np.linalg.norm(vertices-stored,axis=1).max())
        repaired,repair=repair_diagonals(vertices,f,storage_origin,stored_positions=stored,
            plane_tolerance=max(2e-12,4*encoding))
        assert boundary_signature(repaired)==boundary_signature(f)
        for change in repair['changes']:
            assert boundary_signature(np.asarray(change['oldFaces']))==boundary_signature(np.asarray(change['newFaces']))
        if repair['diagonalsChanged']:
            n=igl.per_vertex_normals(stored,repaired)
        f=repaired;vertices=stored
        mesh=Mesh(vertices,f,native.frame,{'normals':n},dict(kind='original-macroscopic-muscle-extent-reference-envelope',
            source_catalogue_sha256=source_sha,whole_tongue_geometry_sha256=native.geometry_sha,
            terminology=terms,exclusive_tissue_volume=False,complete_histological_segmentation=False,
            interwoven_fibres_in_envelope=True,retained_external_muscles_used_as_cutters=False,
            whole_input=binding,
            qualification='Original approximate macroscopic extent reference, not patient-specific muscle-fibre segmentation'))
        path=out/(key+'.npz')
        if path.exists():
            saved=Mesh.load(path)
            assert saved.geometry_sha==mesh.geometry_sha and np.array_equal(saved.attributes['normals'],mesh.attributes['normals'])
        else:mesh.save(path)
        meshes.append((key,mesh))
        if key.startswith('inferior-longitudinal'):
            old=il[key.rsplit('-',1)[1]];assert pre_storage.geometry_sha==old.geometry_sha
        topo=topology(mesh);cpp=diagnostics(mesh)
        q=np.r_[mesh.vertices,mesh.triangles.mean(1),(mesh.triangles[:,0]+mesh.triangles[:,1])/2,
                (mesh.triangles[:,1]+mesh.triangles[:,2])/2,(mesh.triangles[:,2]+mesh.triangles[:,0])/2]
        winding=np.abs(igl.winding_number(native.vertices,native.faces,q))
        _,contacts=pair_contacts(native,mesh,0.,max_pairs=100000,max_tests=1000000)
        flags=dict(closed_oriented=topo['closed_oriented_manifold'],one_component=len(topo['components'])==1,
            self_crossings_zero=cpp['self_intersection_pairs']==0,whole_tongue_contained=bool(np.all(winding>.99)),
            source_budget_complete=contacts['complete_within_pair_budget'],
            source_crossings_zero=contacts['crossing_face_pairs']==0 and contacts['coplanar_overlap_face_pairs']==0)
        checks.append(dict(key=key,topology=topo,cpp_diagnostics=cpp,source_contacts=contacts,
            minimum_absolute_winding=float(winding.min()),checks=flags,passed=all(flags.values())))
        rows.append(dict(key=key,mesh=str(path),sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,
            terms=terms,vertices=len(vertices),triangles=len(f),components=len(topo['components']),color=colors[key],
            representation='macroscopic-muscle-extent-reference-envelope',exclusive_tissue_volume=False,
            encoding_error_m=encoding,pre_storage_mesh=str(pre_path),pre_storage_geometry_sha256=pre_storage.geometry_sha,
            diagonal_repair=repair,repair_shape_deviation_bound_m=repair['maximumPlaneResidualMeters'],anatomical_acceptance=False))
        print(json.dumps(dict(stage=key,triangles=len(f),components=len(topo['components']),checks=flags)),flush=True)
    pair_rows=[]
    for (a,ma),(b,mb)in itertools.combinations(meshes,2):
        _,summary=pair_contacts(ma,mb,0.,max_pairs=50000,max_tests=1000000)
        pair_rows.append(dict(first=a,second=b,scope='intrinsic macro-envelope diagnostic',**summary))
    context=json.loads((work/'intrinsics-v5/retained-context-audit/summary.json').read_text())
    for a,ma in meshes:
        for row in context['context']:
            assert digest_file(row['mesh'])==row['sha256'];mb=Mesh.load(row['mesh'])
            _,summary=pair_contacts(ma,mb,0.,max_pairs=50000,max_tests=1000000)
            pair_rows.append(dict(first=a,second=row['key'],scope='retained nonexclusive macro-envelope interdigitation diagnostic',**summary))
    bone_pairs=[]
    for a,ma in meshes:
        for row in bone_rows:
            mb=Mesh.load(row['meshSavePath'])
            _,summary=pair_contacts(ma,mb,0.,max_pairs=100000,max_tests=1000000)
            bone_pairs.append(dict(first=a,second=row['partId'],scope='mandatory actual fixed-bone exclusion',
                bone_mesh=row['meshSavePath'],bone_mesh_sha256=row['meshSaveSha256'],
                bone_geometry_sha256=row['geometrySha256'],**summary))
    total=sum(r['triangles']for r in rows);raw_bytes=sum(r['vertices']*18+r['triangles']*12 for r in rows)
    resource=dict(mesh_records=len(rows),maximum_mesh_records=4,total_triangles=total,maximum_triangles=30000,
        estimated_packed_raw_bytes=raw_bytes,maximum_raw_bytes=1000000,additional_images=0,
        passed=len(rows)==4 and total<=30000 and raw_bytes<=1000000)
    complete=all(r['complete_within_pair_budget']for r in pair_rows)
    intrinsic_pairs_clear=all(r['crossing_face_pairs']==0 and r['coplanar_overlap_face_pairs']==0
        for r in pair_rows if r['scope']=='intrinsic macro-envelope diagnostic')
    bone_complete=all(r['complete_within_pair_budget']for r in bone_pairs)
    bone_clear=bone_complete and all(r['crossing_face_pairs']==0 and r['coplanar_overlap_face_pairs']==0 for r in bone_pairs)
    passed=all(r['passed']for r in checks)and all(r['passed']for r in simplification)and resource['passed']and complete and intrinsic_pairs_clear and bone_clear
    result=dict(schema='tongue-macroscopic-muscle-reference-envelopes/v1',source_catalogue_sha256=source_sha,
        whole_input=binding,selected_whole_mesh=str(whole_path),selected_whole_mesh_sha256=digest_file(whole_path),
        source_tongue_geometry_sha256=binding['parent_whole_geometry_sha256'],
        construction_whole_geometry_sha256=native.geometry_sha,
        intrinsic_regeneration=dict(all_cells_regenerated_from_selected_whole=whole is not None,
            prior_cells_deformed=False,legacy_native_default_preserved=whole is None,
            inferior_longitudinal_method='Original sweep re-evaluated on whole ray intervals'if whole is not None else'Legacy v5 native-source bellies'),
        frame=native.frame,storage_origin=storage_origin.tolist(),source_geometry_modified=False,primary_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
        vendor_geometry_copied=False,source_dorsum_controls_sha256=digest_file(controls),dorsum_mask=mask,
        dorsal_classifier_depth_m=depth,sl_band_requested_depths_m=[.0007,.0025],erosion_kernel_errors=kernels,
        inner_envelope_evidence_sha256=digest_file(inner_evidence),simplification=simplification,
        construction_separation_m=separation,separation_is_histological_claim=False,
        meshes=rows,checks=checks,pair_diagnostics=pair_rows,bone_diagnostics=bone_pairs,
        intrinsic_pairs_crossing_free=intrinsic_pairs_clear,
        fixed_bone_pairs_complete=bone_complete,fixed_bone_crossings_zero=bone_clear,
        resource=resource,numerical_and_resource_passed=passed,
        inter_object_crossings_are_not_hidden_or_asserted_zero=True,exclusive_tissue_partition_claimed=False,
        shared_tv_representation=dict(one_mesh=True,terms=['A05.1.04.108','A05.1.04.109'],
            directional_controls=dict(transverse=[1,0,0],vertical=[0,1,0]),duplicate_occupied_volumes=False,
            orientation_fields_are_illustrative=True),
        source_reference_analysis_sha256=digest_file(work/'reference-shape-v1/assessment.json'),
        construction_script_sha256=digest_file(__file__),anatomical_acceptance=False,visual_review='pending',installed=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'construction.json',result)
    review=dict(frame=native.frame,resolution=640,wall_seconds=60,
        meshes=[dict(mesh=r['mesh'],name=r['key'],color=r['color'])for r in rows],
        views=[dict(name='top',direction=[0,1,0]),dict(name='right',direction=[1,0,0]),dict(name='front',direction=[0,0,1]),
               dict(name='section',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.625,.055],clear_outer=True))])
    atomic_json(out/'review-spec.json',review)
    atomic_json(out/'review-with-context.json',{**review,'meshes':review['meshes']+[
        dict(mesh=str(whole_path),name='Selected whole tongue envelope',color=[.6,.7,.8,.12])]})
    print(json.dumps(dict(numerical_and_resource_passed=passed,resource=resource,elapsed_seconds=result['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--storage-origin',nargs=3,type=float)
    ap.add_argument('--whole-tongue',type=Path);ap.add_argument('--fit-result',type=Path)
    ap.add_argument('--inner-envelope-evidence',type=Path)
    ap.add_argument('--validate-inputs-only',action='store_true');args=ap.parse_args()
    main(args.output,args.storage_origin,args.whole_tongue,args.fit_result,
         args.inner_envelope_evidence,args.validate_inputs_only)
