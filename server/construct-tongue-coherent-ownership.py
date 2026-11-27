"""Source-bound continuous SL sheet and one shared T/V territory carrier.

This is approximate anatomical territory construction. Intrinsic directional
fields share a carrier; they are not invented exclusive histological slabs.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import itertools
import numpy as np
import igl
import manifold3d as md

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import topology
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import make_solid


def main(out):
    began=time.monotonic();out.mkdir(parents=True,exist_ok=True);work=out.parent
    source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    native=Mesh.load(work/'source-v1/meshes/native-00.npz');source_sha=digest_file(source)
    assert source_sha==native.source['catalogue_sha256']
    eroded=json.loads((work/'intrinsic-territories-v2/erosion-check.json').read_text())
    assert eroded['passed'] and eroded['source_catalogue_sha256']==source_sha
    assert digest_file(eroded['mesh'])==eroded['mesh_sha256']
    mask=json.loads((work/'intrinsic-territories-v2/dorsum-mask-regularization.json').read_text())
    control_path=work/'surface-controls-v1/surface-controls.json'
    assert digest_file(control_path)==mask['source_controls_sha256']
    source_faces=sorted(set(mask['original_source_faces'])|set(mask['added_source_faces']))
    context=json.loads((work/'intrinsics-v5/retained-context-audit/summary.json').read_text())
    origin=np.array([0.,1.625,.055]);clearance=.00005
    def source_solid(mesh):return make_solid(mesh.vertices-origin,mesh.faces)[0]
    native_solid=source_solid(native)
    retained=[]
    for row in context['context']:
        assert digest_file(row['mesh'])==row['sha256']
        retained.append((row['key'],Mesh.load(row['mesh'])))
    occupied=md.Manifold.batch_boolean([source_solid(m)for _,m in retained],md.OpType.Add)
    cube=md.Manifold.cube([2*clearance]*3,center=True)
    occupied_gap=occupied.minkowski_sum(cube)
    kernel_rows=[];erosions={}
    for depth in [.0007,.0025]:
        kernel=md.Manifold.sphere(depth,16);value=kernel.to_mesh64()
        v=np.asarray(value.vert_properties)[:,:3];f=np.asarray(value.tri_verts);t=v[f]
        n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        rmin=float(np.min(abs(np.einsum('ij,ij->i',n,t[:,0]))/np.linalg.norm(n,axis=1)))
        rmax=float(np.linalg.norm(v,axis=1).max())
        assert rmax-rmin<=.00015
        erosions[depth]=native_solid.minkowski_difference(kernel)
        assert erosions[depth].status()==md.Error.NoError and not erosions[depth].is_empty()
        kernel_rows.append(dict(requested_depth_m=depth,kernel_triangles=len(f),inscribed_radius_m=rmin,
            circumscribed_radius_m=rmax,maximum_radial_error_m=rmax-rmin))
    # Literal selected dorsal source triangles define finite downward prisms.
    # Their lower caps lie below the source; only the eroded SL depth band is used.
    columns=[];bottom=native.bounds[0,1]-.005
    prism_faces=np.array([[0,1,2],[5,4,3],[0,3,1],[1,3,4],[1,4,2],[2,4,5],[2,5,0],[0,5,3]])
    for face_id in source_faces:
        top=native.triangles[face_id];lower=top.copy();lower[:,1]=bottom
        v=np.r_[top,lower]-origin
        columns.append(make_solid(v,prism_faces)[0])
    classifier=md.Manifold.batch_boolean(columns,md.OpType.Add)
    sl_band=erosions[.0007]-erosions[.0025]
    sl_raw=sl_band^classifier
    sl=sl_raw-occupied_gap
    assert sl.status()==md.Error.NoError and not sl.is_empty()
    sl_gap=sl.minkowski_sum(cube)
    groups={'superior-longitudinal':sl};before={'superior-longitudinal':float(sl_raw.volume())}
    original_components={'superior-longitudinal':1}
    for side in ['left','right']:
        reference=Mesh.load(work/('intrinsics-v5/inferior-longitudinal-'+side+'.npz'))
        solid=source_solid(reference);groups['inferior-longitudinal-'+side]=(solid-occupied_gap)-sl_gap
        before['inferior-longitudinal-'+side]=float(solid.volume());original_components['inferior-longitudinal-'+side]=1
    protected=md.Manifold.batch_boolean(list(groups.values()),md.OpType.Add).minkowski_sum(cube)
    inner=source_solid(Mesh.load(eroded['mesh']))
    tv=(inner-occupied_gap)-protected
    assert tv.status()==md.Error.NoError and not tv.is_empty()
    groups['transverse-vertical-shared']=tv
    before['transverse-vertical-shared']=float(inner.volume());original_components['transverse-vertical-shared']=1
    rows=[];meshes=[];check_rows=[]
    colors={'superior-longitudinal':[.72,.3,.29,1.],
        'inferior-longitudinal-left':[.68,.34,.31,1.],'inferior-longitudinal-right':[.68,.34,.31,1.],
        'transverse-vertical-shared':[.78,.43,.36,1.]}
    for key,solid in groups.items():
        assert solid.status()==md.Error.NoError and not solid.is_empty()
        value=solid.to_mesh64();v=np.asarray(value.vert_properties)[:,:3]+origin;f=np.asarray(value.tri_verts,dtype=np.int64)
        emitted=(v-origin).astype(np.float32).astype(float)+origin
        mesh=Mesh(emitted,f,native.frame,{'normals':igl.per_vertex_normals(emitted,f)},dict(
            kind='original-coherent-tongue-territory',source_geometry_sha256=native.geometry_sha,
            source_catalogue_sha256=source_sha,territory=key,
            source_dorsum_faces=source_faces if key=='superior-longitudinal'else [],
            representation='one shared interwoven T/V carrier'if key=='transverse-vertical-shared'else 'finite approximate muscle territory',
            complete_histological_partition=False,anatomical_acceptance=False,
            ownership_clearance_m=clearance,all_retained_source_geometry_unchanged=True))
        path=out/(key+'.npz');mesh.save(path);meshes.append((key,mesh))
        topo=topology(mesh);cpp=diagnostics(mesh)
        q=np.concatenate([mesh.vertices,mesh.triangles.mean(1),(mesh.triangles[:,0]+mesh.triangles[:,1])/2,
            (mesh.triangles[:,1]+mesh.triangles[:,2])/2,(mesh.triangles[:,2]+mesh.triangles[:,0])/2])
        winding=np.abs(igl.winding_number(native.vertices,native.faces,q))
        _,contact=pair_contacts(native,mesh,0.,max_pairs=100000,max_tests=1000000)
        checks=dict(closed_oriented=topo['closed_oriented_manifold'],self_crossings_zero=cpp['self_intersection_pairs']==0,
            source_containment=bool(np.all(winding>.99)),source_contacts_complete=contact['complete_within_pair_budget'],
            source_crossings_zero=contact['crossing_face_pairs']==0 and contact['coplanar_overlap_face_pairs']==0)
        check_rows.append(dict(key=key,topology=topo,cpp_diagnostics=cpp,source_contacts=contact,
            minimum_absolute_winding=float(winding.min()),checks=checks,passed=all(checks.values())))
        terms=(['A05.1.04.106']if key=='superior-longitudinal'else ['A05.1.04.108','A05.1.04.109']
            if key=='transverse-vertical-shared'else ['A05.1.04.107'])
        rows.append(dict(key=key,terms=terms,mesh=str(path),sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,
            triangles=len(f),vertices=len(emitted),components=len(topo['components']),color=colors[key],
            before_m3=before[key],after_m3=float(solid.volume()),retained_fraction=float(solid.volume())/before[key],
            representation='shared-interwoven-orientation-carrier'if key=='transverse-vertical-shared'else 'approximate-coherent-territory'))
        print(json.dumps(dict(stage=key,triangles=len(f),components=len(topo['components']),checks=checks)),flush=True)
    pair_rows=[]
    for (name,a),(other,b)in itertools.combinations(meshes,2):
        _,summary=pair_contacts(a,b,0.,max_pairs=50000,max_tests=1000000)
        pair_rows.append(dict(first=name,second=other,kind='new-territory-pair',**summary))
    for name,a in meshes:
        for other,b in retained:
            _,summary=pair_contacts(a,b,0.,max_pairs=50000,max_tests=1000000)
            winding=np.abs(igl.winding_number(b.vertices,b.faces,np.r_[a.vertices,a.triangles.mean(1)]))
            pair_rows.append(dict(first=name,second=other,kind='retained-context',**summary,
                candidate_points_inside_retained=int(np.count_nonzero(winding>.99))))
    pair_pass=all(p['complete_within_pair_budget']and p['crossing_face_pairs']==0 and p['coplanar_overlap_face_pairs']==0
        and p.get('candidate_points_inside_retained',0)==0 for p in pair_rows)
    numerical_pass=all(r['passed']for r in check_rows)and pair_pass
    shape_limits=[]
    if rows[0]['components']!=1:shape_limits.append('SL sheet is no longer one connected territory after occupied-solid exclusion')
    for r in rows[1:3]:
        if r['components']!=1:shape_limits.append(r['key']+' remains fragmented by literal external-muscle envelopes')
    report=dict(schema='tongue-coherent-territories/v1',source_catalogue_sha256=source_sha,
        source_geometry_sha256=native.geometry_sha,source_dorsum_controls_sha256=digest_file(control_path),
        source_mask_regularization=mask,erosion_check_sha256=digest_file(work/'intrinsic-territories-v2/erosion-check.json'),
        erosion_kernels=kernel_rows,meshes=rows,total_meshes=len(rows),total_triangles=sum(r['triangles']for r in rows),
        total_raw_geometry_bytes=sum(r['vertices']*18+r['triangles']*12 for r in rows),
        shared_tv_fields=dict(transverse=dict(direction=[1,0,0]),vertical=dict(direction=[0,1,0]),
            coexistence='Both directions describe one interwoven carrier, never duplicate occupied volumes',
            limitation='Source-model coordinate directions are illustrative, not diffusion-measured fibre directions'),
        checks=check_rows,pair_checks=pair_rows,numerical_passed=numerical_pass,shape_limits=shape_limits,
        literal_retained_muscles_unchanged=True,anatomical_acceptance=False,installed=False,
        complete_muscle_volume_segmentation=False,construction_script_sha256=digest_file(__file__),
        elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'construction.json',report)
    review=dict(frame=native.frame,resolution=640,wall_seconds=60,
        meshes=[dict(mesh=r['mesh'],name=r['key'],color=r['color'])for r in rows],
        views=[dict(name='top',direction=[0,1,0]),dict(name='right',direction=[1,0,0]),dict(name='front',direction=[0,0,1]),
            dict(name='section',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.625,.055],clear_outer=True))])
    atomic_json(out/'review-spec.json',review)
    atomic_json(out/'review-with-context.json',{**review,'meshes':review['meshes']+[
        dict(mesh=str(work/'source-v1/meshes/native-00.npz'),name='Fixed native tongue envelope',color=[.6,.7,.8,.12])]})
    print(json.dumps(dict(numerical_passed=numerical_pass,total_meshes=len(rows),total_triangles=report['total_triangles'],
        shape_limits=shape_limits,elapsed_seconds=report['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve()
    if not out.is_relative_to(ROOT/'.local-models/geometry-workbench/tongue-completion'):
        raise ValueError('Keep original local proposals in the tongue workbench')
    main(out)
