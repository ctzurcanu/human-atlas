"""One finite Boolean ownership proposal, never an anatomical completion claim."""
from pathlib import Path
import argparse
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
from geometry_tools.algorithms.finite_tissue_boolean import make_solid


def solid_mesh(solid,origin,frame,source):
    value=solid.to_mesh64()
    vertices=np.asarray(value.vert_properties)[:,:3]+origin
    faces=np.asarray(value.tri_verts,dtype=np.int64)
    return Mesh(vertices,faces,frame,{'normals':igl.per_vertex_normals(vertices,faces)},source)


def main(previous,out):
    began=time.monotonic();out.mkdir(parents=True,exist_ok=True)
    record=json.loads((previous/'construction.json').read_text())
    checked=json.loads((previous/'checked-geometry.json').read_text())
    context=json.loads((previous/'retained-context-audit/summary.json').read_text())
    assert checked['passed'] and checked['construction_sha256']==digest_file(previous/'construction.json')
    source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(source)==record['source_catalogue_sha256']
    origin=np.array([0.,1.625,.055]);clearance=.00005
    input_rows=[];cutters=[]
    for row in context['context']:
        assert digest_file(row['mesh'])==row['sha256']
        mesh=Mesh.load(row['mesh']);topo=topology(mesh)
        assert topo['closed_oriented_manifold'] and len(topo['components'])==1
        solid,_=make_solid(mesh.vertices-origin,mesh.faces)
        cutters.append(solid)
        input_rows.append(dict(key=row['key'],path=row['mesh'],sha256=row['sha256'],geometry_sha256=mesh.geometry_sha,
            topology=topo,retained_geometry_modified=False))
    atomic_json(out/'inputs.json',dict(schema='tongue-intrinsic-ownership-inputs/v1',
        previous_construction_sha256=digest_file(previous/'construction.json'),
        previous_context_audit_sha256=digest_file(previous/'retained-context-audit/summary.json'),
        source_catalogue_sha256=digest_file(source),cutters=input_rows,clearance_m=clearance,
        clearance_method='Minkowski union dilation with a centered 0.10 mm cube; contains a 0.05 mm Euclidean-radius ball',
        predeclared_rejection_gates=dict(any_original_bundle_erased=True,any_original_bundle_retains_less_than=.50,
            any_group_retains_less_than=.65,components_above_twice_original=True,total_triangles_above=15000),
        anatomical_acceptance=False,installed=False))
    occupied=md.Manifold.batch_boolean(cutters,md.OpType.Add)
    assert occupied.status()==md.Error.NoError
    dilated=occupied.minkowski_sum(md.Manifold.cube([2*clearance]*3,center=True))
    assert dilated.status()==md.Error.NoError
    print(json.dumps(dict(stage='cutter-ready',triangles=dilated.num_tri(),elapsed_seconds=time.monotonic()-began)),flush=True)
    cutter=solid_mesh(dilated,origin,record['frame'],dict(kind='illustrative-ownership-clearance-cutter',
        actual_retained_source_hashes=input_rows,clearance_m=clearance,anatomical_tissue=False))
    cutter.save(out/'clearance-cutter.npz')
    rows=[];group_audits=[];rejection=[]
    for row in record['meshes']:
        assert digest_file(row['mesh'])==row['sha256']
        host=Mesh.load(row['mesh']);pieces=[];bundle_rows=[]
        before=0.;after=0.
        for ordinal,component in enumerate(host.source['components']):
            start=component['vertex_start'];end=start+component['vertex_count']
            faces=host.faces[np.all((host.faces>=start)&(host.faces<end),axis=1)]-start
            solid,_=make_solid(host.vertices[start:end]-origin,faces)
            result=solid-dilated
            assert result.status()==md.Error.NoError
            volume_before=float(solid.volume());volume_after=float(result.volume())
            before+=volume_before;after+=volume_after
            decomposed=result.decompose();component_volumes=[float(p.volume()) for p in decomposed]
            ratio=volume_after/volume_before
            bundle_rows.append(dict(ordinal=ordinal,before_m3=volume_before,after_m3=volume_after,
                removed_m3=volume_before-volume_after,retained_fraction=ratio,components=len(decomposed),
                component_volumes_m3=component_volumes,erased=result.is_empty()))
            if result.is_empty():
                rejection.append(row['key']+f' bundle {ordinal} erased')
            else:
                pieces.append(solid_mesh(result,origin,record['frame'],{}))
            if ratio<.50:rejection.append(row['key']+f' bundle {ordinal} retains less than 50%')
        vv=[];ff=[];index=0
        for mesh in pieces:
            vv.append(mesh.vertices);ff.append(mesh.faces+index);index+=len(mesh.vertices)
        if not vv:
            rejection.append(row['key']+' entirely erased');continue
        v=np.concatenate(vv);f=np.concatenate(ff)
        mesh=Mesh(v,f,record['frame'],{'normals':igl.per_vertex_normals(v,f)},dict(
            kind='original-fascicle-finite-ownership-proposal',source_fascicle_sha256=row['sha256'],
            source_catalogue_sha256=digest_file(source),boolean_method='Manifold64 host minus literal retained envelope union with explicit cube dilation',
            complete_muscle_volume_segmentation=False,anatomical_acceptance=False,retained_geometry_modified=False,
            limitation='Representative surviving fascicle territories only. Finite cutter boundaries do not prove histological tissue boundaries or attachments.'))
        target=out/(row['key']+'.npz');mesh.save(target)
        topo=topology(mesh)
        if after/before<.65:rejection.append(row['key']+' group retains less than 65%')
        if len(topo['components'])>2*row['bundles']:rejection.append(row['key']+' fragments exceed twice the original bundle count')
        if not topo['closed_oriented_manifold']:rejection.append(row['key']+' emitted topology is not closed oriented manifold')
        group=dict(key=row['key'],before_m3=before,after_m3=after,removed_m3=before-after,
            retained_fraction=after/before,original_bundles=row['bundles'],emitted_components=len(topo['components']),
            bundles=bundle_rows,topology=topo)
        atomic_json(out/(row['key']+'-ownership.json'),group);group_audits.append(group)
        rows.append({**row,'mesh':str(target),'sha256':digest_file(target),'geometry_sha256':mesh.geometry_sha,
            'bundles':len(topo['components']),'original_bundles':row['bundles'],'triangles':len(f),'vertices':len(v),
            'representation':'approximate-fascicle-ownership-fragments','anatomical_acceptance':False})
        print(json.dumps(dict(stage=row['key'],retained_fraction=after/before,components=len(topo['components']),triangles=len(f))),flush=True)
    if sum(r['triangles']for r in rows)>15000:rejection.append('Total emitted triangles exceed 15,000')
    if len(rows)!=7:rejection.append('Not all seven muscle records survive')
    construction={**record,'schema':'tongue-intrinsic-finite-ownership/v1','meshes':rows,
        'total_meshes':len(rows),'total_triangles':sum(r['triangles']for r in rows),
        'total_raw_geometry_bytes':sum(r['vertices']*18+r['triangles']*12 for r in rows),
        'construction_script_sha256':digest_file(__file__),'previous_construction_sha256':digest_file(previous/'construction.json'),
        'clearance_m':clearance,'group_ownership_audits':group_audits,'rejection_reasons':rejection,
        'numerical_preservation_gates_passed':not rejection,'anatomical_acceptance':False,
        'complete_muscle_volume_segmentation':False,'source_geometry_modified':False,'installed':False,
        'state':'rejected-preservation-or-resource-gates' if rejection else 'constructed-awaiting-independent-audit-and-rendered-review',
        'elapsed_seconds':time.monotonic()-began}
    atomic_json(out/'construction.json',construction)
    review=json.loads((previous/'review-spec.json').read_text())
    review['meshes']=[dict(mesh=r['mesh'],name=r['key'],color=r['color'])for r in rows]
    atomic_json(out/'review-spec.json',review)
    context_meshes=json.loads((previous/'review-with-context.json').read_text())['meshes'][7:]
    atomic_json(out/'review-with-context.json',{**review,'meshes':review['meshes']+context_meshes})
    print(json.dumps(dict(stage='constructed',state=construction['state'],total_triangles=construction['total_triangles'],
        rejection_reasons=rejection,elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--previous',required=True,type=Path);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();out=args.output.resolve()
    if not out.is_relative_to(ROOT/'.local-models/geometry-workbench/tongue-completion'):
        raise ValueError('Keep local ownership proposals in the tongue workbench')
    main(args.previous.resolve(),out)
