"""Finite C++ Minkowski erosion of the literal local source tongue."""
from pathlib import Path
import argparse
import json
import sys
import time
import importlib.metadata
import numpy as np
import igl
import manifold3d as md

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import topology
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from tongue_whole_input import whole_input,new_output


def main(out,whole=None,fit_result=None,storage_origin=None,validate_inputs_only=False):
    began=time.monotonic();out=new_output(out)
    native,source_path,binding=whole_input(whole,fit_result)
    source=Path(binding['source_catalogue_path'])
    if validate_inputs_only:
        print(json.dumps(dict(validated=True,generation_started=False,whole_input=binding)),flush=True)
        return
    out.mkdir(parents=True)
    origin=np.array([0.,1.625,.055]);radius=.0009
    storage_origin=origin if storage_origin is None else np.asarray(storage_origin,float)
    solid,_=make_solid(native.vertices-origin,native.faces)
    kernel=md.Manifold.sphere(radius,12);km=kernel.to_mesh64()
    kv=np.asarray(km.vert_properties)[:,:3];kf=np.asarray(km.tri_verts,dtype=np.int64);tri=kv[kf]
    normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    ri=float(np.min(np.abs(np.einsum('ij,ij->i',normal,tri[:,0]))/np.linalg.norm(normal,axis=1)))
    ro=float(np.linalg.norm(kv,axis=1).max());error=ro-ri
    assert error<=.00015
    atomic_npz(out/'convex-kernel.npz',vertices=kv,faces=kf)
    evidence=dict(schema='tongue-source-minkowski-erosion/v1',source_catalogue_sha256=digest_file(source),
        source_mesh_sha256=digest_file(source_path),source_geometry_sha256=native.geometry_sha,
        manifold3d_version=importlib.metadata.version('manifold3d'),requested_erosion_m=radius,
        method='Manifold64 morphological erosion via minkowski_difference with one convex geodesic sphere',
        official_installed_api_documentation=md.Manifold.minkowski_difference.__doc__,
        kernel=dict(vertices=len(kv),triangles=len(kf),inscribed_radius_m=ri,circumscribed_radius_m=ro,
            maximum_radial_error_m=error,error_requirement_m=.00015,center=[0.,0.,0.]),
        source_geometry_modified=False,whole_input=binding,anatomical_acceptance=False,installed=False)
    atomic_json(out/'erosion-inputs.json',evidence)
    inner=solid.minkowski_difference(kernel)
    assert inner.status()==md.Error.NoError and not inner.is_empty()
    result=inner.to_mesh64();v=np.asarray(result.vert_properties)[:,:3]+origin
    f=np.asarray(result.tri_verts,dtype=np.int64)
    emitted=(v-storage_origin).astype(np.float32).astype(float)+storage_origin
    encoding_error=float(np.linalg.norm(emitted-v,axis=1).max())
    mesh=Mesh(emitted,f,native.frame,{'normals':igl.per_vertex_normals(emitted,f)},dict(
        kind='actual-source-minkowski-inner-envelope',source_catalogue_sha256=digest_file(source),
        source_geometry_sha256=native.geometry_sha,erosion_input_sha256=digest_file(out/'erosion-inputs.json'),
        kernel_inscribed_radius_m=ri,kernel_circumscribed_radius_m=ro,complete_tissue_segmentation=False,
        whole_input=binding))
    path=out/'native-inner-envelope.npz';mesh.save(path)
    q=np.concatenate([mesh.vertices,mesh.triangles.mean(1),(mesh.triangles[:,0]+mesh.triangles[:,1])/2,
        (mesh.triangles[:,1]+mesh.triangles[:,2])/2,(mesh.triangles[:,2]+mesh.triangles[:,0])/2])
    winding=np.abs(igl.winding_number(native.vertices,native.faces,q))
    distances=np.sqrt(igl.point_mesh_squared_distance(q,native.vertices,native.faces)[0])
    _,contact=pair_contacts(native,mesh,0.,max_pairs=100000,max_tests=1000000)
    topo=topology(mesh);cpp=diagnostics(mesh)
    checks=dict(closed_oriented=topo['closed_oriented_manifold'],one_component=len(topo['components'])==1,
        self_crossings_zero=cpp['self_intersection_pairs']==0,contained=bool(np.all(winding>.99)),
        outer_crossings_zero=contact['crossing_face_pairs']==0 and contact['coplanar_overlap_face_pairs']==0,
        contact_budget_complete=contact['complete_within_pair_budget'],
        sampled_clearance_at_least_kernel_inscribed_radius=float(distances.min())>=ri-encoding_error-1e-12)
    report={**evidence,'mesh':str(path),'mesh_sha256':digest_file(path),'geometry_sha256':mesh.geometry_sha,
        'topology':topo,'cpp_diagnostics':cpp,'outer_contacts':contact,'independent_query_count':len(q),
        'minimum_absolute_winding':float(winding.min()),'minimum_sampled_clearance_m':float(distances.min()),
        'position_encoding_error_m':encoding_error,'before_volume_m3':float(solid.volume()),
        'storage_origin':storage_origin.tolist(),
        'after_volume_m3':float(inner.volume()),'retained_volume_fraction':float(inner.volume()/solid.volume()),
        'checks':checks,'passed':all(checks.values()),'elapsed_seconds':time.monotonic()-began}
    atomic_json(out/'erosion-check.json',report)
    print(json.dumps(dict(passed=report['passed'],triangles=len(f),components=len(topo['components']),
        kernel_radial_error_m=error,retained_volume_fraction=report['retained_volume_fraction'],
        minimum_sampled_clearance_m=report['minimum_sampled_clearance_m'],checks=checks,
        elapsed_seconds=report['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--whole-tongue',type=Path);ap.add_argument('--fit-result',type=Path)
    ap.add_argument('--storage-origin',nargs=3,type=float)
    ap.add_argument('--validate-inputs-only',action='store_true');args=ap.parse_args()
    main(args.output,args.whole_tongue,args.fit_result,args.storage_origin,args.validate_inputs_only)
