"""One finite SL-clearance subtraction from the original curved TV core.

Only TV changes. The exclusion is the actual unchanged SL Minkowski-summed
with a polyhedron containing the planned 0.20 mm Euclidean ball. Its finite
radial conservatism is measured, not hidden as an anatomical thickness change.
This worker saves a candidate; independent original assembly gates still apply.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'server')]
import numpy as np
import manifold3d as md
from scipy.spatial import ConvexHull
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH

spec=importlib.util.spec_from_file_location('tongue_corridor',ROOT/'server/construct-tongue-tv-corridor.py')
corridor=importlib.util.module_from_spec(spec)
spec.loader.exec_module(corridor)
WORK,FRAME,ORIGIN=corridor.WORK,corridor.FRAME,corridor.ORIGIN


def require(ok,message):
    if not ok:raise ValueError(message)


def solid(mesh):
    vertices,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True)
    faces=inverse[mesh.faces]
    result=md.Manifold(md.Mesh64(np.ascontiguousarray(vertices-ORIGIN,dtype=np.float64),
                                np.ascontiguousarray(faces,dtype=np.uint64)))
    require(result.status()==md.Error.NoError,'Native input manifold rejected actual source geometry')
    return result


def export(native):
    require(native.status()==md.Error.NoError,'Native Boolean construction failed')
    raw=native.to_mesh64()
    return Mesh(np.asarray(raw.vert_properties)[:,:3]+ORIGIN,
                np.asarray(raw.tri_verts,dtype=np.int64),FRAME)


def kernel():
    phi=(1+np.sqrt(5.))/2
    vertices=np.array([(0,a,b*phi) for a in (-1,1) for b in (-1,1)]
        +[(a,b*phi,0) for a in (-1,1) for b in (-1,1)]
        +[(a*phi,0,b) for a in (-1,1) for b in (-1,1)],float)
    hull=ConvexHull(vertices)
    faces=hull.simplices.copy()
    cross=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]])
    flip=np.einsum('ij,ij->i',cross,vertices[faces[:,0]])<0
    faces[flip]=faces[flip][:,[0,2,1]]
    normal=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]])
    distance=np.einsum('ij,ij->i',normal,vertices[faces[:,0]])/np.linalg.norm(normal,axis=1)
    vertices*=.0002/distance.min()
    normal=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]])
    distance=np.einsum('ij,ij->i',normal,vertices[faces[:,0]])/np.linalg.norm(normal,axis=1)
    radius=np.linalg.norm(vertices,axis=1)
    require(distance.min()>=.0002-1e-18,'Kernel does not contain the entire planned clearance ball')
    native=md.Manifold(md.Mesh64(np.ascontiguousarray(vertices),np.ascontiguousarray(faces,dtype=np.uint64)))
    require(native.status()==md.Error.NoError,'Clearance kernel is not a closed convex manifold')
    return native,dict(vertices_m=vertices.tolist(),faces=faces.tolist(),
        intended_Euclidean_clearance_m=.0002,minimum_facet_distance_m=float(distance.min()),
        maximum_radial_distance_m=float(radius.max()),maximum_conservative_radial_excess_m=float(radius.max()-.0002),
        contains_planned_clearance_ball=True,geometry='circumscribed regular icosahedron')


def main(source,output,include_mucosa=False):
    began=time.monotonic()
    require(source.is_relative_to(WORK) and source.is_file() and output.is_relative_to(WORK)
            and not output.exists(),'Use an explicit saved construction and one new local stage')
    construction=json.loads(source.read_text())
    rows={r['key']:r for r in construction['meshes']}
    require(len(rows)==4 and construction.get('original_tv_geometry') is True,
            'Use the original curved central core, never the whole-tongue residual')
    require((source.parent/'anterior-cap-correction.json').is_file(),
            'The anterior cap must first query its actual unchanged-Z corridor')
    for row in rows.values():require(digest_file(row['mesh'])==row['sha256'],'Selected source changed')
    sl=Mesh.load(rows['superior-longitudinal']['mesh'])
    tv=Mesh.load(rows['transverse-vertical-shared']['mesh'])
    output.mkdir(parents=True)
    ker,kproof=kernel()
    atomic_json(output/'plan.json',dict(schema='tongue-tv-finite-clearance-conformance-plan/v1',
        source_construction=str(source),source_construction_sha256=digest_file(source),
        source_SL=rows['superior-longitudinal'],source_TV=rows['transverse-vertical-shared'],kernel=kproof,
        changed_anatomy=['TV reference core only'],unchanged=['SL','mucosa','aponeurosis','all fixed skeleton','paired IL','all extrinsics'],
        operation='Original curved TV core minus fixed finite tissue clearance exclusions',
        tissue_exclusions=(['actual literal SL','actual literal finite mucosa'] if include_mucosa else ['actual literal SL']),
        margins=dict(SL_free_gap_m=.0002,mucosal_free_gap_from_actual_finite_coat_m=.0002,
            earlier_mucosal_lower_bound_from_outer_carrier_m=.00062,maximum_modeled_mucosal_thickness_m=.00042,
            existing_coat_thickness_not_added_twice=True),
        broad_distortion_rejection=dict(maximum_removed_volume_fraction=.02,maximum_new_surface_source_distance_m=.001),
        exact_muscle_extent_claim=False,exclusive_histological_tissue_claim=False,
        all_original_assembly_checks_remain_required=True,installed=False,anatomical_acceptance=False))
    native_sl,native_tv=solid(sl),solid(tv)
    expanded=native_sl.minkowski_sum(ker)
    require(expanded.status()==md.Error.NoError,'Finite SL clearance expansion failed')
    expanded_mesh=export(expanded)
    expanded_mesh.save(output/'SL-clearance-exclusion.npz')
    atomic_json(output/'exclusion-checkpoint.json',dict(
        kernel=kproof,exclusion_mesh=str(output/'SL-clearance-exclusion.npz'),
        exclusion_sha256=digest_file(output/'SL-clearance-exclusion.npz'),
        source_SL_geometry_sha256=sl.geometry_sha,source_SL_modified=False,
        faces=len(expanded_mesh.faces),elapsed_seconds=time.monotonic()-began))
    if include_mucosa:
        mucosa_path=corridor.SOURCE/'part-00.npz'
        mucosa=Mesh.load(mucosa_path)
        expanded_mucosa=solid(mucosa).minkowski_sum(ker)
        require(expanded_mucosa.status()==md.Error.NoError,'Finite mucosal clearance expansion failed')
        export(expanded_mucosa).save(output/'mucosa-clearance-exclusion.npz')
        atomic_json(output/'mucosa-exclusion-checkpoint.json',dict(
            kernel=kproof,source_mucosa=str(mucosa_path),source_mucosa_sha256=digest_file(mucosa_path),
            source_mucosa_geometry_sha256=mucosa.geometry_sha,source_mucosa_modified=False,
            exclusion_mesh=str(output/'mucosa-clearance-exclusion.npz'),
            exclusion_sha256=digest_file(output/'mucosa-clearance-exclusion.npz'),
            free_gap_from_finite_coat_m=.0002,elapsed_seconds=time.monotonic()-began))
        result=(native_tv-expanded)-expanded_mucosa
    else:
        result=native_tv-expanded
    require(result.status()==md.Error.NoError,'Finite TV difference failed')
    removed=float(native_tv.volume()-result.volume())
    fraction=removed/native_tv.volume()
    require(len(result.decompose())==1,'Clearance conformance disconnected the TV reference core')
    require(0<=fraction<=.02,'Clearance conformance removed too much of the explicitly chosen TV core')
    generated=export(result)
    generated.save(output/'TV-native-Float64-before-storage.npz')
    # Audit exactly the eventual Float32 coordinates, with the same dyadic Y origin.
    vertices=(generated.vertices-ORIGIN).astype('<f4').astype(float)+ORIGIN
    mesh=Mesh(vertices,generated.faces,FRAME)
    topo=topology(mesh)
    require(topo['closed_oriented_manifold'] and len(topo['components'])==1,
            'Stored conformance topology is not one closed manifold')
    nearest,_=BVH(tv).closest_points(mesh.vertices,backend='bvh')
    deviation=float(nearest['distance_m'].max())
    require(deviation<=.001,'Finite conformance broadly changed the source curved envelope')
    t=mesh.triangles
    raw=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    normals=np.zeros_like(mesh.vertices)
    for j in range(3):np.add.at(normals,mesh.faces[:,j],raw)
    normals/=np.maximum(np.linalg.norm(normals,axis=1)[:,None],1e-30)
    words=np.rint(normals*32767).astype('<i2')
    mesh=Mesh(mesh.vertices,mesh.faces,FRAME,
        dict(normals=words,normalSigns=np.ones(len(vertices),dtype='i1')),
        dict(operation='finite-conformance-to-actual-SL-clearance',histological_tissue_claim=False))
    raw_path=output/'transverse-vertical-shared-before-corner-chart.npz'
    mesh.save(raw_path)
    target=output/'transverse-vertical-shared.npz'
    tvrow=corridor.chart_il(raw_path,target)
    result_rows=[dict(rows[key]) for key in ('superior-longitudinal','inferior-longitudinal-left','inferior-longitudinal-right')]
    result_rows.append(tvrow)
    proof=dict(schema='tongue-tv-finite-clearance-conformance/v1',
        source_construction=str(source),source_construction_sha256=digest_file(source),
        meshes=result_rows,kernel=kproof,removed_volume_m3=removed,removed_volume_fraction=float(fraction),
        source_volume_m3=float(native_tv.volume()),result_volume_m3=float(result.volume()),
        native_Float64_mesh=str(output/'TV-native-Float64-before-storage.npz'),
        native_Float64_sha256=digest_file(output/'TV-native-Float64-before-storage.npz'),
        both_SL_and_finite_mucosa_excluded=include_mucosa,
        source_surface_deviation_maximum_m=deviation,
        triangles_before=len(tv.faces),triangles_after=len(mesh.faces),added_facets=len(mesh.faces)-len(tv.faces),
        components=topo['components'],topology=topo,
        unchanged_SL_geometry_sha256=sl.geometry_sha,
        source_files_modified=False,skeleton_modified=False,SL_and_covering_shared_charts_unchanged=True,
        tissue_ownership_claim=False,exact_muscle_extent_claim=False,numerical_checks_pending=True,
        installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'construction.json',proof)
    print(json.dumps({k:proof[k] for k in ('removed_volume_fraction','source_surface_deviation_maximum_m','triangles_after','elapsed_seconds')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--include-mucosa',action='store_true')
    args=parser.parse_args()
    main(args.source.resolve(),args.output.resolve(),args.include_mucosa)
