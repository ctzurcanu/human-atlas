"""Retriangulate one proven Float32 sliver star with literal boundary controls.

Only a closed local disk around the explicit subnanometre edge is replaced.
Every source coordinate and all faces outside that disk remain literal. No
anatomical radius, extent, clearance, fixed covering or skeleton is changed.
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
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH

spec=importlib.util.spec_from_file_location('tongue_conform',ROOT/'server/conform-tongue-tv-clearance.py')
conform=importlib.util.module_from_spec(spec)
spec.loader.exec_module(conform)
WORK,FRAME=conform.WORK,conform.FRAME


def require(ok,message):
    if not ok:raise ValueError(message)


def main(source,output,edge):
    began=time.monotonic()
    require(source.is_relative_to(WORK) and source.is_file() and output.is_relative_to(WORK)
            and not output.exists(),'Use an explicit saved source and new stage')
    construction=json.loads(source.read_text())
    diagnosis_path=source.parent/'literal-storage-diagnosis.json'
    diagnosis=json.loads(diagnosis_path.read_text())
    require(diagnosis['native_audit']['passed'] and not diagnosis['stored_audit']['passed']
            and diagnosis['literal_face_indices'] and diagnosis['maximum_actual_float32_position_error_m']<=1e-8,
            'Need literal same-face native-pass/stored-fail evidence, not an unexplained anatomy defect')
    source_path=source.parent/'transverse-vertical-shared-before-corner-chart.npz'
    mesh=Mesh.load(source_path)
    native=Mesh.load(construction['native_Float64_mesh'])
    require(np.array_equal(mesh.faces,native.faces),'Native/stored face correspondence changed')
    length=float(np.linalg.norm(mesh.vertices[edge[0]]-mesh.vertices[edge[1]]))
    require(length<=1e-9,'Only the explicit subnanometre numerical edge may be repaired here')
    patch_ids=np.flatnonzero(np.any(np.isin(mesh.faces,edge),axis=1))
    pairs=diagnosis['stored_audit']['self_intersections']['examples']
    require(len(pairs)==3 and set(sum(pairs,[]))<=set(patch_ids.tolist()),
            'All observed self defects must belong to this one explicitly bounded disk')
    patch=mesh.faces[patch_ids]
    directed=np.r_[patch[:,[0,1]],patch[:,[1,2]],patch[:,[2,0]]]
    keys,inverse,count=np.unique(np.sort(directed,axis=1),axis=0,return_inverse=True,return_counts=True)
    boundary_edges=directed[count[inverse]==1]
    next_vertex={int(a):int(b) for a,b in boundary_edges}
    require(len(next_vertex)==len(boundary_edges) and len(set(next_vertex.values()))==len(next_vertex),
            'Patch outer boundary is not one simple directed polygon')
    ring=[]
    start=int(boundary_edges.min())
    current=start
    while current not in ring:
        ring.append(current)
        current=next_vertex[current]
    require(current==start and len(ring)==len(boundary_edges),
            'Patch has multiple boundary loops or an open boundary')
    require(not set(edge)&set(ring),'The tiny edge is not entirely internal to the preserved patch boundary')
    vertices=mesh.vertices[ring]
    center=vertices.mean(0)
    _,_,axes=np.linalg.svd(vertices-center,full_matrices=False)
    projection=(vertices-center)@axes[:2].T
    nxt=np.roll(projection,-1,axis=0)
    signed=float((projection[:,0]*nxt[:,1]-projection[:,1]*nxt[:,0]).sum()/2)
    if signed<0:projection[:,1]*=-1
    local=np.asarray(md.triangulate([np.asarray(projection,dtype=float)]),dtype=np.int64)
    require(local.shape==(len(ring)-2,3),'Triangulator did not preserve the simple patch polygon')
    new_faces=np.asarray(ring,dtype=np.int64)[local]
    old_raw=np.cross(mesh.triangles[patch_ids,1]-mesh.triangles[patch_ids,0],
                     mesh.triangles[patch_ids,2]-mesh.triangles[patch_ids,0])
    proposed=mesh.vertices[new_faces]
    new_raw=np.cross(proposed[:,1]-proposed[:,0],proposed[:,2]-proposed[:,0])
    if np.einsum('ij,j->i',new_raw,old_raw.sum(0)).sum()<0:new_faces=new_faces[:,[0,2,1]]
    keep=np.ones(len(mesh.faces),bool)
    keep[patch_ids]=False
    faces=np.r_[mesh.faces[keep],new_faces]
    result=Mesh(mesh.vertices,faces,FRAME)
    topo=topology(result)
    require(topo['closed_oriented_manifold'] and len(topo['components'])==1
            and topo['components'][0]['closed_orientable_genus']==0,
            'Literal-coordinate patch did not preserve the one closed genus-zero boundary')
    require(np.array_equal(result.vertices,mesh.vertices) and np.array_equal(result.faces[:keep.sum()],mesh.faces[keep]),
            'Patch repair moved coordinates or changed unrelated face words')
    tri=result.triangles
    raw=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    normals=np.zeros_like(result.vertices)
    for j in range(3):np.add.at(normals,result.faces[:,j],raw)
    normals/=np.maximum(np.linalg.norm(normals,axis=1)[:,None],1e-30)
    words=np.rint(normals*32767).astype('<i2')
    raw_mesh=Mesh(result.vertices,result.faces,FRAME,
        dict(normals=words,normalSigns=np.ones(len(result.vertices),dtype='i1')),
        dict(local_patch_retriangulation=True,all_source_positions_literal=True))
    output.mkdir(parents=True)
    raw_path=output/'transverse-vertical-shared-before-corner-chart.npz'
    raw_mesh.save(raw_path)
    tvrow=conform.corridor.chart_il(raw_path,output/'transverse-vertical-shared.npz')
    rows=[dict(r) for r in construction['meshes'] if r['key']!='transverse-vertical-shared']+[tvrow]
    source_core=Mesh.load(construction['source_construction'] and
        json.loads(Path(construction['source_construction']).read_text())['meshes'][-1]['mesh'])
    samples=np.r_[result.vertices[np.unique(new_faces)],result.triangles[-len(new_faces):].mean(1)]
    nearest,_=BVH(source_core).closest_points(samples,backend='bvh')
    maximum=max(construction['source_surface_deviation_maximum_m'],float(nearest['distance_m'].max()))
    source_volume=construction['source_volume_m3']
    volume_delta=topo['oriented_triangle_integral_m3']-topology(mesh)['oriented_triangle_integral_m3']
    removed_fraction=construction['removed_volume_fraction']-volume_delta/source_volume
    require(maximum<=.001 and 0<=removed_fraction<=.02,'Literal local patch exceeded the unchanged source deviation/volume gates')
    proof=dict(schema='tongue-TV-boundary-preserving-storage-patch/v1',source_construction=str(source),
        source_construction_sha256=digest_file(source),storage_diagnosis=str(diagnosis_path),
        storage_diagnosis_sha256=digest_file(diagnosis_path),source_stored_mesh=str(source_path),
        source_stored_sha256=digest_file(source_path),explicit_internal_edge=edge,stored_edge_length_m=length,
        replaced_faces=patch_ids.tolist(),preserved_boundary_vertex_ids=ring,
        preserved_boundary_coordinates=vertices.tolist(),new_faces=new_faces.tolist(),
        source_positions_all_literal=True,unrelated_faces_all_literal=True,
        removed_internal_controls=edge,controls_moved_or_added=0,
        patch_projection_maximum_residual_m=float(np.abs((vertices-center)@axes[2]).max()),
        faces_before=len(mesh.faces),faces_after=len(result.faces),volume_delta_m3=float(volume_delta),
        geometry_sha256=result.geometry_sha,remaining_check_scope='All original self/interface/covering/bone/containment gates plus planned clearance',
        installed=False,anatomical_acceptance=False)
    atomic_json(output/'patch-proof.json',proof)
    updated=dict(construction, schema='tongue-TV-storage-repaired-construction/v1',meshes=rows,
        predecessor_construction=str(source),predecessor_construction_sha256=digest_file(source),
        patch_proof=str(output/'patch-proof.json'),patch_proof_sha256=digest_file(output/'patch-proof.json'),
        source_surface_deviation_maximum_m=maximum,removed_volume_fraction=removed_fraction,
        triangles_after=len(result.faces),numerical_checks_pending=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'construction.json',updated)
    print(json.dumps(dict(triangles=len(result.faces),replaced_faces=len(patch_ids),new_faces=len(new_faces),
        positions_moved=0,removed_volume_fraction=removed_fraction,source_deviation_m=maximum,
        elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--edge',type=int,nargs=2,required=True)
    args=parser.parse_args()
    main(args.source.resolve(),args.output.resolve(),args.edge)
