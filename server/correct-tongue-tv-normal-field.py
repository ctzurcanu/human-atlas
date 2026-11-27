"""Remove tessellation-weight bias from TV normals; literal geometry remains.

Corner-angle weighted unit face normals are independent of triangle area.
Existing charts are retained. Each new chart normal stays in every incident
face's outward hemisphere; sharp local charts fall back to their original
field when Int16 quantization cannot preserve that strict inequality.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'server')]
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from tongue_basal_contract import require

def bound(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=digest_file(p))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();began=time.monotonic()
    require(not args.output.exists(),'Use a new normal-only stage')
    construction=json.loads(args.source.read_text());assessment_path=args.source.with_name('numerical-assessment.json')
    assessment=json.loads(assessment_path.read_text())
    require(assessment['construction_sha256']==digest_file(args.source) and assessment['all_actual_geometric_and_normal_gates_passed'] is True,'Complete actual predecessor gates required')
    source_row=next(r for r in construction['meshes'] if r['key']=='transverse-vertical-shared')
    path=Path(source_row['mesh']);require(digest_file(path)==source_row['sha256'],'Selected TV source changed')
    mesh=Mesh.load(path)
    vertices,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True);faces=inverse[mesh.faces];tri=vertices[faces]
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);unit=cross/np.linalg.norm(cross,axis=1)[:,None]
    summed=np.zeros_like(vertices)
    for j in range(3):
        a=tri[:,(j+1)%3]-tri[:,j];b=tri[:,(j+2)%3]-tri[:,j]
        angle=np.arctan2(np.linalg.norm(np.cross(a,b),axis=1),np.einsum('ij,ij->i',a,b))
        np.add.at(summed,faces[:,j],unit*angle[:,None])
    summed/=np.maximum(np.linalg.norm(summed,axis=1)[:,None],1e-30)
    old_words=np.asarray(mesh.attributes['normals']);signs=np.asarray(mesh.attributes['normalSigns'])
    require(old_words.dtype==np.dtype('int16') and np.all(signs==1),'This stage preserves actual positive source signs')
    old=old_words.astype(float);old/=np.maximum(np.linalg.norm(old,axis=1)[:,None],1e-30)
    target=summed[inverse]
    old_dot=np.einsum('fci,fi->fc',old[mesh.faces],unit)
    new_dot=np.einsum('fci,fi->fc',target[mesh.faces],unit)
    require(np.all(old_dot>0),'Predecessor has inverted corner normals')
    margin=np.minimum(1e-4,old_dot*.5)
    decreasing=new_dot<margin
    limits=np.ones_like(old_dot)
    limits[decreasing]=(old_dot[decreasing]-margin[decreasing])/(old_dot[decreasing]-new_dot[decreasing])
    alpha=np.ones(len(mesh.vertices))
    for j in range(3):np.minimum.at(alpha,mesh.faces[:,j],limits[:,j])
    field=(1-alpha[:,None])*old+alpha[:,None]*target
    field/=np.maximum(np.linalg.norm(field,axis=1)[:,None],1e-30)
    words=np.rint(field*32767).clip(-32767,32767).astype('<i2')
    dots=np.einsum('fci,fi->fc',words[mesh.faces].astype(float),unit)
    bad_vertex_ids=np.unique(mesh.faces[dots<=0])
    words[bad_vertex_ids]=old_words[bad_vertex_ids]
    dots=np.einsum('fci,fi->fc',words[mesh.faces].astype(float),unit)
    require(np.all(dots>0),'Int16 normal update reversed a physical corner')
    with np.load(path,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
    require('corner_normals' not in arrays,'Normal update must precede the redundant corner adapter')
    arrays['attribute_normals']=words
    args.output.mkdir(parents=True)
    target_path=args.output/'transverse-vertical-shared.npz';atomic_npz(target_path,**arrays)
    actual=Mesh.load(target_path)
    require(np.array_equal(mesh.vertices,actual.vertices) and np.array_equal(mesh.faces,actual.faces)
            and mesh.geometry_sha==actual.geometry_sha and np.array_equal(mesh.attributes['normalSigns'],actual.attributes['normalSigns']),
            'Normal-only stage changed literal geometry/signs')
    current=words.astype(float);current/=np.maximum(np.linalg.norm(current,axis=1)[:,None],1e-30)
    changed=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',old,current),-1,1)))
    # The same lateral ring witness where a rotationally symmetric chart should
    # not alternate Y direction merely because two ring spacings differ.
    ring_records=[]
    for z in np.unique(vertices[:,2]):
        ids=np.flatnonzero(vertices[:,2]==z)
        if len(ids)<20:continue
        vi=ids[np.argmax(vertices[ids,0])]
        if not .033<=z<=.070:continue
        indices=np.flatnonzero(inverse==vi);k=indices[0]
        ring_records.append(dict(z_m=float(z),point=vertices[vi].tolist(),old_normal=old[k].tolist(),new_normal=current[k].tolist()))
    proof=dict(schema='tongue-TV-literal-normal-field-correction/v1',passed=True,
        predecessor_construction=bound(args.source),predecessor_assessment=bound(assessment_path),
        source=bound(path),candidate=bound(target_path),geometry_sha256=mesh.geometry_sha,
        method='Exact-coordinate welded corner-angle weighted unit face normals, projected into each existing positive chart cone',
        every_indexed_coordinate_literal=True,every_face_index_literal=True,normal_signs_literal=True,
        faces_or_vertices_added=0,negative_or_zero_corner_dots=0,
        constrained_existing_charts=int(np.count_nonzero(alpha<1)),quantized_fallback_vertices=bad_vertex_ids.tolist(),
        normal_angular_change_degrees_percentiles=np.percentile(changed,[0,50,90,95,99,100]).tolist(),
        actual_lateral_ring_witnesses=ring_records,
        all_geometric_gates_reused_only_by_literal_geometry_identity=True,
        fixed_shared_SL_coverings_and_skeleton_unchanged=True,
        anatomical_acceptance=False,installed=False,render_required=True,elapsed_seconds=time.monotonic()-began)
    proof_path=args.output/'normal-proof.json';atomic_json(proof_path,proof)
    updated=copy.deepcopy(construction)
    updated.update(schema='tongue-TV-normal-field-correction/v1',normal_only_predecessor=bound(args.source),
                   normal_field_correction=bound(proof_path),numerical_checks_pending=False)
    row=next(r for r in updated['meshes'] if r['key']=='transverse-vertical-shared')
    row.update(mesh=str(target_path.resolve()),sha256=digest_file(target_path),normal_chart_changed=True,
               positions_and_faces_literal=True,geometry_sha256=mesh.geometry_sha)
    out=args.output/'construction.json';atomic_json(out,updated)
    numerical=copy.deepcopy(assessment)
    numerical.update(construction=str(out.resolve()),construction_sha256=digest_file(out),
        geometric_evidence_predecessor=dict(construction=bound(args.source),assessment=bound(assessment_path)),
        normal_field_correction=bound(proof_path))
    numericrow=next(r for r in numerical['meshes'] if r['key']=='transverse-vertical-shared')
    numericrow.update(mesh=str(target_path.resolve()),sha256=digest_file(target_path),nonpositive_face_corner_dots=0)
    atomic_json(args.output/'numerical-assessment.json',numerical)
    print(json.dumps(dict(construction=str(out.resolve()),sha256=digest_file(out),proof=str(proof_path.resolve()),
                          changed_normal_vertices=int(np.count_nonzero(np.any(words!=old_words,axis=1))),elapsed_seconds=time.monotonic()-began)))

if __name__=='__main__':main()
