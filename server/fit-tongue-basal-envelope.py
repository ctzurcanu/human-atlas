"""Explicit protected basal tongue cage; fixed bones and sources untouched."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.fit import cage_mesh
from geometry_tools.audit import audit,topology
from geometry_tools.contact import pair_contacts

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'


def main(output):
    began=time.monotonic();output=output.resolve()
    if output.exists() or not output.is_relative_to(WORK) or output==WORK:
        raise ValueError('Use a new local tongue workbench stage')
    source=WORK/'source-v1/meshes/native-00.npz';native=Mesh.load(source)
    if digest_file(native.source['catalogue_path'])!=native.source['catalogue_sha256']:
        raise ValueError('Current catalogue changed')
    # Resolve exact sources from the audited inventory, never guessed indices.
    source_audit=json.loads((WORK/'source-v1/source-audit.json').read_text())
    rows=[r for r in source_audit['contextParts'] if r['sourcePart'].get('system')=='skeletal']
    bones=[]
    for row in rows:
        if digest_file(row['meshSavePath'])!=row['meshSaveSha256']:
            raise ValueError('Stale audited bone export')
        bone=Mesh.load(Path(row['meshSavePath']))
        if bone.geometry_sha!=row['geometrySha256']:
            raise ValueError('Bone geometry differs from inventory')
        bones.append(bone)
    hyoid=[m for m in bones if m.source['record']['id'].startswith('TA98HYOID:')]
    if len(bones)!=8 or len(hyoid)!=5:
        raise ValueError('Expected the eight declared bone sources, including five fixed hyoid parts')
    z0,z1=native.bounds[:,2]
    centres=np.zeros((2,2));cov=np.tile(np.diag([.000225,.0001]),(2,1,1))
    spec=dict(source_bounds=[float(z0),float(z1)],target_bounds=[float(z0),float(z1)],
        fractions=[0.,1.],source_centers=centres.tolist(),target_centers=centres.tolist(),
        source_covariances=cov.tolist(),target_covariances=cov.tolist(),minimum_stretch=0.,
        protected_basal_lift=dict(amplitude_m=.006,chart_axes=['Z','X','Y'],z_fall_bounds_m=[.040,.055],y_fall_bounds_m=[1.610,1.633]))
    output.mkdir(parents=True)
    plan=dict(schema='tongue-protected-basal-cage/v1',source=str(source),source_sha256=digest_file(source),
        source_catalogue_sha256=native.source['catalogue_sha256'],source_geometry_sha256=native.geometry_sha,
        source_bank_sha256=native.source['bank_sha256'],model='local-ta98',whole_organ='ZA:Tongue',
        operation='Proper rigid chart to Z/X/Y, explicit LongitudinalOrganCage, inverse rigid chart',
        cage=spec,hyoid_upper_y_m=max(float(m.vertices[:,1].max()) for m in hyoid),
        moving=['Derived whole-tongue envelope; dependent region/layer candidates require regeneration'],
        fixed=[m.source['record']['id'] for m in bones],
        adjoining_interfaces=[dict(parts=['Genioglossus','Hyoglossus','Styloglossus'],status='Retain bony origins; intralingual insertion checks pending after whole correction')],
        anterior_lock_z_m=.055,dorsal_lock_y_m=1.633,
        assessment=str(WORK/'native-basal-cage-assessment-v2/input-plan.json'),
        bone_geometry_modified=False,source_file_modified=False,anatomical_acceptance=False,installed=False,wall_seconds=45)
    atomic_json(output/'fit-input.json',plan)
    # Cyclic coordinate permutation is a proper rotation, determinant +1.
    order=[2,0,1];back=[1,2,0]
    attrs=dict(native.attributes);attrs['normals']=attrs['normals'][:,order]
    chart=Mesh(native.vertices[:,order],native.faces,'tongue-basal-ZXY-metres',attrs,dict(native.source))
    fitted,metrics=cage_mesh(chart,spec)
    positions=fitted.vertices[:,back]
    origin=np.array([0.,1.5,0.]);encoded=(positions-origin).astype('<f4').astype(float)+origin
    normals=fitted.attributes['normals'][:,back]
    locked=(native.vertices[:,1]>=1.633)|(native.vertices[:,2]>=.055)
    normals[locked]=native.attributes['normals'][locked]
    result=Mesh(encoded,native.faces,native.frame,{**native.attributes,'normals':normals},dict(native.source,
        kind='protected-basal-cage-candidate',derived_from_geometry=native.geometry_sha,cage_input_sha256=digest_file(output/'fit-input.json')))
    path=output/'whole-tongue.npz';result.save(path)
    checks=audit(result,require_closed=True,max_pairs=100000,max_tests=1000000)
    contacts=[]
    for bone in bones:
        _,s=pair_contacts(result,bone,0,max_pairs=100000,max_tests=1000000)
        contacts.append(dict(part_id=bone.source['record']['id'],geometry_sha256=bone.geometry_sha,summary=s))
    if time.monotonic()-began>=45:raise TimeoutError('Basal fit exceeded45-second worker budget')
    bone_clear=all(r['summary']['complete_within_pair_budget'] and r['summary']['crossing_face_pairs']==r['summary']['coplanar_overlap_face_pairs']==0 for r in contacts)
    front=native.vertices[:,2]>=.055
    report=dict(**plan,mesh=str(path),mesh_sha256=digest_file(path),geometry_sha256=result.geometry_sha,
        passed=bool(checks['passed'] and bone_clear),audit=checks,bone_contacts=contacts,
        cage_metrics=metrics,global_jacobian_lower_bound=1-1.5*.006/.023,
        maximum_vertex_motion_m=float(np.linalg.norm(encoded-native.vertices,axis=1).max()),
        maximum_position_encoding_error_m=float(np.linalg.norm(encoded-positions,axis=1).max()),
        anterior_position_error_m=float(np.linalg.norm(encoded[front]-native.vertices[front],axis=1).max()),
        protected_position_error_m=float(np.linalg.norm(encoded[locked]-native.vertices[locked],axis=1).max()),
        before_volume_m3=topology(native)['oriented_triangle_integral_m3'],after_volume_m3=topology(result)['oriented_triangle_integral_m3'],
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'fit-result.json',report)
    print(dict(passed=report['passed'],self_passed=checks['passed'],bone_crossings=sum(r['summary']['crossing_face_pairs'] for r in contacts),global_jacobian_lower_bound=report['global_jacobian_lower_bound'],maximum_motion_m=report['maximum_vertex_motion_m'],elapsed_seconds=report['elapsed_seconds']),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    main(parser.parse_args().output)
