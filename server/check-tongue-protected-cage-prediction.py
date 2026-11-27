"""Predict the exact supported 6 mm cage on native controls; emit evidence only."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
ORIGIN=np.array([0.,1.5,0.]);LIFT=.006


def down(a,full,stop):
    t=np.clip((a-full)/(stop-full),0.,1.)
    return 1-t*t*(3-2*t)


def main(output):
    began=time.monotonic();deadline=began+25;output=output.resolve()
    assert output.is_relative_to(WORK) and not output.exists()
    source_path=WORK/'source-v1/source-audit.json';source=json.loads(source_path.read_text())
    catalogue_path=Path(source['sourceCatalogues']['male']['path']);catalogue_sha=digest_file(catalogue_path)
    assert catalogue_sha==source['sourceCatalogues']['male']['sha256']
    atlas=json.loads(catalogue_path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
    native=reader.read_part(atlas,parts['ZA:Tongue'],frame=frame,catalogue_path=catalogue_path,catalogue_sha=catalogue_sha)
    assert len(native.vertices)==343 and len(native.faces)==682 and native.geometry_sha==source['contextParts'][0]['geometrySha256']
    moved=native.vertices.copy();moved[:,1]+=LIFT*down(moved[:,2],.040,.055)*down(moved[:,1],1.610,1.633)
    actual_positions=(moved-ORIGIN).astype('<f4').astype(float)+ORIGIN
    predicted=Mesh(actual_positions,native.faces.copy(),frame,source=dict(kind='in-memory-protected-cage-prediction',not_an_authoritative_source=True))
    protected=(native.vertices[:,2]>=.055)|(native.vertices[:,1]>=1.633)
    protected_error=float(np.linalg.norm(actual_positions[protected]-native.vertices[protected],axis=1).max())
    xz_error=float(np.max(abs(actual_positions[:,[0,2]]-native.vertices[:,[0,2]])))
    storage_error=float(np.linalg.norm(actual_positions-moved,axis=1).max())
    check=audit(predicted,require_closed=True,expected_components=1,expected_genus=[0],max_pairs=250000,max_tests=2000000)
    before=np.cross(native.triangles[:,1]-native.triangles[:,0],native.triangles[:,2]-native.triangles[:,0])
    after=np.cross(predicted.triangles[:,1]-predicted.triangles[:,0],predicted.triangles[:,2]-predicted.triangles[:,0])
    cosine=np.einsum('ij,ij->i',before,after)/(np.linalg.norm(before,axis=1)*np.linalg.norm(after,axis=1))
    output.mkdir(parents=True);tree=BVH(predicted);rows=[]
    for row in source['contextParts']:
        if row['sourcePart'].get('system')!='skeletal':continue
        if time.monotonic()>=deadline:raise TimeoutError('Protected cage prediction reached its 25 s internal deadline')
        assert parts[row['partId']]==row['sourcePart'] and digest_file(row['bankPath'])==row['bankSha256']
        bone=reader.read_part(atlas,parts[row['partId']],frame=frame,catalogue_path=catalogue_path,catalogue_sha=catalogue_sha)
        assert bone.geometry_sha==row['geometrySha256'];bone_tree=BVH(bone)
        arrays,summary=pair_contacts(predicted,bone,0.,first_bvh=tree,second_bvh=bone_tree,max_pairs=50000,max_tests=1000000)
        bad=np.isin(arrays['relation'],['crossing','coplanar-overlap']);pairs=arrays['pairs'][bad]
        name=row['partId'].replace(':','_').replace(' ','_').replace('/','_');path=output/(name+'-witnesses.npz')
        atomic_npz(path,**arrays,defect_pairs=pairs,defect_first_triangles=predicted.triangles[pairs[:,0]],defect_second_triangles=bone.triangles[pairs[:,1]])
        nearest=tree.closest_mesh(bone_tree)
        rows.append(dict(bone=row['partId'],bone_bank_path=bone.source['bank_path'],bone_bank_sha256=bone.source['bank_sha256'],bone_geometry_sha256=bone.geometry_sha,contact=summary,minimum_distance_m=float(nearest['distance_m']),nearest=nearest,witness_path=str(path),witness_sha256=digest_file(path)))
    assert len(rows)==8 and digest_file(catalogue_path)==catalogue_sha
    report=dict(schema='tongue-protected-cage-predictive-check/v1',source_catalogue_path=str(catalogue_path),source_catalogue_sha256=catalogue_sha,source_audit_path=str(source_path),source_audit_sha256=digest_file(source_path),source_part_id='ZA:Tongue',source_bank_path=native.source['bank_path'],source_bank_sha256=native.source['bank_sha256'],source_geometry_sha256=native.geometry_sha,predicted_geometry_sha256=predicted.geometry_sha,frame=frame,
        field=dict(lift_m=LIFT,z_full_through_m=.040,z_zero_from_m=.055,y_full_through_m=1.610,y_zero_from_m=1.633,transition='C1 cubic smoothstep',analytic_dy_derivative_lower_bound=1-LIFT*1.5/.023),
        actual_storage=dict(origin_m=ORIGIN.tolist(),position_type='Float32',maximum_position_roundtrip_error_m=storage_error,protected_vertices=int(protected.sum()),maximum_protected_coordinate_error_m=protected_error,all_vertices_maximum_xz_error_m=xz_error,all_faces_and_indices_literal=np.array_equal(predicted.faces,native.faces)),geometry_audit=check,bone_checks=rows,minimum_source_face_normal_cosine=float(np.min(cosine)),nonpositive_source_face_normal_cosines=int(np.count_nonzero(cosine<=0)),
        passed=check['passed'] and protected_error==0 and xz_error==0 and all(r['contact']['complete_within_pair_budget'] and r['contact']['crossing_face_pairs']==r['contact']['coplanar_overlap_face_pairs']==0 for r in rows),
        scope='Predictive numerical check of the proposed supported field on actual source vertices and actual Float32-decoded triangle geometry. No authoritative mesh or normal field is emitted; no surface/section or anatomical acceptance is inferred.',mesh_emitted=False,atlas_modified=False,skeleton_modified=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'check.json',report)
    print(json.dumps(dict(output=str(output),passed=report['passed'],actual_storage=report['actual_storage'],self_defects=check['self_intersections']['defect_pairs'],bone_checks=[dict(bone=r['bone'],crossing=r['contact']['crossing_face_pairs'],coplanar=r['contact']['coplanar_overlap_face_pairs'],minimum_distance_mm=r['minimum_distance_m']*1000,complete=r['contact']['complete_within_pair_budget']) for r in rows],minimum_source_face_normal_cosine=report['minimum_source_face_normal_cosine'],elapsed_seconds=report['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();main(args.output)
