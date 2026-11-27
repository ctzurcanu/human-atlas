"""Independent saved-array gates; mirror transfer never repeats the fit."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit,winding_contains
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.annular_pipe_envelopes import close_with_shared_disks
from geometry_tools.algorithms.simple_compartment_section import simple_compartment_section

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1'
OUT=(ROOT/sys.argv[1]).resolve() if len(sys.argv)>1 else BASE/'bilateral-course-proof-v1'
def tangents(c):
    t=np.gradient(c,axis=0);return t/np.linalg.norm(t,axis=1)[:,None]
def words(n,rotation):
    q=n.astype(float)@rotation.T;q/=np.linalg.norm(q,axis=1)[:,None];return np.rint(q*32767).astype('<i2')

def main():
    start=time.monotonic();plan=json.loads((OUT/'plan.json').read_bytes());path=ROOT/plan['left_candidate'];assert digest_file(path)==plan['left_candidate_sha256'];c=np.load(path);proof=dict(schema='native-lcfa-bilateral-saved-array-proof/v1',plan_sha256=digest_file(OUT/'plan.json'),left_candidate_sha256=digest_file(path),passed=False,anatomical_acceptance=False,installed=False)
    try:
        frame=str(c['frame']);v=c['vertices'];f=c['faces'];rings=c['rings'];fixed=c['fixed_station_indices'];old=c['source_vertices'];rot=c['station_rotations']
        catalogue=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';atlas=json.loads(catalogue.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
        def source_attributes(side,positions):
            attributes={}
            for pid in ['ZA:Lateral circumflex femoral artery.'+side,'TA98NATIVE:male:LCFA:ascending:'+side]:
                mesh=reader.read_part(atlas,parts[pid],frame=frame)
                for point,n,u in zip(mesh.vertices,mesh.attributes['normals'],mesh.attributes['uvs']):
                    key=tuple(point)
                    if key in attributes:assert np.array_equal(attributes[key][0],n)and np.array_equal(attributes[key][1],u)
                    attributes[key]=(n,u)
            rows=[attributes[tuple(point)]for point in positions]
            return np.array([r[0]for r in rows],dtype='<i2'),np.array([r[1]for r in rows],dtype='<f4')
        original_normal,uv=source_attributes('l',old);normal=original_normal.copy()
        for j,ring in enumerate(rings):
            if j not in fixed:normal[ring]=words(original_normal[ring],rot[j])
        assert np.array_equal(normal[rings[fixed]],original_normal[rings[fixed]])
        proof['attribute_source_mapping']='Exact current-source world-coordinate identity; original normals/UVs remapped into welded geometry order without averaging.'
        proof['attribute_catalogue_sha256']=digest_file(catalogue)
        oldclosed,oldfaces,_=close_with_shared_disks(old,f,expected_mouths=2);newclosed,newfaces,_=close_with_shared_disks(v,f,expected_mouths=2)
        classifiers=[Mesh(oldclosed,oldfaces,frame),Mesh(newclosed,newfaces,frame)]
        assert all(audit(m,require_closed=True,expected_components=1,expected_genus=[0])['passed']for m in classifiers)
        oldcenters=old[rings].mean(1);newcenters=v[rings].mean(1);normals=[tangents(oldcenters),tangents(newcenters)];sections=[]
        for station in plan['section_stations']:
            values=[simple_compartment_section(m,centers[station],n[station])for m,centers,n in zip(classifiers,[oldcenters,newcenters],normals)]
            change=abs(values[1]['area_m2']-values[0]['area_m2'])/values[0]['area_m2'];sections.append(dict(station=station,source_actual_cut_area_m2=values[0]['area_m2'],candidate_actual_cut_area_m2=values[1]['area_m2'],relative_change=change))
            assert change<=plan['maximum_actual_outer_section_area_change_relative']
        proof['actual_outer_calibre_sections']=sections
        mirror_path=ROOT/plan['mirror_contract'];assert digest_file(mirror_path)==plan['mirror_sha256'];mirror=np.load(mirror_path);right=np.load(BASE/'source-gate-v1/r-component-0.npz');rv=right['vertices'].copy();rn,ru=source_attributes('r',rv);S=np.array([-1.,1,1]);lookup={tuple(p):i for i,p in enumerate(rv)}
        left_used=np.unique(f);right_ids=np.array([lookup[tuple(old[i]*S)]for i in left_used]);rv[right_ids]=v[left_used]*S;rn[right_ids]=(normal[left_used].astype(np.int32)*S.astype(int)).astype('<i2')
        sides=[]
        for side,p,faces,n,u,bonefile in [('l',v,f,normal,uv,BASE/'femur-whole-reference-v1/whole-femur-classifier.npz'),('r',rv,right['faces'],rn,ru,BASE/'right-mirror-contract-v1/right-whole-femur-classifier.npz')]:
            mesh=Mesh(p,faces,frame,dict(normals=n,uvs=u));bone_data=np.load(bonefile);bone=Mesh(bone_data['vertices'],bone_data['faces'],frame);q,auditproof=pair_contacts(bone,mesh,0.);assert auditproof['complete_within_pair_budget']and not len(q['pairs']);used=np.unique(faces);membership=winding_contains(bone,p[used],distance_backend='open3d');assert not any(membership['inside'])and not any(membership['boundary_ambiguous']);near,_=SurfaceIndex(bone).closest_points(p[used],certify=True);assert near['distance_m'].min()>=plan['actual_serialized_bone_clearance_metres'];quality=audit(mesh,require_closed=False,expected_components=1);assert quality['passed'];tri=mesh.triangles;geometric=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);geometric/=np.linalg.norm(geometric,axis=1)[:,None];cosine=np.einsum('fci,fi->fc',n[faces].astype(float)/32767,geometric);assert cosine.min()>0
            atomic_npz(OUT/(side+'-course.npz'),vertices=p,faces=faces,frame=np.array(frame),normal_words=n,uvs=u)
            sides.append(dict(side=side,geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,actual_bone_contacts=auditproof,minimum_actual_vertex_bone_clearance_m=float(near['distance_m'].min()),minimum_forward_corner_normal_cosine=float(cosine.min()),self_quality=quality))
        proof.update(passed=True,sides=sides,source_root_positions_and_normals_literal=True,mirror_transfer_exact=True,new_source_face_count=0,new_anatomical_caps=0)
    except Exception as error:
        import traceback
        proof['failure']=dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc())
    proof['elapsed_seconds']=time.monotonic()-start;atomic_json(OUT/'proof.json',proof);print(json.dumps(dict(passed=proof['passed'],elapsed_seconds=proof['elapsed_seconds'],sections=proof.get('actual_outer_calibre_sections'),failure=proof.get('failure'),sides=[dict(side=x['side'],minimum_bone_gap=x['minimum_actual_vertex_bone_clearance_m'],minimum_forward_normal=x['minimum_forward_corner_normal_cosine'])for x in proof.get('sides',[])])))

if __name__=='__main__':main()
