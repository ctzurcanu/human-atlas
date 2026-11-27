"""One explicit adjacent-unjoined, Float32-safe GG attachment variant.

Only new cortical contact/cut controls move, by a nominal 0.5 micrometre
tissue-side offset. All source bones remain literal. The finite precision
clearance is recorded; this does not claim a shared histological enthesis.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.distance import point_triangle
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import diagnostics

BASE=ROOT/'.local-models/geometry-workbench/tongue-completion'
SOURCE=BASE/'source-v1/meshes';PREVIOUS=BASE/'genioglossus-v2-stage-b';OUT=BASE/'genioglossus-v2-encoding-safe'
START=time.monotonic();DEADLINE=START+35
ORIGIN=np.array([0.,1.625,0.]);NOMINAL=.5e-6;TOLERANCE=1e-6


def check():
    if time.monotonic()>DEADLINE:raise TimeoutError('GG encoding variant exceeded its 35 second native deadline')


def minimum_offset(normals,signed):
    # Dimensionless active-set QP, at most three independent constraints in
    # three dimensions. The displacement must progress outside every actual
    # nearby cortical facet, not only along one averaged shading normal.
    b=(NOMINAL-signed)/NOMINAL;solutions=[]
    for k in range(1,min(3,len(normals))+1):
        for ids in combinations(range(len(normals)),k):
            a=normals[list(ids)];gram=a@a.T
            if np.linalg.cond(gram)>1e10:continue
            weights=np.linalg.solve(gram,b[list(ids)])
            if np.any(weights<-1e-10):continue
            u=a.T@weights
            if np.all(normals@u>=b-1e-10):solutions.append(u)
    if not solutions:raise ValueError('No finite tissue-side incident-normal offset cone')
    value=min(solutions,key=lambda x:x@x)*NOMINAL
    if np.linalg.norm(value)>.9e-6:raise ValueError('Cortical offset exceeds the explicit submicrometre budget')
    return value


def main():
    OUT.mkdir(exist_ok=False)
    bones=[Mesh.load(SOURCE/f'native-{i:02d}.npz')for i in [5,13,14,15,16,17]]
    bv=[];bf=[];owners=[];offset=0
    for bone in bones:
        bv.append(bone.vertices);bf.append(bone.faces+offset);owners.extend([bone.source['record']['id']]*len(bone.faces));offset+=len(bone.vertices)
    fixed=Mesh(np.vstack(bv),np.vstack(bf),bones[0].frame);tree=BVH(fixed)
    tri=fixed.triangles;area=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normals=area/np.linalg.norm(area,axis=1)[:,None]
    report=[]
    for side,index in [('left',1),('right',2)]:
        check();parent=Mesh.load(PREVIOUS/f'genioglossus-{side}.npz');original=Mesh.load(SOURCE/f'native-{index:02d}.npz')
        points,alias=np.unique(parent.vertices,axis=0,return_inverse=True)
        near,evidence=tree.closest_points(points,backend='open3d')
        controls=np.flatnonzero(near['distance_m']<1e-9);proposed=points.copy();moves=[]
        for control in controls:
            check();p=points[control];ids=tree.query_bounds(p,p,distance_m=1e-9)
            d,_,_=point_triangle(np.repeat(p[None],len(ids),axis=0),tri[ids]);ids=ids[d<1e-18]
            if not len(ids):raise ValueError('New contact control has no literal cortical face witness')
            signed=np.einsum('ij,ij->i',p-tri[ids,0],normals[ids])
            delta=minimum_offset(normals[ids],signed);proposed[control]+=delta
            moves.append(dict(compact_control=int(control),cortical_source_face_ids=ids.tolist(),
                cortical_owners=[owners[i]for i in ids],displacement_m=float(np.linalg.norm(delta)),
                direction=delta.tolist()))
        encoded=(proposed-ORIGIN).astype('<f4');stored=encoded.astype(float)+ORIGIN
        vertex_values=stored[alias]
        candidate=Mesh(vertex_values,parent.faces,parent.frame,parent.attributes,source={
            'record':parent.source['record'],'parent_geometry_sha256':parent.geometry_sha,
            'geometryOrigin':ORIGIN.tolist(),'stored_position_dtype':'Float32',
            'geometric_attachment':'adjacent-unjoined','geometric_attachment_tolerance_m':TOLERANCE,
            'nominal_tissue_side_clearance_m':NOMINAL,'literal_shared_coordinates':False,
            'histological_enthesis_completion':False,'anatomical_acceptance':False})
        path=OUT/f'genioglossus-{side}.npz';candidate.save(path)
        stored_path=OUT/f'genioglossus-{side}-stored-float32.npz'
        np.savez_compressed(stored_path,positions=encoded[alias],geometryOrigin=ORIGIN,indices=candidate.faces.astype('<u4'),
            normals=candidate.attributes['normals'],uvs=candidate.attributes['uvs'])
        tests=[]
        for bone in bones:
            check();contacts,summary=pair_contacts(candidate,bone,distance_m=0.,max_pairs=25000,max_tests=200000)
            np.savez_compressed(OUT/f'{side}-{bone.source["record"]["id"].replace(":","_")}-contacts.npz',**contacts)
            tests.append(dict(bone_part_id=bone.source['record']['id'],**summary))
        diag=diagnostics(candidate,self_intersections=True);topo=topology(candidate)
        key=lambda t:min(tuple(map(tuple,np.roll(t,i,axis=0)))for i in range(3))
        emitted={key(t)for t in candidate.triangles};upper=np.flatnonzero(np.all(original.triangles[:,:,1]>=1.630,axis=1))
        preserved=sum(key(original.triangles[i])in emitted for i in upper)
        actual,_=tree.closest_points(stored[controls],backend='open3d')
        position_only_rounding=float(np.linalg.norm(stored-proposed,axis=1).max())
        row=dict(side=side,mesh=str(path),stored_float32_bank=str(stored_path),
            mesh_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),stored_bank_sha256=hashlib.sha256(stored_path.read_bytes()).hexdigest(),
            geometric_attachment='adjacent-unjoined',literal_shared_coordinates=False,histological_enthesis_completion=False,
            nominal_offset_m=NOMINAL,maximum_allowed_attachment_gap_m=TOLERANCE,
            contact_controls=len(controls),contact_control_moves=moves,
            actual_stored_contact_gap_min_m=float(actual['distance_m'].min()),actual_stored_contact_gap_max_m=float(actual['distance_m'].max()),
            maximum_position_encoding_error_m=position_only_rounding,
            unmodified_controls_only_encoded=True,upper_source_triangles=len(upper),upper_source_triangles_retained_literal=preserved,
            topology=topo,open3d=diag,bone_contacts=tests,
            numerical_geometry_passed=bool(topo['closed_oriented_manifold']and diag['self_intersection_pairs']==0 and
                all(t['crossing_face_pairs']==0 for t in tests)and preserved==len(upper)and actual['distance_m'].max()<=TOLERANCE),
            source_bones_modified=False,installed=False,anatomical_acceptance=False)
        report.append(row)
    result=dict(schema='native-gg-encoding-safe-adjacency/v1',elapsed_seconds=time.monotonic()-START,
        native_deadline_seconds=35,source_bones_modified=False,predicates_relaxed=False,candidates=report,
        installed=False,anatomical_acceptance=False)
    (OUT/'candidate.json').write_text(json.dumps(result,indent=2)+'\n')
    spec=json.loads((PREVIOUS/'review-focused-spec.json').read_text())
    for row in spec['meshes']:
        for side in ['left','right']:
            if row['mesh'].endswith(f'genioglossus-{side}.npz'):row['mesh']=str(OUT/f'genioglossus-{side}.npz');row['name']=f'Float32 adjacent-unjoined GG attachment ({side})'
    (OUT/'review-spec.json').write_text(json.dumps(spec,indent=2)+'\n')
    print(json.dumps(dict(report=str(OUT/'candidate.json'),elapsed_seconds=result['elapsed_seconds'],
        candidates=[{k:r[k]for k in ['side','contact_controls','actual_stored_contact_gap_min_m','actual_stored_contact_gap_max_m','maximum_position_encoding_error_m','upper_source_triangles_retained_literal','numerical_geometry_passed']}for r in report]),indent=2))


if __name__=='__main__':main()
