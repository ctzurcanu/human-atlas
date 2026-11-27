"""Retain failed source incidence and assess a literal-sheet inner classifier.

This is a source-bound feasibility report only. It emits no anatomy/candidate,
never changes the installed model and does not continue the failed Boolean.
"""
import argparse,importlib.util,json,sys,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,boundary_signature
from geometry_tools.algorithms.exact_source_intersections import contains,_sub,_dot,_cross
_spec=importlib.util.spec_from_file_location('pg_exact_diagnostic',Path(__file__).with_name('diagnose-palatoglossal-exact-incidence.py'));_diagnostic=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_diagnostic)
exact_solution,qvec,rational_string=_diagnostic.exact_solution,_diagnostic.qvec,_diagnostic.rational_string
W=ROOT/'.local-models/geometry-workbench/tongue-completion';O=np.array([0.,1.625,0.]);FRAME='model:local-ta98:assembled-metres'
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def words(mesh):
    if 'normals'not in mesh.attributes:
        t=mesh.triangles;n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];return np.repeat(np.rint(n*32767).astype('<i2')[:,None,:],3,axis=1)
    return (mesh.attributes['normals'][mesh.faces].astype(np.int32)*mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),dtype='i1'))[mesh.faces][:,:,None]).astype('<i2')
def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use a new read-only diagnostic directory')
    plan_path=W/'palatoglossal-original-parent-v2/plan.json';plan=read(plan_path);meshes={name:Mesh.load(path)for name,path in plan['source_files'].items()};addon=read(plan['addon']);reader=AtlasReader()
    for row in addon['parts']:meshes['PG_'+row['side']]=reader.read_part(addon,row,frame=FRAME)
    solids={};source={};family={}
    for name,m in meshes.items():
        solid,oid=make_solid(m.vertices,m.faces,properties=words(m));solids[name]=solid;source[oid]=m;family[oid]=name
    raw=(solids['whole']-solids['mucosa']).to_mesh64();v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts);oids,pids=provenance(raw)
    point=np.array([-.010896470184692813,1.624745261768657,.05745634185342098]);ids=np.flatnonzero(np.all(v==point,axis=1));require(len(ids)>0,'Retained failed control missing from literal reproduction')
    parents=sorted(set((int(o),int(p))for face,o,p in zip(f,oids,pids)if set(map(int,face))&set(map(int,ids))))
    triangles=[source[o].triangles[p]for o,p in parents];qt=[tuple(qvec(x)for x in t)for t in triangles];require(len(qt)==3,'Expected retained three-source-plane witness')
    hit,det=exact_solution(qt);failed=dict(original_parents=[dict(source=family[o],face=p)for o,p in parents],point=point.tolist(),
        determinant=str(det),determinant_float=float(det),unique_exact_plane_intersection=rational_string(hit),exact_hit_finite_triangle_containment=[contains(hit,t)for t in qt],
        exact_hit_distance_from_boolean_point_m=float(np.linalg.norm(np.array([float(x)for x in hit])-point)),original_ids_available=True,
        missing_triple_plane_candidate=False,producer_rejected=True)
    shared=W/'shared-coverings-exposed-apo-sl-v3';construction_path=shared/'construction.json';construction=read(construction_path);exposed_path=Path(construction['exposedSource']);exposed=read(exposed_path)
    fields_path=exposed_path.with_name('fields.npz');excluded_path=exposed_path.with_name('excluded-basal-carrier.npz');fields=np.load(fields_path);outer=fields['fields'][0];carrier_faces=fields['faces'];excluded=Mesh.load(excluded_path)
    require(np.array_equal(outer,excluded.vertices),'Excluded basal carrier does not retain literal clipped carrier control identities')
    mu=meshes['mucosa'];count=exposed['carrierTriangles'];require(count==348,'Actual mucosal carrier count changed')
    actual_outer=(outer-O).astype('<f4').astype(float)+O
    require(np.array_equal(actual_outer[carrier_faces],mu.triangles[:count]),'Actual outer chart differs from its declared clipped controls/encoding')
    shared_original={r['mucosaOriginalFace']for r in construction['canonicalSourceFaceMapping']};require(len(shared_original)==62,'Expected literal62facet APO chart')
    # The producer explicitly emits: outer; remaining inner; reflection rim;
    # then reciprocal APO chart. Preserve that source identity, never infer it
    # from normal direction or distances.
    inner_original=[i for i in range(count,2*count)if i not in shared_original];rim_original=list(range(2*count,len(mu.faces)))
    free_original=[i for i in range(count,len(mu.faces))if i not in shared_original];require(len(free_original)==414,'Actual finite producer free face ordering changed')
    original_to_actual={old:count+i for i,old in enumerate(free_original)}
    inner_actual=[original_to_actual[i]for i in inner_original]+list(range(count+len(free_original),len(mu.faces)))
    rim_actual=[original_to_actual[i]for i in rim_original]
    require(len(inner_actual)==348 and len(rim_actual)==128,'Inner sheet/reflection rim identity is incomplete')
    triangles_for_closure=np.concatenate([mu.triangles[inner_actual][:,[0,2,1]],mu.triangles[rim_actual][:,[0,2,1]],actual_outer[excluded.faces]])
    # This exact edge census assesses closure feasibility without emitting a
    # reconstructed solid. All corners remain the literal source-sheet/rim
    # or the already declared outer-carrier encoding of original control IDs.
    lookup={};pos=[];faces=[]
    for tri in triangles_for_closure:
        face=[]
        for xyz in tri:
            key=tuple(xyz)
            if key not in lookup:lookup[key]=len(pos);pos.append(xyz)
            face.append(lookup[key])
        faces.append(face)
    boundary=boundary_signature(np.asarray(faces));edges=Counter(tuple(sorted((face[i],face[(i+1)%3])))for face in faces for i in range(3))
    possible=not boundary and all(n==2 for n in edges.values())
    feasibility=dict(feasible_literal_closure=possible,actual_inner_sheet_faces=inner_actual,actual_reflection_rim_faces=rim_actual,
        inner_sheet_count=len(inner_actual),literal_reflection_rim_count=len(rim_actual),excluded_basal_carrier_count=len(excluded.faces),
        planned_orientations=dict(inner_sheet='Opposite current mucosal cell',reflection_rim='Opposite current mucosal cell',basal_classifier_patch='Original whole exterior direction'),
        uncancelled_directed_edges=sum(boundary.values()),non_two_sided_edges=sum(n!=2 for n in edges.values()),literal_outer_identity_verified=True,
        original_source_face_normals_and_coordinates_changed=False,whole_minus_mucosa_boolean_required=False,inner_clipping_omitted=False,
        cap_semantics='Excluded basal patch closes a computational negative volume only. It must be proved outside every exposed arch/mucosal host before any final cut; never emitted as anatomy.',
        pending=['Audit closed classifier topology/self-crossings and tissue occupied volume','Prove classifier-only basal faces outside all target host cells','Retain literal inner-sheet/reflection/APO chart source controls and normal words','Use original-parent arrangement for PG/palate clipping and final reciprocal owner charts'])
    output.mkdir(parents=True);atomic_npz(output/'failed-whole-minus-M-witness.npz',point=point,triangles=np.asarray(triangles));atomic_json(output/'evidence.json',dict(schema='PG-inner-negative-source-feasibility/v1',failed_raw_route=failed,direct_inner_domain=feasibility,
        source_bindings={str(p):digest_file(p)for p in [plan_path,construction_path,exposed_path,fields_path,excluded_path,Path(__file__).resolve(),Path(plan['source_files']['mucosa']),Path(plan['source_files']['whole'])]},
        emitted_anatomical_meshes=0,qualifying_candidate=False,installed=False,elapsed_seconds=time.monotonic()-began));print(json.dumps(dict(failed=failed,feasibility=feasibility,output=str(output))))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
