"""Retain the stopped direct-arch single-parent control and its complete star.

This audit reproduces the same Boolean once solely to save the failed witness.
It emits no candidate and never retries a different construction or tolerance.
"""
import argparse,importlib.util,json,sys,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
from geometry_tools.algorithms.exact_source_intersections import contains,_sub,_dot,_cross
from palatoglossal_incidence_fixture import qvec
from palatoglossal_original_provenance import interior_steiner_control
import manifold3d as md
_spec=importlib.util.spec_from_file_location('pg_direct_arch',Path(__file__).with_name('construct-palatoglossal-direct-arch.py'));producer=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(producer)
W=ROOT/'.local-models/geometry-workbench/tongue-completion';FRAME='model:local-ta98:assembled-metres'
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def main(output,source_stage,point,full_audit):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use a new witness directory');plan_path=W/source_stage/'plan.json';p=read(plan_path)
    for path,sha in p['bindings'].items():require(digest_file(path)==sha,'Source changed: '+path)
    meshes={name:Mesh.load(path)for name,path in p['sources'].items()};addon=read(p['addon']);reader=AtlasReader()
    for row in addon['parts']:meshes['PG_'+row['side']]=reader.read_part(addon,row,frame=FRAME)
    solids={};source={};family={};ids={}
    for name,m in meshes.items():
        solid,oid=make_solid(m.vertices,m.faces,properties=producer.words(m));solids[name]=solid;source[oid]=m;family[oid]=name;ids[name]=oid
    nm=read(p['direct_classifier']);negative=Mesh.load(nm['negative_mesh']);basal_id=md.Manifold.reserve_ids(1);actual_m=np.r_[nm['original_inner_source_faces'],nm['original_rim_source_faces']];core=producer.direct_classifier(negative,ids['mucosa'],basal_id,actual_m)
    core.to_mesh64();pg_union=solids['PG_left']+solids['PG_right'];pg_union.to_mesh64();pg_pal=pg_union+solids['palate'];pg_pal.to_mesh64();cutter=pg_pal+core;cutter.to_mesh64()
    left=solids['arch_left']-solids['mucosa'];left.to_mesh64();right=solids['arch_right']-solids['mucosa'];right.to_mesh64();result=left-cutter;raw=result.to_mesh64();v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts);oid,pid=provenance(raw)
    point=np.asarray(point);controls=np.flatnonzero(np.all(v==point,axis=1));require(len(controls)>0,'Literal failed control differs from same Boolean reproduction')
    face_ids=np.flatnonzero(np.any(np.isin(f,controls),axis=1));families=set((int(oid[i]),int(pid[i]))for i in face_ids);require(len(families)==1,'Failed control has additional original parents')
    original,parent=next(iter(families));require(family[original]=='mucosa'and parent==405,'Failed original source changed');triangle=source[original].triangles[parent];qt=tuple(qvec(p)for p in triangle);a,b,c=qt;n=_cross(_sub(b,a),_sub(c,a));qp=qvec(point);scale=_dot(n,_sub(qp,a))/_dot(n,n);on_plane=tuple(qp[i]-scale*n[i]for i in range(3));require(contains(on_plane,qt),'Single-parent control projects outside its finite original parent')
    # Literal geometry-coordinate vertex links avoid chart-attribute duplicates.
    neighbors=[]
    for face_id in face_ids:
        face=f[face_id];index=next(i for i in range(3)if np.array_equal(v[face[i]],point));neighbors.append((tuple(v[face[(index+1)%3]]),tuple(v[face[(index+2)%3]])))
    degree=Counter(p for edge in neighbors for p in edge);closed_link=bool(degree)and all(d==2 for d in degree.values())
    same_parent=True;round_trip=np.array([float(x)for x in on_plane]);distance=float(np.linalg.norm(round_trip-point));require(distance<=2e-10,'Single-plane analytical point exceeds unchanged replay bound')
    try:
        actual_rule=interior_steiner_control(point,v[f[face_ids]],[(original,parent)],triangle);rule_result=dict(passed=True,control=[str(x)for x in actual_rule])
    except ValueError as error:rule_result=dict(passed=False,guard=str(error))
    output.mkdir(parents=True);atomic_npz(output/'actual-star.npz',point=point,source_triangle=triangle,star_triangles=v[f[face_ids]],output_face_ids=face_ids,analytical_single_plane_control=round_trip)
    full_result=None
    if full_audit:
        actual=Mesh(v.copy(),f.copy(),FRAME);actual.save(output/'rejected-raw-left-arch.npz');full_result=audit(actual,require_closed=True,max_pairs=250000,max_tests=2000000);atomic_json(output/'rejected-raw-left-audit.json',full_result)
    report=dict(schema='PG-single-parent-coplanar-star-witness/v1',qualifying_candidate=False,installed=False,source_family='mucosa',source_face=parent,raw_output_control=point.tolist(),incident_faces=face_ids.tolist(),incident_face_count=len(face_ids),
        every_incident_face_has_identical_original_parent=True,literal_geometry_link_is_closed=closed_link,link_neighbor_degrees=list(degree.values()),
        strict_identical_encoding_steiner_rule=rule_result,
        failed_raw_geometry_audit=full_result,
        original_finite_parent_incidence_established=contains(on_plane,qt),analytical_control=[str(x)for x in on_plane],analytical_control_float=round_trip.tolist(),distance_to_reported_float_m=distance,
        diagnosis='Same-parent interior control qualifies only if its full oriented link is one connected cycle; the retained strict rule result takes precedence over the weaker degree-two census.',
        possible_producer_rule='Only for a complete closed same-original-face star, reconstruct an exact rational interior barycentric/plane control; never apply to boundary stars or multiple distinct original planes. Emit once with analytical encoding bound.',
        producer_rule_applied=False,raw_failed_stage=str(W/source_stage/'failed-stage.json'),source_bindings={str(plan_path):digest_file(plan_path),str(Path(__file__).resolve()):digest_file(__file__)},elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'evidence.json',report);print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--source-stage',default='palatoglossal-direct-arch-v1');p.add_argument('--point',type=float,nargs=3,default=[.015807216472070562,1.6338496977065908,.04814157652421331]);p.add_argument('--full-audit',action='store_true');a=p.parse_args();main(a.output.resolve(),a.source_stage,a.point,a.full_audit)
