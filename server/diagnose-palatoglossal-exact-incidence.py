"""Retain the exact finite-triangle incidence proof for the rejected PG cut.

Read-only with respect to anatomy and producer/core code. This writes one new
diagnostic directory and never emits a qualifying or installable candidate.
"""
from __future__ import annotations
import argparse, importlib.util, json, sys, time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.algorithms.exact_source_intersections import contains,edge_plane,_sub,_dot,_cross
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance,subtract_occupied
from palatoglossal_incidence_fixture import rational_controls,common_parent_arrangement,oriented_area
W=ROOT/'.local-models/geometry-workbench/tongue-completion'

def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def qvec(v):return tuple(Q.from_float(float(x))for x in v)
def rational_string(v):return [str(x)for x in v]
def determinant(rows):return _dot(rows[0],_cross(rows[1],rows[2]))
def plane(tri):
    n=_cross(_sub(tri[1],tri[0]),_sub(tri[2],tri[0]));return n,_dot(n,tri[0])
def exact_solution(triangles):
    eq=[plane(t)for t in triangles];a=[r[0]for r in eq];b=[r[1]for r in eq];d=determinant(a)
    require(bool(d),'Three source planes are not independent')
    answer=[]
    for column in range(3):
        matrix=[tuple(b[i]if j==column else a[i][j]for j in range(3))for i in range(3)]
        answer.append(determinant(matrix)/d)
    return tuple(answer),d
def words(mesh):
    return mesh.attributes['normals'][mesh.faces].astype(np.int32)*mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),dtype='i1'))[mesh.faces][:,:,None]

def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use a new diagnostic directory')
    checkpoint=W/'palatoglossal-mucosal-entry-v2-owner-checkpoint';ev=read(checkpoint/'evidence.json');z=np.load(checkpoint/'failed-incidence-witness.npz')
    tri=z['triangles'];point=z['point'];qtri=[tuple(qvec(p)for p in t)for t in tri];hit,det=exact_solution(qtri)
    common=set(qtri[1])&set(qtri[2]);require(len(common)==2,'Expected the recorded exact shared cutter edge')
    a,b=sorted(common);edge_hit=edge_plane(a,b,qtri[0]);require(edge_hit==hit,'Edge-plane candidate differs from exact triple-plane solution')
    inside=[contains(hit,t)for t in qtri];require(inside==[False,True,True],'Recorded finite-source rejection changed')
    plan=read(W/'palatoglossal-mucosal-entry-v2/plan.json');mu=Mesh.load(plan['mucosa']);whole=Mesh.load(plan['whole']);palate=Mesh.load(plan['palate']);addon=read(plan['addon']);reader=AtlasReader()
    pg={p['side']:reader.read_part(addon,p,frame='model:local-ta98:assembled-metres')for p in addon['parts']}
    m,mi=make_solid(mu.vertices,mu.faces,properties=words(mu));w,wi=make_solid(whole.vertices,whole.faces);pal,pi=make_solid(palate.vertices,palate.faces)
    solids={};ids={mi:('tongue-mucosa',mu),wi:('whole-tongue',whole),pi:('native-palate',palate)}
    for side in ['left','right']:
        solid,pid=make_solid(pg[side].vertices,pg[side].faces,properties=words(pg[side]));solids[side]=solid;ids[pid]=('PG:'+side,pg[side])
    cutter=solids['left']+solids['right']+pal+(w-m);raw=cutter.to_mesh64();cp=np.asarray(raw.vert_properties)[:,:3];cf=np.asarray(raw.tri_verts);oid,pid=provenance(raw)
    saved=np.load(checkpoint/'owner-arrangement-and-cutter.npz')
    # Original IDs change run ordering. Match literal finite triangles, never
    # nearby points, so source recovery remains independent of row ordering.
    def triangle_key(t):return tuple(sorted(tuple(p)for p in t))
    new_faces={triangle_key(cp[face]):i for i,face in enumerate(cf)}
    old_tri=saved['composite_cutter_positions'][saved['composite_cutter_faces']]
    require(all(triangle_key(old_tri[i])in new_faces for i in [517,519]),'Literal cutter witness triangles absent from reproduction')
    parent_rows=[];source_arrays={}
    for old_fi in [517,519]:
        fi=new_faces[triangle_key(old_tri[old_fi])]
        name,source=ids[int(oid[fi])];parent=int(pid[fi]);st=source.vertices[source.faces[parent]];source_arrays['cutter_source_'+str(old_fi)]=st
        sq=tuple(qvec(p)for p in st);incidences=[contains(qvec(p),sq)for p in cp[cf[fi]]]
        parent_rows.append(dict(cutter_face=old_fi,reproduced_cutter_face=fi,original_source=name,source_face=parent,all_cutter_vertices_exactly_in_ancestor_triangle=all(incidences),vertex_exact_incidence=incidences,
            ancestor_plane_residual_float=[float(_dot(plane(sq)[0],_sub(qvec(p),sq[0])))for p in cp[cf[fi]]]))
    # The first-source cell is flattened from its actual named shared surfaces.
    meta=ev['surfaces'];surface={r['id']:saved[r['array']]for r in meta};left_uses=ev['cells']['left'];flat=[];flat_parent=[]
    for use in left_uses:
        fs=surface[use['surface']];flat.extend(fs if use['side']==1 else fs[:,[0,2,1]]);flat_parent.extend((use['surface'],i,use['side'])for i in range(len(fs)))
    sf=638;require(np.array_equal(saved['positions'][np.asarray(flat)[sf]],tri[0]),'First source witness differs from flattened cell')
    # Reproduce the already successful owner partition to recover its original
    # anatomical source face, rather than making coplanarity inferences by distance.
    hosts={s:Mesh.load(W/'palatoglossal-mucosal-entry-v2'/(s+'-uncut-expanded-arch.npz'))for s in ['left','right']}
    pos=np.vstack([hosts[s].vertices for s in ['left','right']]);surfaces={};rows={};cells={};parents={};uses={};offset=0
    for side in ['left','right']:
        h=hosts[side];sid='arch:'+side;surfaces[sid]=h.faces+offset;rows[sid]=dict(id=sid,kind='outer',array=sid);cells[side]=surfaces[sid];parents[side]=[(sid,i,1)for i in range(len(h.faces))];uses[side]=[dict(surface=sid,side=1)];offset+=len(h.vertices)
    partition=subtract_occupied(pos,surfaces,rows,cells,parents,uses,mu.vertices,mu.faces,occupied_face_owners=['tongue-mucosa']*len(mu.faces),partition_free_outer=True,exact_source_replay=True,reconcile_shared_subdivisions=True)
    require(np.array_equal(partition['positions'],saved['positions']),'Owner arrangement reproduction differs')
    sid,index,sign=flat_parent[sf];chart=partition['charts'][sid];parent=int(chart['parentFaces'][index]);source_arrays['first_source_ancestor']=mu.vertices[mu.faces[parent]]
    first_parent=dict(surface=sid,surface_local_face=index,side=sign,original_source='tongue-mucosa',source_face=parent)
    original_parent_planes_equal=all(r['original_source']=='tongue-mucosa'and r['source_face']==parent for r in parent_rows)
    require(original_parent_planes_equal,'The local common-parent fixture cannot apply across distinct original source faces')
    hs,hi=make_solid(hosts['left'].vertices,hosts['left'].faces);ms,mi2=make_solid(mu.vertices,mu.faces)
    owner_raw=(hs-ms).to_mesh64();owner_q,owner_ev=rational_controls(owner_raw,{hi:hosts['left'],mi2:mu})
    exact_owner_by_float={tuple(float(x)for x in q):q for q in owner_q.values()}
    require(all(tuple(p)in exact_owner_by_float for p in tri[0]),'Literal rounded owner controls absent from original-parent replay')
    first_exact=tuple(exact_owner_by_float[tuple(p)]for p in tri[0])
    cut_selected={tuple(p)for i in [517,519]for p in old_tri[i]};cut_q,cut_ev=rational_controls(raw,{i:row[1]for i,row in ids.items()},selected_points=cut_selected)
    corrected=[first_exact]+[tuple(cut_q[tuple(p)]for p in old_tri[i])for i in [517,519]]
    ancestor=tuple(qvec(p)for p in source_arrays['first_source_ancestor']);require(all(contains(p,ancestor)for t in corrected for p in t),'Reconstructed source corners are outside their true source facet')
    axes,common_poly=common_parent_arrangement(corrected,[('tongue-mucosa',parent)]*3,ancestor)
    require(bool(common_poly),'The exact ancestor arrangement has no shared finite area/edge')
    original_n,original_d=plane(ancestor);drop=next(i for i in range(3)if i not in axes)
    def lift(p):
        q=[Q(0)]*3
        for i,x in zip(axes,p):q[i]=x
        q[drop]=(original_d-sum(original_n[i]*q[i]for i in axes))/original_n[drop]
        return tuple(q)
    common_vertices=[lift(p)for p in common_poly]
    require(all(contains(p,t)for p in common_vertices for t in corrected),'Exact common arrangement is not incident to all recovered source triangles')
    require(len(common_vertices)==2,'This fixture expects one exact finite common segment')
    ca,cb=common_vertices;direction=_sub(cb,ca);target=qvec(point);parameter=_dot(_sub(target,ca),direction)/_dot(direction,direction);parameter=max(Q(0),min(Q(1),parameter))
    chosen=tuple(ca[i]+parameter*direction[i]for i in range(3));encoded=np.array([float(x)for x in chosen]);error=np.array([float(abs(Q.from_float(float(x))-q))for x,q in zip(encoded,chosen)])
    require(all(contains(chosen,t)for t in corrected),'Chosen common source control is outside reunited fragments')
    bound=np.spacing(np.abs(encoded))/2;require(np.all(error<=bound),'Single common rounded encoding exceeds half-ULP coordinate bound')
    correction_fixture=dict(original_parent_source='tongue-mucosa',original_parent_face=parent,exact_original_parent_vertices=[rational_string(p)for p in ancestor],
        corrected_fragment_triangles=[[rational_string(p)for p in t]for t in corrected],all_corrected_fragment_corners_in_literal_original_parent=True,
        common_exact_polygon=[rational_string(p)for p in common_vertices],common_polygon_signed_chart_area=str(oriented_area(common_poly)),
        maximum_original_ancestry_replay_m=max([r['replay_m']for r in owner_ev.values()]+[r['replay_m']for r in cut_ev.values()],default=0.),
        selected_common_control=rational_string(chosen),selected_common_control_float=encoded.tolist(),selected_control_distance_from_boolean_point_m=float(np.linalg.norm(encoded-point)),
        encoded_coordinate_error_m=error.tolist(),analytical_encoding_bound_per_coordinate_m=bound.tolist(),single_common_encoding=True,
        original_distinct_planes_collapsed=False,distance_weld=False,source_points_moved=False,
        method='Recover exact rational controls from original finite triangle ancestry, then intersect fragment polygons in the single original source-face chart; no fragment-plane proximity test')
    output.mkdir(parents=True);atomic_npz(output/'witness-and-ancestors.npz',point=point,triangles=tri,exact_hit_float=np.array([float(x)for x in hit]),**source_arrays)
    report=dict(schema='palatoglossal-exact-incidence-diagnostic/v1',qualifying_candidate=False,installed=False,
        failed_function='geometry_tools.algorithms.exact_source_intersections.replay_boolean_vertices',failed_predicate='all(contains(q,tri) for group in source for tri in group)',
        maximum_allowed_replay_m=2e-10,point=point.tolist(),exact_triple_plane_hit=rational_string(hit),exact_triple_plane_hit_float=[float(x)for x in hit],
        exact_determinant=str(det),determinant_float=float(det),hit_finite_triangle_containment=inside,shared_edge_plane_hit_equals_triple_plane_hit=True,
        triple_plane_hit_distance_from_boolean_point_m=float(np.linalg.norm(np.array([float(x)for x in hit])-point)),
        exact_hit_dyadic_coordinates=[not(x.denominator&(x.denominator-1))for x in hit],
        proposed_float_point_exactly_contained=[contains(qvec(point),t)for t in qtri],
        first_source_parent=first_parent,cutter_source_parents=parent_rows,all_three_facets_have_same_original_source_face=original_parent_planes_equal,
        reunited_original_parent_fixture=correction_fixture,
        missing_triple_plane_candidate=False,correct_rejection=True,
        correction='Reconstruct one canonical arrangement on original parent source planes before further clipping. Keep exact rational intersection/provenance until one common rounded encoding. Do not use rounded fragment triangles as independent authoritative source planes.',
        source_bindings={**plan['bindings'],str(checkpoint/'evidence.json'):digest_file(checkpoint/'evidence.json'),str(checkpoint/'failed-incidence-witness.npz'):digest_file(checkpoint/'failed-incidence-witness.npz'),str(checkpoint/'owner-arrangement-and-cutter.npz'):digest_file(checkpoint/'owner-arrangement-and-cutter.npz'),str(Path(__file__).resolve()):digest_file(__file__),str(Path(__file__).with_name('palatoglossal_incidence_fixture.py')):digest_file(Path(__file__).with_name('palatoglossal_incidence_fixture.py'))},
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'evidence.json',report);print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
