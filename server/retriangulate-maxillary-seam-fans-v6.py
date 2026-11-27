"""One deterministic source-plane diagonal flip in each diagnosed CSG quad.

Positions and vertex fields stay literal; only four two-face fan incidences may
change. A projected-domain/envelope argument bounds surface Hausdorff error.
"""
import argparse,hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops
H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,default=lambda a:a.tolist()if hasattr(a,'tolist')else str(a))+'\n')
def run(work,out):
 start=time.monotonic();assert not out.exists();out.mkdir(parents=True);raw=work/'candidate-v4';v5=work/'quantized-junction-v5';source=json.loads((raw/'proof.json').read_text());failure=json.loads((v5/'compiled-cohort-audit.json').read_text());cleanup=json.loads((v5/'storage-cleanup-proof.json').read_text());assert source['numerical_gate_pass']and failure['failures']==['self_intersection_check']
 mesh=Mesh.load(v5/'compiled-cohort.npz');a=np.load(raw/'ancestry.npz',allow_pickle=False);keep=np.asarray(cleanup['retained_candidate_face_ids'],int);owner=a['source_indices'][keep];parent=a['source_face_ids'][keep];kinds=a['face_kinds'][keep];faces=mesh.faces.copy();original_faces=faces.copy();gv,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True);gf=inverse[faces];groups=set()
 for pair in failure['self_intersections']['examples']:
  for fi in pair:groups.add((int(owner[fi]),int(parent[fi])))
 assert len(groups)==4,'Only the four independently diagnosed source-parent quads are authorized'
 changes=[]
 for own,par in sorted(groups):
  ids=np.flatnonzero((owner==own)&(parent==par));assert len(ids)==2 and np.all(kinds[ids]=='original')
  loops=boundary_loops(gf[ids]);assert len(loops)==1 and len(loops[0])==4
  loop=loops[0];point_ids=np.unique(gf[ids]);assert set(point_ids)==set(loop)
  shared=set(gf[ids[0]])&set(gf[ids[1]]);assert len(shared)==2
  offset=1 if shared=={int(loop[0]),int(loop[2])}else 0
  assert shared=={int(loop[0]),int(loop[2])}or shared=={int(loop[1]),int(loop[3])}
  ordered=np.roll(loop,-offset);new=np.asarray([[ordered[0],ordered[1],ordered[2]],[ordered[0],ordered[2],ordered[3]]])
  chart=a['parent_triangles'][keep[ids[0]]];assert np.array_equal(chart,a['parent_triangles'][keep[ids[1]]])
  normal=np.cross(chart[1]-chart[0],chart[2]-chart[0]);normal/=np.linalg.norm(normal);u=chart[1]-chart[0];u/=np.linalg.norm(u);v=np.cross(normal,u)
  points=gv[loop];xy=np.column_stack([(points-chart[0])@u,(points-chart[0])@v])
  oldtri=gv[gf[ids]];newtri=gv[new]
  oldareas=np.einsum('ij,j->i',np.cross(oldtri[:,1]-oldtri[:,0],oldtri[:,2]-oldtri[:,0]),normal)/2
  newareas=np.einsum('ij,j->i',np.cross(newtri[:,1]-newtri[:,0],newtri[:,2]-newtri[:,0]),normal)/2
  assert np.all(newareas>0)
  polygon_area=np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1))/2
  assert abs(oldareas.sum()-polygon_area)<1e-17 and abs(newareas.sum()-polygon_area)<1e-17
  delta=float(np.abs((points-chart[0])@normal).max());strip_width=0.;negative=np.flatnonzero(oldareas<0)
  assert len(negative)<=1 and np.any(oldareas>0)
  if len(negative):
   # One quantized fan triangle folded across its almost-collinear edge.
   # Its replacement is the positive quad inside the positive old triangle.
   # The removed projected strip has width at most its minimum altitude.
   big=oldtri[np.argmax(oldareas)];bxy=np.column_stack([(big-chart[0])@u,(big-chart[0])@v]);matrix=np.column_stack([bxy[1]-bxy[0],bxy[2]-bxy[0]]);weights=np.linalg.solve(matrix,(xy-bxy[0]).T).T;bary=np.column_stack([1-weights.sum(1),weights]);assert bary.min()>=-1e-10
   thin=oldtri[negative[0]];edge_lengths=np.linalg.norm(thin-np.roll(thin,-1,axis=0),axis=1);strip_width=float(2*abs(oldareas[negative[0]])/edge_lengths.max())
  else:assert np.all(oldareas>0)
  bound=float(np.hypot(strip_width,2*delta));assert bound<=1e-8,'Source-plane graph/strip envelope exceeds the ten-nanometre budget'
  chart_indices={}
  for vi,gi in zip(faces[ids].ravel(),gf[ids].ravel()):
   if int(gi)in chart_indices:
    old=chart_indices[int(gi)];assert np.array_equal(mesh.attributes['normals'][old],mesh.attributes['normals'][vi])and np.array_equal(mesh.attributes['uvs'][old],mesh.attributes['uvs'][vi])
   else:chart_indices[int(gi)]=int(vi)
  faces[ids]=np.vectorize(chart_indices.__getitem__)(new)
  changes.append(dict(source_owner=own,source_parent_face=par,changed_v5_face_ids=ids.tolist(),original_v4_face_ids=keep[ids].tolist(),before_corner_indices=original_faces[ids].tolist(),after_corner_indices=faces[ids].tolist(),all_four_boundary_controls_literal=True,source_plane_max_vertex_residual_m=delta,maximum_surface_Hausdorff_bound_m=bound,projected_polygon_area_m2=float(polygon_area),old_projected_inverted_strip_width_m=strip_width,positive_replacement_projected_triangles=True,proof='All triangle graph heights lie in one ancestry-plane slab of halfwidth delta. Positive old/new fans cover the same projected polygon. Where the quantized old fan has one inverted thin triangle, the new positive quad lies inside the old positive triangle; the removed projected strip width is bounded by the thin triangle minimum altitude h. Bidirectional surface Hausdorff distance is bounded by hypot(h,2delta).'))
 candidate=Mesh(mesh.vertices,faces,mesh.frame,{k:v.copy()for k,v in mesh.attributes.items()},dict(source_master_sha256=source['source_master_sha256'],storage_diagonal_fans_only=True));assert np.array_equal(candidate.vertices,mesh.vertices)and all(np.array_equal(candidate.attributes[k],mesh.attributes[k])for k in mesh.attributes)
 candidate.save(out/'cohort.npz');np.savez_compressed(out/'ownership.npz',source_indices=owner,original_v4_face_ids=keep);check=audit(candidate,require_closed=False,expected_components=2,max_pairs=600000,max_tests=3000000);write(out/'actual-Float32-audit.json',check)
 gv,inv=np.unique(candidate.vertices,axis=0,return_inverse=True);afterloops={frozenset(map(tuple,gv[x]))for x in boundary_loops(inv[faces])};beforeloops={frozenset(map(tuple,gv[x]))for x in boundary_loops(inv[original_faces])};assert len(beforeloops)==len(afterloops)==22 and beforeloops==afterloops
 result=dict(schema='maxillary-storage-source-plane-fan-retriangulation/v6',complete=True,passed=check['passed'],source_master_sha256=source['source_master_sha256'],producer_v4_sha256=H(raw/'cohort.npz'),v5_input_sha256=H(v5/'compiled-cohort.npz'),cohort_sha256=H(out/'cohort.npz'),all_positions_normal_UV_arrays_literal_from_v5=True,changed_fans=changes,changed_face_count=sum(len(c['changed_v5_face_ids'])for c in changes),maximum_surface_Hausdorff_bound_m=max(c['maximum_surface_Hausdorff_bound_m']for c in changes),all22external_profile_fields_literal=True,original_source_prefix_protected=True,one_deterministic_flip_per_diagnosed_quad=True,threshold_search=False,actual_geometry_failures=check['failures'],installed=False,elapsed_seconds=time.monotonic()-start)
 write(out/'proof.json',result);print(json.dumps({k:v for k,v in result.items()if k!='changed_fans'}));assert check['passed'],check['failures']
 spec=importlib.util.spec_from_file_location('pack_maxillary',ROOT/'server/pack-maxillary-native-junction-v4.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.run(raw,out/'packed-review',prepared_storage_packet=out)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.work.resolve(),a.output.resolve())
