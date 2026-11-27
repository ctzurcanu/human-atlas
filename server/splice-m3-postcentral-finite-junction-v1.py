"""Literal source-face replacement proof, still uninstalled and outer-only away from root."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank,barycentric
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';C=P/'finite-junction-educational-v1';OUT=C/'literal-splice-v1';FRAME='model:local-ta98:assembled-metres'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(P/'source-port-charts.npz');g=json.loads((C/'generation-complete.json').read_text());q=np.load(C/'outer-union.npz');sources={r:Mesh.load(C/(r+'-source-untouched.npz'))for r in ['parent','branch']};prefix=np.concatenate([sources['parent'].vertices,sources['branch'].vertices]);bank=CommonPositionBank(prefix,tolerance=2e-12);mapped=np.array([bank.insert(v)for v in q['vertices']],int)
 nativeids={r['original_id']:r for r in g['classifiers']if r['role'].endswith('outer')};keep=[];lineage=[];normals=[];uvs=[]
 for i,(oid,fid)in enumerate(zip(q['original_id'],q['source_face_id'])):
  row=nativeids[int(oid)]
  if int(fid)>=row['source_faces']:continue
  role=row['role'].split('-')[0];source_id=int(z[role+'_local_patch_face_ids'][int(fid)]);source=sources[role];tri=source.triangles[source_id];points=bank.positions[mapped[q['faces'][i]]];weights=barycentric(points,tri)
  assert np.max(np.linalg.norm(weights@tri-points,axis=1))<2e-12 and weights.min()>-1e-8 and weights.max()<1+1e-8
  keep.append(i);lineage.append([0 if role=='parent'else 1,source_id]);normals.append(weights@source.attributes['normals'][source.faces[source_id]])
  if 'uvs'in source.attributes:uvs.append(weights@source.attributes['uvs'][source.faces[source_id]])
 localfaces=mapped[q['faces'][keep]];parts=[];fixed=[];offset=0;fixed_attributes={}
 for role,source in sources.items():
  ids=np.setdiff1d(np.arange(len(source.faces)),z[role+'_local_patch_face_ids']);fs=source.faces[ids]+offset;parts.append(fs);fixed.append(dict(role=role,source_face_ids=ids.tolist(),literal_original_face_indices_retained=True,literal_original_vertices_retained=True,literal_original_packed_attributes_retained=True))
  for name,value in source.attributes.items():fixed_attributes[role+'_source_'+name]=value
  offset+=len(source.vertices)
 faces=np.concatenate(parts+[localfaces]);v=bank.positions;mesh=Mesh(v,faces,FRAME);mesh.save(OUT/'source-outer-assembly-with-repaired-root.npz');local=Mesh(v,localfaces,FRAME);local.save(OUT/'new-local-outer.npz');old=Mesh(np.r_[z['parent_vertices'],z['branch_vertices']],np.r_[z['parent_faces'][z['parent_local_patch_face_ids']],z['branch_faces'][z['branch_local_patch_face_ids']]+len(z['parent_vertices'])],FRAME);old.save(OUT/'old-local-outer.npz')
 atomic_npz(OUT/'literal-splice-ancestry.npz',vertices=v,faces=faces,local_faces=localfaces,local_source_part_and_face=np.array(lineage,int),local_corner_interpolated_source_normal_words=np.array(normals),**fixed_attributes,frame=np.array(FRAME))
 assert np.array_equal(v[:len(prefix)],prefix);topo=topology(mesh);assert len(topo['components'])==2 and topo['open_edges']==60 and not topo['nonmanifold_edges'] and not topo['inconsistent_winding_edges'] and not topo['duplicate_faces'] and not topo['zero_area_faces']
 # Compare only new local outer versus fixed continuations/secondary course;
 # declared reciprocal port edges may touch but may not cross.
 fixedmesh=Mesh(v,np.concatenate(parts),FRAME);contact,summary=pair_contacts(local,fixedmesh,0.);assert summary['complete_within_pair_budget'] and summary['crossing_face_pairs']==0
 # A literal source interpolation transports normals without changing the
 # old normal words. Bound its length and retain the true outer-only scope.
 normal_lengths=np.linalg.norm(normals,axis=2);assert np.isfinite(normal_lengths).all()and normal_lengths.min()>0
 proof=dict(schema='literal-native-M3-postcentral-local-outer-splice/v1',installed=False,passed=True,source_face_removal=dict(parent=192,branch=96),new_local_outer_faces=len(localfaces),fixed_source_faces=len(np.concatenate(parts)),literal_source_prefix_vertices=len(prefix),literal_source_prefix_preserved=True,fixed_source_rows=fixed,all_new_outer_triangles_are_subsets_of_literal_source_triangles=True,maximum_source_roundoff_snap_m=bank.maximum_snap,topology=topo,new_local_vs_fixed_contacts=summary,interpolated_corner_normal_word_length_range=[float(normal_lengths.min()),float(normal_lengths.max())],scope='Only a source-provenance outer-surface splice plus separately derived local finite wall. Complete source continuations remain outer sheets. Secondary branch/main-branch 58 pre-existing crossing pairs are not repaired.',remaining=['Joined packed-normal and geometric tangent continuity review across all three ports.','Actual rendering/sections and local recipient context.','Source secondary course self-crossing and whole-vessel finite walls.'],elapsed_seconds=time.monotonic()-start)
 atomic_json(OUT/'splice-proof.json',proof);print(json.dumps({k:proof[k]for k in ['passed','new_local_outer_faces','fixed_source_faces','maximum_source_roundoff_snap_m','elapsed_seconds']}))

if __name__=='__main__':main()
