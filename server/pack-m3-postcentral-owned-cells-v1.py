"""Encode the explicit two-cell arrangement once; reject bad Float32 charts."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.tissue import decode_packed
C=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3/m3-postcentral-port-feasibility-v1/three-course-finite-educational-v1';A=C/'owned-cells-v1';OUT=C/'packed-owned-cells-v1'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);m=json.loads((A/'owned-arrangement.json').read_text());z=np.load(A/'owned-arrangement.npz');positions=z['positions'];origin=np.array(m['geometry_origin']);surfaces={r['id']:z[r['array']]for r in m['surfaces']};words={r['id']:z[r['corner_words']]for r in m['surfaces']};uvs={r['id']:z[r['corner_uvs']]for r in m['surfaces']};encoded=(positions-origin).astype('<f4').astype(float)+origin;preflight=[]
 for key,f in surfaces.items():
  t=encoded[f];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);dot=np.einsum('tki,ti->tk',words[key].astype(float),n);bad=np.argwhere(dot<=0);preflight.append(dict(surface=key,faces=len(f),nonpositive_corners=len(bad),bad_face_corner_ids=bad.tolist(),minimum_dot=float(dot.min())))
 atomic_json(OUT/'Float32-outward-preflight.json',dict(origin=origin.tolist(),surfaces=preflight));assert not any(r['nonpositive_corners']for r in preflight),'An emitted Float32 chart has a collapsed or non-outward corner; stop without packing'
 packed=pack_finite_cells(positions,surfaces,words,m['cells'],m['part_ids'],origin,corner_uvs=uvs);atomic_bytes(OUT/'cells.bin',packed.data);atomic_json(OUT/'packed-layout.json',dict(frame=m['frame'],geometry_origin=origin.tolist(),cells=packed.cells,surfaces={k:v.tolist()for k,v in packed.surfaces.items()},source_arrangement_sha256=digest_file(A/'owned-arrangement.npz'),source_manifest_sha256=digest_file(A/'owned-arrangement.json'),identity_boundary=m['identity_boundary']))
 decoded=decode_packed(packed,m['frame'],origin)
 for key,mesh in decoded.items():mesh.save(OUT/(key+'-decoded.npz'))
 # The stored surface bank has literal words/UVs, not regenerated shading.
 for key,f in surfaces.items():
  indices=packed.surfaces[key];assert np.array_equal(packed.normal_words[indices],words[key]);uv_offset=next(iter(packed.cells.values()))['uvs'];stored_uv=np.ndarray((len(packed.positions),2),'<f4',buffer=packed.data,offset=uv_offset);assert np.array_equal(stored_uv[indices],uvs[key]);assert np.array_equal(packed.positions[indices].astype(float)+origin,encoded[f])
 shared=np.unique(packed.surfaces['shared:station176']);signs=[]
 for key in ['M3','postcentral']:
  row=packed.cells[key];signs.append(np.ndarray((row['vertexCount'],),'i1',buffer=packed.data,offset=row['normalSigns'])[shared]);assert len([g for g in row['tissueDrawGroups']if g.get('neighbor')])==1
 assert np.array_equal(signs[0],-signs[1]);report=dict(schema='packed-two-owned-M3-postcentral-cells/v1',installed=False,bytes=len(packed.data),stored_attribute_controls=len(packed.positions),cell_triangles={k:len(v.faces)for k,v in decoded.items()},single_shared_interface_faces=len(packed.surfaces['shared:station176']),shared_coordinate_normal_words_identical=True,reciprocal_opposed_normal_signs=True,native_outer_normal_words_and_UVs_preserved=True,no_visible_assembled_internal_interface='Existing tissueDrawGroups neighbor rule; publication adapter must enable paired normal signs on arterial material.',identity_boundary=m['identity_boundary'],remaining=['Actual emitted cell topology/self-crossing, disjointness, aggregate boundary/volume and blood path proof.','Actual packed-normal renders/sections and fixed current context before installation.'],elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'packing-proof.json',report);print(json.dumps({k:report[k]for k in ['bytes','stored_attribute_controls','cell_triangles','single_shared_interface_faces','elapsed_seconds']}))

if __name__=='__main__':main()
