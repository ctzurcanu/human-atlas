"""Separate bounded emitted-cell or exact shared-boundary/union checks."""
from pathlib import Path
import json,sys,time
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,topology
from geometry_tools.algorithms.finite_tissue_boolean import same_oriented_triangles
from geometry_tools.algorithms.shared_cell_disjoint import shared_cell_disjoint
C=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3/m3-postcentral-port-feasibility-v1/three-course-finite-educational-v1';P=C/'packed-owned-cells-v1';FRAME='model:local-ta98:assembled-metres'

def cell(name):
 start=time.monotonic();out=P/(name+'-intrinsic-v1');out.mkdir(exist_ok=False);mesh=Mesh.load(P/(name+'-decoded.npz'));result=audit(mesh,require_closed=True,expected_components=1,expected_genus=[2]);result['elapsed_seconds']=time.monotonic()-start;atomic_json(out/'audit.json',result);print(json.dumps(dict(cell=name,passed=result['passed'],failures=result['failures'],defect_pairs=result['self_intersections'].get('defect_pairs'),elapsed_seconds=result['elapsed_seconds'])));assert result['passed']

def union():
 start=time.monotonic();out=P/'shared-union-v1';out.mkdir(exist_ok=False);a=Mesh.load(P/'M3-decoded.npz');b=Mesh.load(P/'postcentral-decoded.npz');layout=json.loads((P/'packed-layout.json').read_text());origin=np.array(layout['geometry_origin']);bank=(P/'cells.bin').read_bytes();record=next(iter(layout['cells'].values()));positions=np.ndarray((record['vertexCount'],3),'<f4',buffer=bank,offset=record['positions']).astype(float)+origin;external=[];interfaces=[]
 for name,row in layout['cells'].items():
  indices=np.ndarray((row['indexCount'],),'<u4',buffer=bank,offset=row['indices'])
  for group in row['tissueDrawGroups']:
   faces=indices[group['start']:group['start']+group['count']].reshape(-1,3)
   (interfaces if group.get('neighbor')else external).append(faces)
 ext=np.concatenate(external);inter=np.concatenate(interfaces);whole=Mesh.load(C/'assembly-v1/wall.npz');encoded=(whole.vertices-origin).astype('<f4').astype(float)+origin;_,inverse=np.unique(np.r_[positions,encoded],axis=0,return_inverse=True);left=inverse[:len(positions)];right=inverse[len(positions):];assert same_oriented_triangles(left[ext],right[whole.faces]);interface=left[inter]
 def oriented(t):
  r=tuple(map(int,t));return min(r,r[1:]+r[:1],r[2:]+r[:2])
 forward=Counter(oriented(t)for t in interface[:len(interface)//2]);reverse=Counter(oriented(t[[0,2,1]])for t in interface[len(interface)//2:]);assert forward==reverse
 emittedwhole=Mesh(positions,ext,FRAME);emittedwhole.save(out/'emitted-assembled-wall-no-internal-caps.npz');t=topology(emittedwhole);assert t['closed_oriented_manifold']and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==3
 # The independent disjointness helper includes native libigl winding. Save
 # all exact input meshes and pending state before that C++ operation.
 atomic_json(out/'disjoint-check-pending.json',dict(inputs=[str(P/'M3-decoded.npz'),str(P/'postcentral-decoded.npz')],method='Complete triangle contacts and bidirectional nonshared-face generalized winding'))
 disjoint=shared_cell_disjoint(a,b);atomic_json(out/'disjoint-proof.json',disjoint);assert disjoint['passed'],disjoint.get('reason')
 volumes=[topology(m)['oriented_triangle_integral_m3']for m in [a,b]];wholevolume=t['oriented_triangle_integral_m3'];residual=abs(sum(volumes)-wholevolume);assert residual<=max(abs(wholevolume)*1e-10,1e-20)
 result=dict(schema='emitted-owned-cells-exact-union/v1',passed=True,installed=False,cell_volumes_m3=volumes,assembled_whole_volume_m3=wholevolume,volume_sum_residual_m3=residual,external_boundary_equals_encoded_validated_whole_triangle_for_triangle=True,assembled_wall_faces=len(ext),single_canonical_interface_faces=len(interface)//2,interface_uses_exactly_opposed=True,disjointness=disjoint,assembled_topology=t,global_native_terminal_annuli=4,no_internal_interface_in_assembled_display=True,scope='Actual emitted Float32 cells, shared interface and external union. Educational inner thickness and provisional station-176 identity boundary remain qualified.',elapsed_seconds=time.monotonic()-start);atomic_json(out/'union-proof.json',result);print(json.dumps({k:result[k]for k in ['passed','assembled_wall_faces','single_canonical_interface_faces','volume_sum_residual_m3','elapsed_seconds']}))

if __name__=='__main__':
 if len(sys.argv)!=2 or sys.argv[1]not in ['M3','postcentral','union']:raise SystemExit('Choose one separately capped stage: M3, postcentral or union')
 union()if sys.argv[1]=='union'else cell(sys.argv[1])
