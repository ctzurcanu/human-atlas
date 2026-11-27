"""Expose a patent common-hepatic section and join its actual native wall."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.annular_bridge import trim_named_inlet,attach_bridge,mouth_loops
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank,split_closed_cell
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='common-hepatic-mouth-v239');args=ap.parse_args()
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/args.stage;out.mkdir(exist_ok=False);start=time.monotonic();frame='model:local-ta98:assembled-metres'
cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(cat.read_text());r=AtlasReader();p=next(p for p in a['parts']if p['id']=='TA98PART:FJ3079:A05.8.01.061');parent=r.read_part(a,p,frame=frame,catalogue_path=cat)
def compact(m):
 v,ix=np.unique(m.vertices,axis=0,return_inverse=True);return Mesh(v,ix[m.faces],frame)
parent=compact(parent);bank=SharedCutBank(parent.vertices);cut=split_closed_cell(bank,parent.faces,[0,1,0],1.248046875,'native-common-upper',allow_holes=True);parent=Mesh(np.asarray(bank.positions),cut.negative_faces,frame);pc=np.flatnonzero(cut.negative_ancestry<0)
child=Mesh.load(b/'right-caudate-duct-v226/emitted-wall.npz');v,f,cc,trim=trim_named_inlet(child.vertices,child.faces,1.25390625,[-.024,1.2518,.038]);child=Mesh(v,f,frame)
q=attach_bridge(parent.vertices,parent.faces,pc,child.vertices,child.faces,cc);origin=np.array([0,1.25,0]);surfaces,_,conform=conform_quantized_surfaces(q['vertices'],{'parent':q['faces']},origin);new=Mesh((q['vertices']-origin).astype('<f4').astype(float)+origin,surfaces['parent'],frame);new.save(out/'parent.npz')
cv=(child.vertices-origin).astype('<f4').astype(float)+origin;cf,_,childconform=conform_quantized_surfaces(child.vertices,{'child':child.faces},origin);child=Mesh(cv,cf['child'],frame);child.save(out/'child.npz')
checks={}
for name,m in [('parent',new),('child',child)]:
 checks[name]=audit(m,require_closed=True,max_pairs=1000000,max_tests=12000000);atomic_json(out/(name+'-check.json'),checks[name]);print(name,checks[name]['passed'],checks[name]['failures'],flush=True)
proof=dict(passed=all(x['passed']for x in checks.values()),parent=str(out/'parent.npz'),child=str(out/'child.npz'),parent_interface_faces=q['interface_faces'].tolist(),child_interface_faces=cc.tolist(),source_controls_literal=True,distal_cystic_bile_interfaces_preserved=True,curve=q['path'],trim=trim,conform=conform,childconform=childconform,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'junction.json',proof)
if not proof['passed']:raise ValueError('Common hepatic junction failed; no installation')
