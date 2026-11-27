"""Apply only point identifications already forced by literal Float32 storage."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--retire-zero-area-cut-charts',action='store_true');args=ap.parse_args()
h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.input;out=h/args.stage;out.mkdir(exist_ok=False);start=time.monotonic()
m=json.loads((src/'ownership.json').read_text());z=np.load(src/'arrangement.npz');arrays={k:z[k] for k in z.files};p=z['positions'];s={r['id']:z[r['faces']] for r in m['surfaces']}
retired=[]
if args.retire_zero_area_cut_charts:
 for sid,fs in list(s.items()):
  q=p[fs].astype('<f4').astype(float)
  if np.any(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])!=0):continue
  assert '/mask-' in sid or '/mask-interface' in sid, 'Only newly partitioned cut charts may be retired'
  tri=p[fs];retired.append(dict(surface=sid,source_area_m2=float(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()/2),stored_area_m2=0.,faces=len(fs)))
  del s[sid]
 gone={r['surface']for r in retired};m['surfaces']=[r for r in m['surfaces']if r['id']not in gone]
 m['cells']={pid:[u for u in uses if u['surface']not in gone]for pid,uses in m['cells'].items()}
 assert all(m['cells'].values()), 'A complete tissue cell vanished'
surface,parents,proof=conform_quantized_surfaces(p,s,np.zeros(3))
zero=0
for r in m['surfaces']:
 sid=r['id'];f=surface[sid];q=p[f].astype('<f4').astype(float);keep=np.any(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])!=0,axis=1);zero+=int((~keep).sum())
 assert keep.any(), 'A complete chart vanished; needs explicit anatomical review'
 arrays[r['faces']]=f[keep]
 for k in ['normals','uvs']:arrays[r[k]]=z[r[k]][parents[sid][keep]]
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'ownership.json',m)
proof.update(positions_literal=True,source_stage=args.input,source_metadata_sha256=digest_file(src/'ownership.json'),
 collinear_stored_faces_removed=zero,elapsed_seconds=time.monotonic()-start,
 numerical_and_anatomical_acceptance_pending=True,installed=False)
proof['explicit_zero_area_cut_chart_retirements']=retired
atomic_json(out/'storage-conformity.json',proof);print(json.dumps(proof),flush=True)
