"""One measured similarity for the entire coordinated female wall template."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
ap=argparse.ArgumentParser();ap.add_argument('--source',default='regional-clean-v101');ap.add_argument('--stage',default='female-wall-v105');ap.add_argument('--condition-split-stars',action='store_true');args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/args.source;out=b/args.stage;out.mkdir(exist_ok=False);z=np.load(source/'cells.npz');metadata=json.loads((source/'construction.json').read_text());fit=json.loads((b/'female-fit-v68/fit.json').read_text());r=np.asarray(fit['rotation']);s=fit['scale'];t=np.asarray(fit['translation_m']);oldorigin=np.array([0,1.375,.03125]);origin=np.array([0.,1.3125,.03125])if args.condition_split_stars else s*(r@oldorigin)+t;ideal=s*(z['positions']@r.T)+t-origin;p=ideal.astype('<f4').astype(float)+origin;names=metadata['names'];remap=np.arange(len(p));collapses=[]
if args.condition_split_stars:
 from heart_lattice_cleanup import collapse_shared_stars
 base=np.load(b/'chart-coats-v82/cells.npz')['positions'];dtype=np.dtype((np.void,12));fixed=np.isin((z['positions']-oldorigin).astype('<f4').view(dtype).ravel(),(base-oldorigin).astype('<f4').view(dtype).ravel());original_cells=[z['cell_'+str(i)]for i in range(len(names))];refs=[]
 for faces in original_cells:
  tri=z['positions'][faces]-oldorigin;refs.append(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])@r.T*s*s)
 collapsed,cleanup=collapse_shared_stars(p,original_cells,refs,fixed,seconds=20.);collapses=cleanup['collapses'];print('Link-guarded calculation-star collapses',len(collapses),flush=True)
else:collapsed=[z['cell_'+str(i)]for i in range(len(names))]
q=normalize_shared_boundary(p,{str(i):collapsed[i]for i in range(len(names))},origin);arrays={'positions':q['positions']};rows=[]
for i,name in enumerate(names):
 faces=q['cells'][str(i)];arrays['cell_'+str(i)]=faces;used,inv=np.unique(faces,return_inverse=True);m=Mesh(arrays['positions'][used],inv.reshape(-1,3),'model:local-female-ta98:assembled-metres');m.save(out/('cell-'+str(i)+'.npz'));proof=audit(m,self_intersections=False,require_closed=True);atomic_json(out/('topology-'+str(i)+'.json'),proof);rows.append(dict(name=name,passed=proof['passed'],failures=proof['failures']));print(name,proof['passed'],proof['failures'],flush=True)
atomic_npz(out/'cells.npz',**arrays);atomic_json(out/'encoding.json',q['evidence']);atomic_json(out/'construction.json',dict(passed=all(v['passed']for v in rows),source_sha256=digest_file(source/'cells.npz'),names=names,origin=origin.tolist(),fit=fit,maximum_position_encoding_error_m=float(np.linalg.norm(p-origin-ideal,axis=1).max()),calculation_star_collapses=collapses,all_layers_and_septa_use_one_transform=True,skeleton_moved=False,full_emitted_audits_pending=True,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False))
