"""Only coplanar interior diagonal changes; retain every boundary/control word."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.float32_diagonal_repair import repair_diagonals
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--fresh-owned-interface-normals',action='store_true');args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.input;out=h/args.stage;out.mkdir(exist_ok=False);meta=json.loads((src/'ownership.json').read_text());z=np.load(src/'arrangement.npz');arrays={k:z[k].copy()for k in z.files};p=arrays['positions'];q=p.astype('<f4').astype(float);rows=[];start=time.monotonic()
for r in meta['surfaces']:
 f=arrays[r['faces']];n=arrays[r['normals']];a=p[f];b=q[f];nn=np.cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);qn=np.cross(b[:,1]-b[:,0],b[:,2]-b[:,0]);bad=(np.einsum('fi,fi->f',nn,qn)<=0)|np.any(np.einsum('fki,fi->fk',n.astype(float),qn)<=0,axis=1)
 if not np.any(bad):continue
 if time.monotonic()-start>25:raise TimeoutError('Diagonal stage deadline')
 try:
  fresh=args.fresh_owned_interface_normals and r['id'].startswith('owned-interface:')
  result,proof=repair_diagonals(p,f,np.zeros(3),corner_normals=None if fresh else n)
  for change in proof['changes']:
   ids=change['faceIds'];sources=change['newCornerSources']
   for key in [r['normals'],r['uvs']]:
    before=arrays[key][ids].copy();arrays[key][ids]=np.array([[before[s[0],s[1]]for s in row]for row in sources])
  arrays[r['faces']]=result;row=dict(surface=r['id'],passed=True,source_bad_faces=int(bad.sum()),**proof)
  if fresh:
   from geometry_tools.tissue import normal_words
   arrays[r['normals']]=normal_words(q,{r['id']:result})[r['id']]
   row['new_internal_chart_normals_recomputed_on_stored_geometry']=True
 except ValueError as exc:row=dict(surface=r['id'],passed=False,source_bad_faces=int(bad.sum()),error=str(exc))
 rows.append(row);print(r['id'],row['passed'],row.get('diagonalsChanged',row.get('error')),flush=True)
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'ownership.json',meta);atomic_json(out/'diagonals.json',dict(passed=all(r['passed']for r in rows),source_metadata_sha256=digest_file(src/'ownership.json'),source_stage=args.input,positions_literal=True,source_boundary_edges_literal=True,rows=rows,elapsed_seconds=time.monotonic()-start,installed=False))
