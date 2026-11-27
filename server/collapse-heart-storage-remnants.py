"""Remove bounded cut slivers with every affected tissue link checked."""
from pathlib import Path
import argparse,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.shared_edge_collapse import collapse_short_edges
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--reference',required=True);ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args()
h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.input;out=h/args.stage;out.mkdir(exist_ok=False)
m=json.loads((src/'ownership.json').read_text());z=np.load(src/'arrangement.npz');p=z['positions'];desired=np.load(h/args.reference/'arrangement.npz')['positions'];prefix=len(np.load(h/args.source/'source-arrangement.npz')['positions'])
surfaces={r['id']:z[r['faces']] for r in m['surfaces']};normals={r['id']:z[r['normals']] for r in m['surfaces']}
r=collapse_short_edges(p,surfaces,m['cells'],normals,desired_positions=desired,maximum_displacement_m=2e-8,
 maximum_steps=12,wall_seconds=20.,candidate_face_rings=0,preserve_material_membership=False,locked=np.arange(prefix))
arrays={k:z[k] for k in z.files}
for row in m['surfaces']:
 sid=row['id'];assert len(r['surfaces'][sid]), 'A complete chart vanished'
 arrays[row['faces']]=r['surfaces'][sid]
 for key in ['normals','uvs']:arrays[row[key]]=z[row[key]][r['source_faces'][sid]]
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'ownership.json',m)
proof=r['evidence'];proof.update(source_metadata_sha256=digest_file(src/'ownership.json'),immutable_reference=args.reference,
 whole_geometry_audits_required=True,installed=False)
atomic_json(out/'contractions.json',proof);print(json.dumps(proof),flush=True)
