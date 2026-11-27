"""Reproduce every originally audited buffer before orientation compaction.

No mesh is newly accepted here. All 97 full geometry and attribute hashes
must match their saved prechecks. The installer can then prove literal used
corner equality between those buffers and the corrected compact output.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.tissue import normal_words,decode_packed
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--prechecks',required=True);ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');origin=np.array([0.,1.35,-.05]);local=(z['positions']-origin).astype('<f4');encoded,mapping=np.unique(local.astype(float),axis=0,return_inverse=True);surfaces={}
for sid in d['surfaces']:
    f=mapping[z[sid]];q=encoded[f];keep=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0;surfaces[sid]=f[keep]
words=normal_words(encoded,surfaces)
for sid,f in surfaces.items():
    q=encoded[f];fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);bad=np.einsum('tki,ti->tk',words[sid].astype(float),fn)<=0
    if bad.any():
        facet=np.rint(fn/np.linalg.norm(fn,axis=1)[:,None]*32767).astype('<i2');words[sid][bad]=np.broadcast_to(facet[:,None],words[sid].shape)[bad]
ids={n:'TA98CORD:female:'+n for n in d['cells']};frame='model:local-female-ta98:assembled-metres'
packed=pack_finite_cells(encoded+origin,surfaces,words,d['cells'],ids,origin,share_orientation_buffers=False)
checks={}
prefix=Path(args.prechecks)
for path in prefix.parent.glob(prefix.name+'-precheck-*.json'):
    record=json.loads(path.read_text());assert record['source_sha256']==digest_file(source/'cells.npz');checks.update(record['checks']['cells'])
meshes=decode_packed(packed,frame,origin);assert set(checks)==set(meshes)
for name,mesh in meshes.items():
    assert mesh.geometry_sha==checks[name]['geometry_sha256'],name
    assert mesh.attributes_sha==checks[name]['attributes_sha256'],name
bank=out/'original-audited-buffers.bin';atomic_bytes(bank,packed.data)
atomic_json(out/'audit-layout.json',dict(frame=frame,bank=str(bank),bank_sha256=digest_file(bank),parts={p['id']:p for p in packed.cells.values()},origin=origin.tolist(),
    reproduced_precheck_attributes_and_geometry=True,source_cell_bank_sha256=digest_file(source/'cells.npz'),
    reference_only=True,not_an_installable_candidate=True,original_failed_cells_retained_as_unaccepted=[n for n,c in checks.items()if not c['passed']],elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(reproduced_cells=len(meshes),bytes=len(packed.data),reference_only=True,elapsed_seconds=time.monotonic()-start)),flush=True)
