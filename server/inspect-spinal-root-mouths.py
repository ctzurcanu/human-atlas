"""Record actual native nerve mouths against the fitted pial surface.

Boundary loops are source geometry, not inferred dorsal/ventral identities or
certified joins. Peripheral nerves and the native skeleton remain unchanged.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';bank=base/'female-31-segment-bank-v3';d=json.loads((bank/'tissue.json').read_text());assert d['passed'];raw=(bank/'cells.bin').read_bytes()
faces=[]
for name,p in d['parts'].items():
    if not name.startswith('pia@'):continue
    positions=np.ndarray((p['vertexCount'],3),'<f4',buffer=raw,offset=p['positions']).astype(float)+[0.,1.35,-.05]
    ff=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices'])
    for g in p['tissueDrawGroups']:
        if not g.get('neighbor'):faces.append(ff[g['start']//3:(g['start']+g['count'])//3])
frame=d['frame'];pia=Mesh(positions,np.concatenate(faces),frame);index=BVH(pia)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';atlas=json.loads(path.read_text());by_id={p['id']:p for p in atlas['parts']};reader=AtlasReader()
source=json.loads((base/'female-root-regions-v2/root-regions.json').read_text());rows=[]
for record in source['root_patches']:
    p=by_id[record['part_id']];mesh=reader.read_part(atlas,p,frame=frame,catalogue_path=path)
    vertices,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces];ee=np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]
    keys,which,counts=np.unique(np.sort(ee,axis=1),axis=0,return_inverse=True,return_counts=True);boundary=ee[counts[which]==1]
    nxt=defaultdict(list)
    for a,b in boundary:nxt[int(a)].append(int(b))
    branched=any(len(x)!=1 for x in nxt.values())or set(nxt)!={int(x)for x in boundary[:,1]};loops=[]
    if not branched:
        todo=set(nxt)
        while todo:
            first=min(todo);ids=[];i=first
            while i in todo:todo.remove(i);ids.append(i);i=nxt[i][0]
            assert i==first
            points=vertices[ids];center=points.mean(0);q,_=index.closest_points(center[None],backend='bvh')
            loops.append(dict(vertex_ids=ids,center_m=center.tolist(),radius_m=float(np.linalg.norm(points-center,axis=1).mean()),
                axial_span_m=float(np.ptp(points[:,1])),candidate_proximal_mouth=bool(abs(center[0])<.02),
                nearest_pial_point_m=q['points'][0].tolist(),pial_distance_m=float(q['distance_m'][0])))
    filename=out/(p['id'].replace(':','_')+'.npz');atomic_npz(filename,vertices=vertices,faces=f,frame=np.array(frame))
    rows.append(dict(part_id=p['id'],level=record['level'],source_geometry_sha256=mesh.geometry_sha,source_mesh=str(filename),
        triangles=len(f),boundary_edges=len(boundary),branched_boundary=branched,nonmanifold_edges=int((counts>2).sum()),loops=loops,certified_interface=False))
atomic_json(out/'root-mouths.json',dict(source_catalogue_sha256=digest_file(path),pial_bank_sha256=d['bank_sha256'],roots=rows,
    geometry_changed=False,skeleton_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(root_sources=len(rows),proximal_mouths=sum(l['candidate_proximal_mouth']for r in rows for l in r['loops']),
    branched_sources=[r['part_id']for r in rows if r['branched_boundary']],elapsed_seconds=time.monotonic()-start)),flush=True)
