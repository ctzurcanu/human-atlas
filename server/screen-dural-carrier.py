"""Actual fixed-bone crossings and central-root witnesses for a dural candidate."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import SurfaceIndex

ap=argparse.ArgumentParser();ap.add_argument('stage');args=ap.parse_args();stage=Path(args.stage).resolve();start=time.monotonic()
z=np.load(stage/'surface.npz');carrier=Mesh(z['vertices'],z['faces'],'model:local-female-ta98:assembled-metres');surface=SurfaceIndex(carrier)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(path.read_text());reader=AtlasReader();rows=[]
for p in a['parts']:
    if p.get('suppressed')or not p['id'].startswith('LOCAL:female:'):continue
    if 'lumbar_vertebra_'not in p['id']and not p['id'].startswith('LOCAL:female:sacrum:'):continue
    bone=reader.read_part(a,p,frame=carrier.frame,catalogue_path=path)
    result,summary=pair_contacts(carrier,bone,0.,max_pairs=150000,max_tests=1000000)
    assert not summary['truncated'],'Bone crossing screen exceeded its pair budget'
    bad=np.isin(result['relation'],['crossing','coplanar-overlap']);pairs=result['pairs'][bad]
    rows.append(dict(part_id=p['id'],geometry_sha256=bone.geometry_sha,crossings=int(bad.sum()),
        examples=pairs[:12].tolist(),candidate_face_y_range_m=None if not len(pairs)else[float(x)for x in np.quantile(carrier.triangles[pairs[:,0],:,1],[0,1])],summary=summary))
    print(json.dumps(dict(bone=p['name'],crossings=int(bad.sum()),elapsed_seconds=time.monotonic()-start)),flush=True)
root_path=stage.parent/'female-covering-source-v1/cauda_equina_base-845.npz';root=np.load(root_path);points=np.unique(root['vertices'],axis=0);points=points[(np.abs(points[:,0])<.003)&(points[:,1]<1.178)]
query,evidence=surface.closest_points(points);tri=carrier.triangles[query['face_ids']];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
signed=np.einsum('ni,ni->n',points-query['points'],normal);outside=signed>.00005
atomic_json(stage/'fixed-context-screen.json',dict(candidate_sha256=digest_file(stage/'surface.npz'),bone_screens=rows,
    central_root_clues=dict(points=len(points),outside_oriented_surface=int(outside.sum()),
        maximum_oriented_distance_m=float(signed.max()),outside_y_m=points[outside,1].tolist(),
        scope='Closest oriented surface witness, not a closed-volume containment certificate'),
    passed=not any(r['crossings']for r in rows)and not outside.any(),skeleton_changed=False,
    anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
