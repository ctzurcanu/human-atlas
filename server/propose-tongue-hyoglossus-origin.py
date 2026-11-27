"""Save literal native body/greater-horn surface guides for HG reconstruction.

This deliberately exports an origin SURFACE, never labels an open bone patch
as a finished muscle. No fixed skeleton controls, normals or faces are moved.
"""
from pathlib import Path
from collections import defaultdict
import heapq
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

SOURCE = ROOT/'.local-models/geometry-workbench/tongue-completion/source-v1/meshes'
OUT = ROOT/'.local-models/geometry-workbench/tongue-completion/native-attachment-candidate-v1'
FRAME = 'model:local-ta98:assembled-metres'
START = time.monotonic()


def load(i):
    return Mesh.load(SOURCE/f'native-{i:02d}.npz')


def save_patch(v, f, indices, owners, side):
    used, inv = np.unique(f[indices], return_inverse=True)
    patch = Mesh(v[used], inv.reshape(-1, 3), FRAME, source={
        'classification': 'Anatomical origin surface proposal on unchanged native bone',
        'bone_modified': False, 'anatomical_acceptance': False,
    })
    patch.save(OUT/f'hyoglossus-{side}-bone-origin-surface.npz')
    return dict(mesh=str(OUT/f'hyoglossus-{side}-bone-origin-surface.npz'),
        actual_source_face_owners=[owners[i] for i in indices], topology=topology(patch))


def main():
    rows=[]
    for side, sign, hi, mi in [('left', 1, 16, 3), ('right', -1, 17, 4)]:
        body, horn, muscle = load(13), load(hi), load(mi)
        v, inv = np.unique(np.vstack([body.vertices, horn.vertices]), axis=0, return_inverse=True)
        f = np.vstack([inv[body.faces], inv[horn.faces+len(body.vertices)]])
        owners = [dict(part=body.source['record']['id'], face=int(i)) for i in range(len(body.faces))]
        owners += [dict(part=horn.source['record']['id'], face=int(i)) for i in range(len(horn.faces))]
        tri=v[f]; c=tri.mean(1); n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        # Source facet witnesses, mirrored by actual native surface proximity.
        target0=np.array([sign*.00866382,1.6072648,.03739575])
        target1=np.array([sign*.02124292,1.6125284,.01410411])
        allowed=(sign*c[:,0]>.004)&(n[:,1]>.05)
        seed=int(np.argmin(np.where(allowed,np.linalg.norm(c-target0,axis=1),np.inf)))
        terminal=int(np.argmin(np.where(allowed,np.linalg.norm(c-target1,axis=1),np.inf)))
        edges=defaultdict(list)
        for i,t in enumerate(f):
            for a,b in zip(t,np.roll(t,-1)):edges[tuple(sorted([int(a),int(b)]))].append(i)
        adj=defaultdict(set)
        for values in edges.values():
            for i in values:adj[i].update(j for j in values if j!=i)
        dist={seed:0.}; parent={}; queue=[(0.,seed)]
        while queue:
            if time.monotonic()-START>15:raise TimeoutError('HG origin proposal exceeded 15 second native cap')
            cost,i=heapq.heappop(queue)
            if cost!=dist[i]:continue
            if i==terminal:break
            for j in adj[i]:
                if not allowed[j]:continue
                value=cost+np.linalg.norm(c[j]-c[i])*(1+2*(1-n[j,1]))
                if value<dist.get(j,np.inf):dist[j]=value;parent[j]=i;heapq.heappush(queue,(value,j))
        if terminal not in dist:raise ValueError('No connected upper body-to-greater-horn surface path in literal source')
        path=[terminal]
        while path[-1]!=seed:path.append(parent[path[-1]])
        path=path[::-1]
        # One-ring around the actual connected upper-surface path forms an
        # explicit footprint proposal; patch topology is reported, not assumed.
        selected=set(path)
        for i in path:selected.update(j for j in adj[i] if allowed[j] and n[j]@n[i]>.25)
        indices=sorted(selected)
        patch=save_patch(v,f,indices,owners,side)
        preserve=np.flatnonzero(np.all(muscle.triangles[:,:,1]>=1.621,axis=1))
        Mesh(muscle.vertices,muscle.faces[preserve],FRAME,muscle.attributes,
             {'classification':'Literal native superior HG surface preservation constraint, open region',
              'parent_geometry_sha256':muscle.geometry_sha}).save(OUT/f'hyoglossus-{side}-preserved-upper.npz')
        rows.append(dict(side=side, origin_surface=patch,
            actual_bone_path=[dict(**owners[i], centroid=c[i].tolist(),outward_normal=n[i].tolist())for i in path],
            path_length_m=float(sum(np.linalg.norm(c[b]-c[a])for a,b in zip(path,path[1:]))),
            origin_surface_bounds=np.array([tri[indices].min((0,1)),tri[indices].max((0,1))]).tolist(),
            fixed_native_superior_constraint_faces=preserve.tolist(),
            reconstruction_required='Replace narrow medial lower origin by a broad curved sheet rooted along actual upper hyoid body and greater horn surface; merge into literal superior native lateral-tongue insertion.',
            reasons_native_local_snap_is_inadequate=['86/87 body crossing pairs and 1.100 mm sampled bone penetration.',
                'Actual closest greater-horn gap is 4.4185 mm; no greater-horn contact exists.',
                'Bone path spans approximately 30 mm in a different plane from the current short medial tip.'],
            full_muscle_candidate_constructed=False, anatomical_acceptance=False))
    report=dict(schema='hyoglossus-native-origin-proposal/v1', elapsed_seconds=time.monotonic()-START,
        frame=FRAME,skeleton_modified=False,installed=False,proposals=rows)
    (OUT/'hyoglossus-origin-proposal.json').write_text(json.dumps(report,indent=2)+'\n')
    spec=json.loads((OUT/'review-spec.json').read_text())
    for row in rows:
        spec['meshes'].append(dict(mesh=row['origin_surface']['mesh'],name=f"HG proposed origin surface ({row['side']})",color=[.95,.72,.18,1]))
    (OUT/'review-spec.json').write_text(json.dumps(spec,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items()if k!='proposals'}))
    print(json.dumps([dict(side=r['side'],surface=r['origin_surface']['mesh'],path_length_m=r['path_length_m'],bounds=r['origin_surface_bounds'])for r in rows],indent=2))


if __name__=='__main__':main()
