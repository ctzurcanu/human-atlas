"""Exploratory native-nerve solid classifiers; most native sources fail.

These results are not used to cut meningeal ports. The measurements in
measure-spinal-root-ports.py instead retain literal crossing witnesses.

Caps close classifier cuts at 10 mm from the midline. They are never emitted
as nerve anatomy. Multiple source records of one native nerve are assembled
before exact-coordinate welding. Every clipped classifier is independently
audited; an invalid source cannot authorize a meningeal port.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
import manifold3d as md
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,topology
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.source_face_complement import conform_surface_boundary_edges

ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--output',required=True);ap.add_argument('--begin',type=int,default=0);ap.add_argument('--count',type=int,default=10);args=ap.parse_args()
start=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
data=json.loads((source/'root-mouths.json').read_text());groups=defaultdict(list)
for r in data['roots']:
    side='left'if ':l_'in r['part_id']else'right';groups[(r['level'],side)].append(r)
selected=sorted(groups)[args.begin:args.begin+args.count];rows=[];frame='model:local-female-ta98:assembled-metres'
for level,side in selected:
    records=groups[(level,side)];vv=[];ff=[]
    for r in records:
        z=np.load(r['source_mesh']);v=z['vertices'];f=z['faces'];ff.append(f+sum(len(q)for q in vv));vv.append(v)
    v,inv=np.unique(np.concatenate(vv),axis=0,return_inverse=True);f=inv[np.concatenate(ff)]
    # Cancel only literally coincident opposite faces; retain one instance of
    # a same-oriented duplicate. This is a classifier, not a repaired nerve.
    parity=((f[:,0]>f[:,1]).astype(int)+(f[:,0]>f[:,2])+(f[:,1]>f[:,2]))%2
    keys,which=np.unique(np.sort(f,axis=1),axis=0,return_inverse=True);balance=np.bincount(which,weights=1-2*parity)
    nonzero=balance!=0;f=keys[nonzero].copy();negative=balance[nonzero]<0;f[negative]=f[negative][:,[0,2,1]]
    sign=1 if side=='left'else-1;points=list(v);edges={};cut=.010
    def crossing(a,b):
        key=tuple(sorted((int(a),int(b))))
        if key not in edges:
            p,q=v[list(key)];t=(cut-sign*p[0])/(sign*(q[0]-p[0]));point=p+t*(q-p);point[0]=sign*cut
            edges[key]=len(points);points.append(point)
        return edges[key]
    faces=[]
    for tri in f:
        poly=[]
        for a,b in zip(tri,np.roll(tri,-1)):
            inside_a=sign*v[a,0]<=cut;inside_b=sign*v[b,0]<=cut
            if inside_a:poly.append(int(a))
            if inside_a!=inside_b:poly.append(crossing(a,b))
        for j in range(1,len(poly)-1):faces.append([poly[0],poly[j],poly[j+1]])
    positions=np.asarray(points);faces=np.asarray(faces,int).reshape(-1,3);q=positions[faces]
    faces=faces[np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0]
    all_edges=np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]];keys,which,counts=np.unique(np.sort(all_edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    boundary=all_edges[counts[which]==1];nxt=defaultdict(list)
    for a,b in boundary:nxt[int(a)].append(int(b))
    branched=any(len(p)!=1 for p in nxt.values())or set(nxt)!={int(x)for x in boundary[:,1]};caps=[];cap_records=[]
    if not branched:
        todo=set(nxt)
        while todo:
            first=min(todo);loop=[];i=first
            while i in todo:todo.remove(i);loop.append(i);i=nxt[i][0]
            assert i==first
            ids=np.asarray(loop,int);p=positions[ids];center=p.mean(0);_,_,axes=np.linalg.svd(p-center,full_matrices=False);xy=(p-center)@axes[:2].T
            desired=-np.cross(p-center,np.roll(p,-1,axis=0)-center).sum(0)
            if np.linalg.norm(desired)<=1e-12*np.max(np.ptp(p,axis=0))**2:
                # A collinear source-boundary visit needs edge conformity,
                # not a zero-area disk that merely hides an open-edge count.
                continue
            cap=ids[np.asarray(md.triangulate([xy*1e6],epsilon=0.,allow_convex=False),int)]
            q=positions[cap];flip=np.einsum('ij,j->i',np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),desired)<0;cap[flip]=cap[flip][:,[0,2,1]]
            caps.extend(cap);cap_records.append(dict(center_m=center.tolist(),vertices=len(ids),plane_residual_m=float(np.abs((p-center)@axes[2]).max()),
                on_declared_classifier_cut=bool(np.max(np.abs(sign*p[:,0]-cut))<1e-12),not_an_anatomical_cap=True))
    surfaces={'native':faces,'classifier-caps':np.asarray(caps,int).reshape(-1,3)}
    surfaces={k:f for k,f in surfaces.items()if len(f)};cells={'mask':[dict(surface=s,side=1,neighbor=None)for s in surfaces]}
    surfaces,conformance=conform_surface_boundary_edges(positions,surfaces,cells,tolerance=2e-12)
    faces=np.concatenate(list(surfaces.values()));edge_owners=defaultdict(list)
    for i,tri in enumerate(faces):
        for a,b in zip(tri,np.roll(tri,-1)):
            edge_owners[tuple(sorted((int(a),int(b))))].append((i,int(a)<int(b)))
    adjacency=defaultdict(list)
    for owners in edge_owners.values():
        if len(owners)!=2:continue
        (i,a),(j,b)=owners;adjacency[i].append((j,a==b));adjacency[j].append((i,a==b))
    flips={};components=[];orientation_conflicts=0
    for first in range(len(faces)):
        if first in flips:continue
        flips[first]=False;pending=[first];component=[]
        while pending:
            i=pending.pop();component.append(i)
            for j,relative in adjacency[i]:
                value=flips[i]^relative
                if j in flips:orientation_conflicts+=int(flips[j]!=value);continue
                flips[j]=value;pending.append(j)
        components.append(component)
    negative=np.array([flips[i]for i in range(len(faces))]);faces[negative]=faces[negative][:,[0,2,1]]
    for component in components:
        q=positions[faces[component]]
        if np.einsum('ij,ij->i',q[:,0],np.cross(q[:,1],q[:,2])).sum()<0:faces[component]=faces[component][:,[0,2,1]]
    mesh=Mesh(positions,faces,frame);canonicalization=None
    if topology(mesh)['closed_oriented_manifold']:
        solid,_=make_solid(mesh.vertices,mesh.faces);simplified=solid.simplify(2e-7);encoded=simplified.to_mesh64()
        volume_before=float(solid.volume());volume_after=float(simplified.volume());change=abs(volume_after-volume_before)/max(abs(volume_before),1e-20)
        assert change<=.002,'Classifier simplification exceeded the explicit 0.2% volume budget'
        mesh=Mesh(encoded.vert_properties[:,:3],encoded.tri_verts,frame)
        canonicalization=dict(surface_tolerance_m=2e-7,relative_volume_change=change,source_nerve_not_modified=True)
    proof=audit(mesh,require_closed=True,max_pairs=200000,max_tests=2000000)
    path=out/(level+'-'+side+'.npz');mesh.save(path)
    row=dict(level=level,side=side,source_part_ids=[r['part_id']for r in records],source_geometry_sha256=[r['source_geometry_sha256']for r in records],
        mask=str(path),mask_sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,passed=proof['passed'],caps=cap_records,
        branched_boundary_before_caps=branched,conformance=conformance,classifier_orientation_conflicts=orientation_conflicts,classifier_canonicalization=canonicalization,proof=proof,not_installed_anatomy=True)
    atomic_json(path.with_suffix('.json'),row);rows.append(row)
    print(json.dumps(dict(level=level,side=side,passed=row['passed'],triangles=len(mesh.faces),failures=proof['failures'],elapsed_seconds=time.monotonic()-start)),flush=True)
atomic_json(out/'masks.json',dict(source_report_sha256=digest_file(source/'root-mouths.json'),roots=rows,
    classification_only=True,anatomy_changed=False,skeleton_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
