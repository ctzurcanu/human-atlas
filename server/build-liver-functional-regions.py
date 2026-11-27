"""Save a conforming calculation mesh and native functional guide fields.

Construction stages are separately bounded; no atlas is modified here.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.harmonic_wall import tetrahedralize
from geometry_tools.open3d_backend import library
from geometry_tools.algorithms.labelled_tetrahedra import partition_labelled_tetrahedra

ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['field','partition','regularize']);ap.add_argument('--version',default='v144');args=ap.parse_args()
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion'
out=base/('functional-territories-'+args.version);out.mkdir(exist_ok=True)
if args.stage=='field':
    archive=base/'closed-liver-cells-v102/arrangement.npz';z=np.load(archive)
    manifest=json.loads((archive.parent/'assembly.json').read_text());v=z['positions']
    # Remove only the reciprocal internal gross-lobe interfaces. The actual
    # capsule-facing skin and computational hilar closures remain literal.
    owners=defaultdict(list)
    for name in ['right','left','quadrate','caudate']:
        for ref in manifest['cells'][name]:owners[ref['surface']].append(ref['side'])
    ff=[]
    for sid,sides in owners.items():
        if len(sides)==2:
            assert sorted(sides)==[-1,1];continue
        assert len(sides)==1
        f=z[manifest['surfaces'][sid]];ff.append(f if sides[0]==1 else f[:,[0,2,1]])
    f=np.concatenate(ff);used,inv=np.unique(f,return_inverse=True);v=v[used];f=inv.reshape(-1,3)
    bulk=Mesh(v,f,'model:local-ta98:assembled-metres');q=audit(bulk,require_closed=True,max_pairs=1000000,max_tests=8000000)
    atomic_json(out/'bulk-audit.json',q);assert q['passed'],q['failures'];bulk.save(out/'bulk.npz')
    p,t,*_,tet_evidence=tetrahedralize(v,f)
    atomic_json(out/'tetrahedral-evidence.json',tet_evidence)
    labels=['VI','V','III','IV','II','VII','VIII','I'];sources=[0,1,2,3,4,7,8,9];scores=[];bindings={}
    o3d=library()
    for label,ordinal in zip(labels,sources):
        path=base/f'common-field-v13/{ordinal:03d}.npz';m=Mesh.load(path)
        scene=o3d.t.geometry.RaycastingScene(nthreads=1)
        scene.add_triangles(o3d.core.Tensor(m.vertices.astype(np.float32)),o3d.core.Tensor(m.faces.astype(np.uint32)))
        hit=scene.compute_closest_points(o3d.core.Tensor(p.astype(np.float32)),nthreads=1)
        closest=hit['points'].numpy().astype(float);normal=hit['primitive_normals'].numpy().astype(float)
        delta=p-closest;distance=np.linalg.norm(delta,axis=1)
        signed=distance*np.where(np.einsum('ij,ij->i',delta,normal)>=0,1.,-1.)
        scores.append(-signed);bindings[label]=dict(file=str(path),sha256=digest_file(path),geometry_sha256=m.geometry_sha)
    atomic_npz(out/'field.npz',positions=p,tetrahedra=t,scores=np.array(scores).T,labels=np.array(labels))
    atomic_json(out/'field-evidence.json',dict(source_bulk_sha256=digest_file(out/'bulk.npz'),guides=bindings,
        method='Closest oriented native segment surface distance, piecewise affine on a boundary-conforming calculation complex',
        boundaries_provisional=True,source_exterior_literal=True,calculation_tetrahedra=len(t),elapsed_seconds=time.monotonic()-start))
    print(json.dumps(dict(stage=args.stage,tetrahedra=len(t),points=len(p),elapsed_seconds=time.monotonic()-start)),flush=True)
elif args.stage=='regularize':
    from scipy.sparse import coo_matrix,diags
    from scipy.sparse.linalg import spsolve
    from geometry_tools.audit import winding_contains
    src=base/'functional-territories-v144';z=np.load(src/'field.npz');p=z['positions'];t=z['tetrahedra'];scores=z['scores'].copy()
    edges=np.unique(np.sort(np.concatenate([t[:,[i,j]]for i in range(4)for j in range(i+1,4)]),axis=1),axis=0)
    weights=1/np.maximum(np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1),1e-8)
    W=coo_matrix((np.r_[weights,weights],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(p),len(p))).tocsr()
    L=diags(np.asarray(W.sum(1)).ravel())-W
    # One documented 6 mm diffusion of the supplied source fields suppresses
    # tiny disconnected closest-face sign islands; no exterior control moves.
    mass=np.asarray(W.sum(1)).ravel();A=diags(mass)+.35*L
    scores[:,:7]=spsolve(A,mass[:,None]*scores[:,:7])
    caudate=Mesh.load(base/'closed-liver-cells-v102/caudate.npz');o3d=library();scene=o3d.t.geometry.RaycastingScene(nthreads=1)
    scene.add_triangles(o3d.core.Tensor(caudate.vertices.astype(np.float32)),o3d.core.Tensor(caudate.faces.astype(np.uint32)))
    signed=scene.compute_signed_distance(o3d.core.Tensor(p.astype(np.float32)),nthreads=1,nsamples=3).numpy().astype(float)
    # The validated gross caudate is segment I. The remaining seven fields
    # compete only outside its own source-boundary distance field.
    scores[:,:7]-=scores[:,:7].max(1)[:,None];scores[:,:7]+=signed[:,None];scores[:,7]=-signed
    atomic_npz(out/'field.npz',positions=p,tetrahedra=t,scores=scores,labels=z['labels'])
    atomic_json(out/'field-evidence.json',dict(source=str(src/'field.npz'),sourceSha256=digest_file(src/'field.npz'),
        caudateGeometrySha256=caudate.geometry_sha,smoothing='One normalized graph diffusion, coefficient 0.35, of the seven native segment distance fields',
        exteriorLiteral=True,boundariesProvisional=True,elapsedSeconds=time.monotonic()-start))
    print(json.dumps(dict(stage=args.stage,elapsedSeconds=time.monotonic()-start)),flush=True)
else:
    z=np.load(out/'field.npz');result=partition_labelled_tetrahedra(z['positions'],z['tetrahedra'],z['scores'],seconds=20,maximum_tetrahedra=160000,
        identical_score_priority=[7,0,1,2,3,4,5,6] if args.version!='v144' else None)
    arrays={'positions':result['positions']};labels=z['labels'].tolist()
    for i,label in enumerate(labels):arrays[label]=result['cells'][i]
    atomic_npz(out/'partition.npz',**arrays)
    atomic_json(out/'partition-evidence.json',dict(evidence=result['evidence'],labels=labels,elapsed_seconds=time.monotonic()-start,installed=False))
    print(json.dumps(dict(stage=args.stage,triangles={label:len(result['cells'][i])for i,label in enumerate(labels)},elapsed_seconds=time.monotonic()-start)),flush=True)
