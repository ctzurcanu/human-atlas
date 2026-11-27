"""Finite anatomically directed translation proposals for intact hand vessels.

Uses an explicit reviewed component declaration, fixed whole-bone envelopes and
the actual outer skin. Native point samples only propose a candidate; actual
triangles, source interfaces and anatomical review must still pass.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.open3d_backend import SurfaceIndex, library
from geometry_tools.course_io import sha, write_json


def signed_proposal(index, points, *, open_oriented_surface=False):
    nearest, _ = index.closest_points(points)
    if open_oriented_surface:
        # A local side-of-surface proposal, never a closed-volume certificate.
        t = index.mesh.triangles[nearest['face_ids']]
        normal = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
        length = np.linalg.norm(normal, axis=1)
        if np.any(length == 0):
            raise ValueError('Degenerate visible skin reference')
        dot = np.einsum('ij,ij->i', points-nearest['points'], normal/length[:, None])
        return nearest['distance_m']*np.where(dot < 0, -1., 1.)
    p = ((points-index.origin)/index.scale).astype(np.float32)
    inside = index.scene.compute_occupancy(library().core.Tensor(p), nthreads=2, nsamples=3).numpy() > .5
    return nearest['distance_m']*np.where(inside, -1., 1.)


def unit(x):
    length = np.linalg.norm(x)
    if length < 1e-10:
        raise ValueError('No anatomical translation direction')
    return x/length


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['candidate', 'binding', 'guides', 'bones', 'catalogue', 'output']:
        p.add_argument('--'+key, required=True)
    p.add_argument('--skin')
    p.add_argument('--skin-reference', help='Hash-bound visible surface diagnostic; open surfaces use local normal proposals only')
    args = p.parse_args(); start = time.monotonic()
    if bool(args.skin) == bool(args.skin_reference):
        p.error('Declare exactly one visible --skin or --skin-reference')
    paths = {k: Path(getattr(args, k)).resolve() for k in ['candidate','binding','guides','bones','catalogue']}
    data, declaration, guides, bones, atlas = [json.loads(q.read_text()) for q in paths.values()]
    side = data['side']; suffix = '.l' if side == 'left' else '.r'
    if bones['catalogueSha256'] != sha(paths['catalogue']):
        raise ValueError('Fixed bone envelopes are bound to another atlas revision')
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new proposal stage')
    out.mkdir()
    meshes = []
    for r in bones['rows']:
        if r['sourceObject'].endswith(suffix):
            if sha(r['path']) != r['sha256']:
                raise ValueError('Fixed bone changed')
            meshes.append(Mesh.load(r['path']))
    triangles = np.concatenate([m.triangles for m in meshes])
    bone_index = SurfaceIndex(Mesh(triangles.reshape(-1,3),np.arange(triangles.size//3).reshape(-1,3),meshes[0].frame))
    if args.skin_reference:
        reference_path = Path(args.skin_reference).resolve()
        reference = json.loads(reference_path.read_text())
        if reference['catalogueSha256'] != sha(paths['catalogue']) or sha(reference['path']) != reference['sha256']:
            raise ValueError('Visible skin source changed')
        if not reference.get('defaultHiddenExcluded'):
            raise ValueError('Skin reference must exclude hidden source variants')
        skin = Mesh.load(reference['path'])
        skin_source = dict(path=str(reference_path), sha256=sha(reference_path))
    else:
        if args.skin in atlas['defaultHidden']:
            raise ValueError('Hidden source skin cannot constrain the displayed hand')
        skin_part = next(r for r in atlas['parts'] if r['id'] == args.skin)
        skin = AtlasReader().read_part(atlas,skin_part,frame=meshes[0].frame,catalogue_path=paths['catalogue'])
        skin_source = dict(partId=args.skin)
    skin_index = SurfaceIndex(skin)
    _, aliases = np.unique(skin.vertices,axis=0,return_inverse=True)
    sf = aliases[skin.faces]; edges=np.sort(sf[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1)
    _, counts = np.unique(edges,axis=0,return_counts=True)
    open_skin = bool(np.any(counts != 2))
    if open_skin and not args.skin_reference:
        raise ValueError('Outer skin is not closed; do not use occupancy as a compartment proposal')
    groups = [r for r in guides['rows'] if r['side'] == side]
    palm = np.array([r['namedBoneCorrespondences'][0]['targetCentre_m'] for r in groups])
    _,_,vh=np.linalg.svd(palm-palm.mean(0)); normal=vh[-1]; normal *= 1 if normal[2]>=0 else -1
    rows=[]
    for row in data['rows']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Source pose changed')
        original=Mesh.load(row['path']); vertices=original.vertices.copy(); proposals=[]
        bindings=[r for r in declaration['rows'] if r['side']==side and r['partId']==row['partId']]
        for binding in bindings:
            ids=np.array(binding['rawVertexIds'],int); f=np.array(binding['sourceFaceIds'],int)
            t=original.triangles[f]; centre=t.mean((0,1)); numbers=binding['proposedCourses']
            sign=-1 if 'dorsal_digital' in row['partId'] else 1
            outward=normal*sign
            if len(numbers)==1:
                group=next(r for r in groups if r['digit']==numbers[0])
                line=np.array([r['targetCentre_m'] for r in group['namedBoneCorrespondences']])
                line=np.r_[line[:1]-(line[1:2]-line[:1]),line,line[-1:]+(line[-1:]-line[-2:-1])]
                a=line[:-1]; e=np.diff(line,axis=0)
                fraction=np.clip(np.einsum('ij,ij->i',centre-a,e)/np.sum(e*e,axis=1),0,1)
                nearest=a+fraction[:,None]*e; k=np.argmin(np.linalg.norm(nearest-centre,axis=1))
                tangent=unit(e[k]); away=centre-nearest[k];away-=tangent*(away@tangent)
                if np.linalg.norm(away)>1e-8:outward=unit(away)
            directions=[unit(outward),unit(normal*sign)]
            combined=outward+normal*sign
            if np.linalg.norm(combined)>1e-8:directions.append(unit(combined))
            offsets=np.array([d*s for s in np.linspace(0,.006,25) for d in directions])
            samples=np.concatenate([original.vertices[ids],t.mean(1),((t[:,0]+t[:,1])/2),((t[:,1]+t[:,2])/2),((t[:,2]+t[:,0])/2)])
            query=(samples[None]+offsets[:,None]).reshape(-1,3)
            bone=signed_proposal(bone_index,query).reshape(len(offsets),len(samples))
            skin_distance=signed_proposal(skin_index,query,open_oriented_surface=open_skin).reshape(len(offsets),len(samples))
            bone_margin=bone.min(1)-.0001;skin_margin=(-skin_distance).min(1)-.0001
            score=np.minimum(bone_margin,skin_margin);feasible=np.flatnonzero(score>=0)
            if len(feasible):
                selected=int(feasible[np.argmin(np.linalg.norm(offsets[feasible],axis=1))])
            else:selected=int(np.argmax(score))
            delta=offsets[selected];vertices[ids]+=delta
            proposals.append(dict(sourceComponent=binding['physicalComponent'],translation_m=delta.tolist(),
                sampledFeasible=bool(len(feasible)),boneMargin_m=float(bone_margin[selected]),skinMargin_m=float(skin_margin[selected]),
                originalWorstSampleMargin_m=float(score[0]),selectedWorstSampleMargin_m=float(score[selected]),
                finiteGridCandidates=len(offsets),samplesPerComponent=len(samples),nativeProposalOnly=True))
            if time.monotonic()-start>18:raise TimeoutError('18-second native component clearance stage')
        candidate=Mesh(vertices,original.faces,original.frame,original.attributes,dict(sourcePartId=row['partId']))
        path=out/f'part-{len(rows):02}.npz';candidate.save(path)
        rows.append(dict(row,path=str(path),sha256=sha(path),geometrySha256=candidate.geometry_sha,
                         parentPosePath=row['path'],parentPoseSha256=row['sha256'],clearanceProposals=proposals))
    report=dict(schema='intact-hand-component-clearance-proposal/v1',side=side,rows=rows,
        inputs={k:dict(path=str(q),sha256=sha(q))for k,q in paths.items()},skinSource=skin_source,skinGeometrySha256=skin.geometry_sha,
        skinProposalMethod='nearest oriented visible surface; not volume containment' if open_skin else 'closed surface occupancy',
        uniformSourceComponentShapesRetained=True,sourceFacesAndNormalWordsRetained=True,skeletonChanged=False,
        allSampledComponentsFeasible=all(c['sampledFeasible']for r in rows for c in r['clearanceProposals']),
        installed=False,installationAllowed=False,actualTriangleAndPortChecksPending=True,
        elapsedSeconds=time.monotonic()-start)
    write_json(out/'report.json',report)
    print(json.dumps(dict(side=side,components=sum(len(r['clearanceProposals'])for r in rows),
        sampledFeasible=sum(c['sampledFeasible']for r in rows for c in r['clearanceProposals']),elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
