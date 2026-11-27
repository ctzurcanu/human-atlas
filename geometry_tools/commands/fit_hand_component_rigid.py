"""Bounded rigid hand-artery proposals against the visible skin and fixed bones.

Each source-connected component stays intact. Signed closest-surface residuals
guide a small rigid correction only; full triangles and ports remain separate
acceptance obligations. No runtime installation occurs here.
"""
import argparse
import copy
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from geometry_tools.commands.clear_hand_components import signed_proposal
from geometry_tools.course_io import sha, write_json
from geometry_tools.fit import transform_mesh
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['candidate', 'binding', 'bones', 'skin-reference', 'output']:
        parser.add_argument('--'+key, required=True)
    parser.add_argument('--seed-proposal', help='Hash-bound translation grid proposal on the same input pose')
    parser.add_argument('--network-part', help='Fit one declared part while retaining the whole connected-network report')
    parser.add_argument('--only-components', help='Limit correction to these original component IDs')
    args = parser.parse_args(); started = time.monotonic()
    component_scope = set(map(int, args.only_components.split(','))) if args.only_components else None
    paths = {k: Path(getattr(args,k.replace('-','_'))).resolve()
             for k in ['candidate','binding','bones','skin-reference']}
    data, binding, bones, skin_ref = [json.loads(p.read_text()) for p in paths.values()]
    seed = None
    if args.seed_proposal:
        seed_path=Path(args.seed_proposal).resolve();seed=json.loads(seed_path.read_text())
        if seed['inputs']['candidate']['sha256'] != sha(paths['candidate']):
            raise ValueError('Translation seed belongs to another pose')
        paths['seed-proposal']=seed_path
    if not skin_ref['defaultHiddenExcluded']:
        raise ValueError('Visible skin reference required')
    out = Path(args.output).resolve()
    if out.exists(): raise ValueError('Use a new stage')
    selected_rows = [r for r in data['rows'] if r['partId'] == args.network_part] if args.network_part else data['rows']
    if not selected_rows:
        raise ValueError('Declared network part not found')
    side = selected_rows[0]['side'] if args.network_part else data['side']
    suffix = '.l' if side=='left' else '.r'
    fixed=[]
    for row in bones['rows']:
        if row['sourceObject'].endswith(suffix):
            if sha(row['path']) != row['sha256']: raise ValueError('Fixed bone changed')
            fixed.append(Mesh.load(row['path']))
    triangles=np.concatenate([m.triangles for m in fixed])
    bone=SurfaceIndex(Mesh(triangles.reshape(-1,3),np.arange(triangles.size//3).reshape(-1,3),fixed[0].frame))
    if sha(skin_ref['path']) != skin_ref['sha256']: raise ValueError('Visible skin changed')
    skin=SurfaceIndex(Mesh.load(skin_ref['path']))
    if skin.mesh.frame != bone.mesh.frame: raise ValueError('Context frames differ')
    out.mkdir(); rows=[]
    for row in selected_rows:
        if sha(row['path']) != row['sha256']: raise ValueError('Source candidate changed')
        original=Mesh.load(row['path']); vertices=original.vertices.copy()
        attributes={k:a.copy() for k,a in original.attributes.items()}; corrections=[]
        for part in [r for r in binding['rows'] if r['side']==side and r['partId']==row['partId']]:
            if component_scope is not None and part['physicalComponent'] not in component_scope:
                continue
            if not part['anatomicallyReviewed']: raise ValueError('Unreviewed source ownership')
            face_ids=np.array(part['sourceFaceIds'],int); ids=np.unique(original.faces[face_ids])
            t=original.triangles[face_ids]; all_samples=np.concatenate([original.vertices[ids],t.mean(1),
                (t[:,0]+t[:,1])/2,(t[:,1]+t[:,2])/2,(t[:,2]+t[:,0])/2])
            centre=original.vertices[ids].mean(0)
            def pose(x, points):
                return (points-centre)@Rotation.from_rotvec(x[3:]).as_matrix().T+centre+x[:3]*.001
            def margins(points):
                target_margin=.00025 if args.network_part else .00015
                a=signed_proposal(bone,points)-target_margin
                b=-signed_proposal(skin,points,open_oriented_surface=True)-target_margin
                return a,b
            initial=np.zeros(6)
            if seed:
                seeded_row=next(r for r in seed['rows'] if r['partId']==row['partId'])
                seeded_component=next(c for c in seeded_row['clearanceProposals'] if c['sourceComponent']==part['physicalComponent'])
                initial[:3]=np.array(seeded_component['translation_m'])*1000
            source_margins=margins(all_samples);before=margins(pose(initial,all_samples))
            chosen=[np.linspace(0,len(all_samples)-1,min(160,len(all_samples))).astype(int)]
            # Include actual worst witnesses; an even vertex stride misses short
            # bone/skin penetrations and leaves the optimizer unaware of them.
            for values in [*source_margins,*before]:
                chosen.append(np.argsort(values)[:48])
            samples=all_samples[np.unique(np.concatenate(chosen))]
            def residual(x):
                if time.monotonic()-started>22: raise TimeoutError('22-second native rigid proposal stage')
                a,b=margins(pose(x,samples))
                tr,rr=(.00001,.00001)if args.network_part else(.015,.03)
                return np.r_[np.minimum(a,0)*1000,np.minimum(b,0)*1000,tr*x[:3],rr*x[3:]]
            result=least_squares(residual,initial,bounds=(np.r_[[-6.]*3,[-.35]*3],np.r_[[6.]*3,[.35]*3]),
                                 max_nfev=28,ftol=1e-4,xtol=1e-4,gtol=1e-4,diff_step=1e-3)
            x=result.x; after=margins(pose(x,all_samples))
            before_score=min(float(v.min()) for v in before); after_score=min(float(v.min()) for v in after)
            # A sampled optimizer improvement cannot override the full sample set.
            retained=after_score>before_score+1e-8
            if not retained: x=initial; after=before
            if min(float(v.min())for v in source_margins)>min(float(v.min())for v in after):
                x=np.zeros(6);after=source_margins;retained=False
            matrix=np.eye(4); matrix[:3,:3]=Rotation.from_rotvec(x[3:]).as_matrix()
            matrix[:3,3]=centre-matrix[:3,:3]@centre+x[:3]*.001
            moved=transform_mesh(original,matrix);vertices[ids]=moved.vertices[ids]
            for name in attributes: attributes[name][ids]=moved.attributes[name][ids]
            corrections.append(dict(sourceComponent=part['physicalComponent'],matrix=matrix.tolist(),
                translationAtComponentCentre_m=(x[:3]*.001).tolist(),rotationVectorRadians=x[3:].tolist(),
                evaluations=int(result.nfev),optimizerSuccess=bool(result.success),
                fullSampleWorstMarginBefore_m=before_score,
                fullSampleWorstMarginAfter_m=min(float(v.min()) for v in after),
                boneMargin_m=float(after[0].min()),skinLocalSideMargin_m=float(after[1].min()),
                improvementRetained=retained,fullSamples=len(all_samples),optimizerSamples=len(samples),
                sampledFeasible=all(float(v.min())>=0 for v in after)))
        mesh=Mesh(vertices,original.faces,original.frame,attributes,dict(sourcePartId=row['partId']))
        path=out/f'part-{len(rows):02}.npz';mesh.save(path)
        rows.append(dict(row,path=str(path),sha256=sha(path),geometrySha256=mesh.geometry_sha,
            parentPosePath=row['path'],parentPoseSha256=row['sha256'],rigidClearanceProposals=corrections))
    report=dict(schema='bounded-rigid-hand-component-proposal/v1',side=side,rows=rows,
        inputs={k:dict(path=str(p),sha256=sha(p))for k,p in paths.items()},sourceCodeSha256=sha(__file__),
        skeletonChanged=False,sourceFacesRetained=True,componentShapeAndCalibreRetained=True,
        skinContainmentEstablished=False,installed=False,installationAllowed=False,
        pending=['Full actual triangle checks, source ports, recipient soft tissues, rendering and storage'],
        elapsedSeconds=time.monotonic()-started)
    if args.network_part:
        replacements={r['partId']:r for r in rows}
        merged=copy.deepcopy(data)
        merged.update(rows=[replacements.get(r['partId'],r)for r in data['rows']],
                      rigidClearanceReport={k:v for k,v in report.items()if k!='rows'},
                      elapsedSeconds=report['elapsedSeconds'],sourceCodeSha256=sha(__file__),
                      complete=all(c['sampledFeasible']for r in rows for c in r['rigidClearanceProposals']),
                      installed=False,installationAllowed=False,portsRequireRecoupling=True)
        for pair in merged['pairs']:
            if args.network_part in [pair['firstPart'],pair['secondPart']]:
                pair.update(literalCoordinatesAndEffectiveNormals=False,packedCoordinatesAndEffectiveNormalWordsLiteral=False,requiresRecoupling=True)
        report=merged
    write_json(out/'report.json',report)
    print(json.dumps(dict(side=side,components=sum(len(r['rigidClearanceProposals'])for r in rows),
        sampledFeasible=sum(c['sampledFeasible']for r in rows for c in r['rigidClearanceProposals']),
        elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__': main()
