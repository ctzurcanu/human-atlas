"""Regenerate six typed controls on an explicit exposed mucosal cell.

Writes reference metadata only, never geometry or a catalogue. Basal/stem
facets removed by the declared floor-of-mouth reflection are explicitly
reported, and the inferior control is deliberately redefined as the genuine
exposed ventral blade. All regional boundaries remain approximate.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from tongue_whole_input import whole_input
from tongue_basal_contract import MUCO,bind_features,points_on_triangles,require

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'


def read(path):
    return json.loads(Path(path).read_text())


def main(patch,previous,output,carrier_manifest=None,carrier_mesh=None):
    began=time.monotonic()
    def remaining():
        if time.monotonic()-began>15:
            raise TimeoutError('Exposed feature metadata exceeded its 15 second worker limit')
    patch,previous,output=patch.resolve(),previous.resolve(),output.resolve()
    if patch.is_dir():
        patch/='candidate.json'
    require(patch.is_relative_to(WORK) and previous.is_relative_to(WORK), 'Use explicit local tongue inputs')
    require(output.is_relative_to(WORK) and output!=WORK and not output.exists(), 'Use a new metadata stage')
    candidate,old=read(patch),read(previous)
    source_path=Path(candidate['source'])
    source,_,binding=whole_input(source_path,Path(old['wholeInput']['fit_result_path']))
    require(digest_file(source_path)==candidate['sourceSha256']==old['carrierMeshSha256']
            and source.geometry_sha==candidate['sourceGeometrySha256'], 'Actual corrected whole source changed')
    require(old.get('schema')=='tongue-construction-features/v1'
            and old.get('wholeInput')==binding, 'Previous features must bind this reviewed whole')
    raw_carrier_path=Path(candidate['carrier'])
    raw_carrier=Mesh.load(raw_carrier_path)
    provenance_path=patch.with_name('source-face-provenance.json')
    provenance=read(provenance_path)
    require(provenance.get('schema')=='tongue-exposed-carrier-ancestry/v1'
            and provenance['sourceGeometrySha256']==source.geometry_sha
            and provenance['carrierGeometrySha256']==raw_carrier.geometry_sha,
            'Exposed patch/source-face provenance is stale')
    require(bool(carrier_manifest)==bool(carrier_mesh), 'Override both actual carrier manifest and mesh together')
    manifest_path=Path(carrier_manifest or candidate['packedManifest']).resolve()
    manifest=read(manifest_path)
    actual_path=Path(carrier_mesh or manifest_path.with_name('TA98TONGUE_male_mucosa.npz')).resolve()
    actual=Mesh.load(actual_path)
    require(manifest_path.is_relative_to(WORK) and actual_path.is_relative_to(WORK)
            and manifest.get('schema')=='shared-finite-tissues/v1'
            and manifest.get('passed') is True
            and digest_file(manifest['bank'])==manifest['bank_sha256'],
            'Use an actual finite-cell bank and matching mucosal export')
    require(actual.geometry_sha==manifest['checks']['cells']['mucosa']['geometry_sha256'],
            'The actual mucosal export does not match its emitted-bank audit')
    count=len(raw_carrier.faces)
    require(count<=len(actual.faces) and count>0, 'Actual bank has no exposed outer carrier')
    triangles=actual.triangles[:count]
    require(np.linalg.norm(triangles-raw_carrier.triangles,axis=2).max()<=2e-8,
            'Actual emitted outer facet order/ancestry changed')
    parents={r['carrierFace']:r['sourceFace'] for r in provenance['faces']+provenance['cutFaces']}
    require(set(parents)==set(range(count)), 'Every exposed outer triangle needs actual source provenance')
    ancestry_rows=[]
    maximum_residual=0.
    for face in range(count):
        remaining()
        parent=parents[face]
        original=source.triangles[parent]
        row=dict(sourceFace=parent,carrierFace=face)
        orders=[[0,1,2],[1,2,0],[2,0,1]]
        literal=next((order for order in orders if np.array_equal(original[order],triangles[face])),None)
        if literal is not None:
            row['cornerOrder']=literal
        else:
            edges=np.stack((original[1]-original[0],original[2]-original[0]),axis=1)
            uv=np.linalg.lstsq(edges,(triangles[face]-original[0]).T,rcond=None)[0].T
            bary=np.column_stack((1-uv.sum(1),uv))
            require(np.min(bary)>=-2e-5 and np.max(bary)<=1+2e-5,
                    'Cut emitted vertex leaves its declared source triangle')
            bary=np.clip(bary,0,1);bary/=bary.sum(1)[:,None]
            error=float(np.linalg.norm(bary@original-triangles[face],axis=1).max())
            require(error<=2e-8,'Actual cut-face provenance exceeds its position encoding bound')
            maximum_residual=max(maximum_residual,error)
            row.update(sourceBarycentric=bary.tolist(),maximumSourcePointErrorM=error)
        ancestry_rows.append(row)
    ancestry=dict(schema='tongue-covering-face-ancestry/v1',carrierPartId=MUCO,
        carrierBankSha256=manifest['bank_sha256'],sourceGeometrySha256=source.geometry_sha,
        carrierGeometrySha256=actual.geometry_sha,faces=ancestry_rows,
        actualCarrierMesh=str(actual_path),actualCarrierMeshSha256=digest_file(actual_path),
        wholeInput=binding,exposureBoundary=str(patch.with_name('reflection-boundary.json')),
        maximumCutPointReconstructionErrorM=maximum_residual,anatomicalAcceptance=False)
    old_groups={g['key']:g for g in old['termGroups']}
    parent_array=np.array([parents[i] for i in range(count)])
    centers=triangles.mean(1)
    geometric=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    geometric/=np.linalg.norm(geometric,axis=1)[:,None]
    groups=[];decisions=[]
    present=set(parent_array.tolist())
    def anchor(face,weights=None):
        weights=np.asarray(weights if weights is not None else [1/3]*3,float)
        return dict(partId=MUCO,face=int(face),barycentric=weights.tolist(),point=(weights@triangles[face]).tolist())
    for key in ('dorsum','anterior-dorsum','posterior-dorsum','inferior-surface'):
        original=old_groups[key]
        old_faces=set(original['surface']['faces'])
        if key=='inferior-surface':
            # Deliberate redefinition: only the actual exposed free ventral
            # tongue surface remains; the computational inferior stem is not
            # a mucosal surface. Source anterior identity stays approximate.
            keep=np.isin(parent_array,list(old_faces))&(geometric[:,1]<-.25)
            decision=('Redefined as true exposed ventral body/free blade above the declared floor-of-mouth '
                      'reflection. Computational basal/stem faces are excluded explicitly. Exact frenulum '
                      'and individual reflection anatomy remain unestablished.')
        else:
            keep=np.isin(parent_array,list(old_faces))
            decision=('Preserves exposed descendant facets of the prior superior/posterior region. '
                      'Prior computational basal posterior faces below the declared reflection are '
                      'explicitly removed from this reference; no donor-specific dorsal segmentation is claimed.')
        faces=np.flatnonzero(keep)
        require(len(faces)>0,'No genuine exposed facets remain for '+key)
        group=copy.deepcopy(original)
        group['surface']=dict(partId=MUCO,faces=faces.tolist(),sourceMeshSha256=digest_file(actual_path))
        centre=centers[faces].mean(0)
        selected=int(faces[np.argmin(np.linalg.norm(centers[faces]-centre,axis=1))])
        group['anchor']=anchor(selected)
        group['coverageLimitation']=decision+' Regional limits are approximate; this is a typed surface reference, not additional tissue.'
        group['exposureDefinitionRevised']=True
        groups.append(group)
        missing=sorted(old_faces-present)
        old_anchor=np.asarray(original['anchor']['point'])
        try:
            points_on_triangles([old_anchor],triangles)
            anchor_still_exposed=True
        except ValueError:
            anchor_still_exposed=False
        decisions.append(dict(key=key,previousSourceFaces=len(old_faces),actualExposedFaces=len(faces),
            excludedSourceFaces=missing,cutSourceFaces=sorted(old_faces&{r['sourceFace'] for r in provenance['cutFaces']}),
            previousAnchorSourceFace=original['anchor']['face'],previousAnchorStillExposed=anchor_still_exposed,
            replacementAnchorFace=selected,decision=decision))
    sets=[set(g['surface']['faces']) for g in groups[:3]]
    require(sets[0]==sets[1]|sets[2] and not sets[1]&sets[2],
            'Anterior/posterior references must partition the exposed dorsal reference exactly')
    tip=source.vertices[np.argmax(source.vertices[:,2])]
    matches=np.argwhere(np.all(triangles==tip,axis=2))
    require(len(matches)>0,'The actual native apex is missing from the exposed carrier')
    face,corner=map(int,matches[0])
    apex=copy.deepcopy(old_groups['apex']);apex['anchor']=anchor(face,np.eye(3)[corner])
    apex['coverageLimitation']='Literal protected anterior-most native vertex, retained exactly on the emitted exposed mucosa; point reference only.'
    groups.append(apex)
    margin=copy.deepcopy(old_groups['margin'])
    margin['anchor']=copy.deepcopy(apex['anchor'])
    for curve in margin['curves']:
        points_on_triangles(curve['points'],triangles)
        curve.pop('sourceBankSha256',None)
        curve['sourceMeshSha256']=digest_file(actual_path)
    margin['coverageLimitation']='Prior lateral-margin reference points all remain on actual exposed outer facets and are retained literally. Curves are approximate anatomical references, not tubes or additional solids.'
    groups.append(margin)
    features=dict(schema='tongue-construction-features/v1',featureFaceSpace='carrier',
        exposureRegenerated=True,carrierGeometrySha256=actual.geometry_sha,
        sourceCatalogueSha256=old['sourceCatalogueSha256'],inventorySha256=old['inventorySha256'],
        hostPartId=MUCO,carrierBankHashes={},termGroups=groups,wholeInput=binding,
        sourceWholeMesh=str(source_path),sourceWholeMeshSha256=digest_file(source_path),
        carrierMesh=str(actual_path),carrierMeshSha256=digest_file(actual_path),carrierBankBindingPending=True,
        sourceGeometryChanged=False,additionalSolidMeshes=0,additionalImageAssets=0,
        sourceFacets=len(source.faces),exposedCarrierFacets=count,terms=6,anatomicalAcceptance=False,
        references=old['references']+['https://pmc.ncbi.nlm.nih.gov/articles/PMC6850428/'])
    bound=bind_features(features,ancestry,source,actual,manifest['bank_sha256'])
    require(len(bound)==6,'Actual runtime-compatible feature binding failed')
    audit=dict(schema='tongue-exposed-controls-audit/v1',source=str(source_path),sourceSha256=digest_file(source_path),
        patch=str(patch),patchSha256=digest_file(patch),previousControls=str(previous),previousControlsSha256=digest_file(previous),
        actualCarrier=str(actual_path),actualCarrierSha256=digest_file(actual_path),bank=str(manifest['bank']),bankSha256=manifest['bank_sha256'],
        decisions=decisions,apexLiteral=True,marginCurvesRetainedLiteral=True,
        marginPoints=sum(len(c['points'])for c in margin['curves']),actualAnchorsAndCurveAttachmentPassed=True,
        literalAncestorFacets=sum('cornerOrder'in r for r in ancestry_rows),cutAncestorFacets=sum('sourceBarycentric'in r for r in ancestry_rows),
        maximumCutPointReconstructionErrorM=maximum_residual,geometryModified=False,anatomicalAcceptance=False,
        installed=False,coveringMuscleInterfacesPassed=False,elapsedSeconds=time.monotonic()-began)
    remaining()
    output.mkdir(parents=True)
    atomic_json(output/'surface-controls.json',features)
    atomic_json(output/'face-ancestry.json',ancestry)
    atomic_json(output/'controls-audit.json',audit)
    print(json.dumps(dict(output=str(output),terms=6,additionalMeshes=0,excludedSourceFaces={r['key']:len(r['excludedSourceFaces'])for r in decisions},
        apexLiteral=True,marginPointsRetained=audit['marginPoints'],elapsedSeconds=audit['elapsedSeconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--patch',type=Path,required=True)
    parser.add_argument('--features',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--carrier-manifest',type=Path)
    parser.add_argument('--carrier-mesh',type=Path)
    args=parser.parse_args()
    main(args.patch,args.features,args.output,args.carrier_manifest,args.carrier_mesh)
