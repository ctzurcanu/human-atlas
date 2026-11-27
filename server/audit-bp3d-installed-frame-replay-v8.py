"""Replay the recorded whole-body height conversion; never estimate a fit."""
from pathlib import Path
import sys,time,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader

SOURCE=ROOT/'public/models/atlas.json';TARGET=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1/bp3d-installed-frame-replay-v8'

def main():
    began=time.monotonic();assert digest_file(TARGET)==EXPECTED;OUT.mkdir(exist_ok=False)
    a=json.loads(SOURCE.read_text());b=json.loads(TARGET.read_text());ap={p['id']:p for p in a['parts']};bp={p['id']:p for p in b['parts']};height=b['heightAdjustment'];s=float(height['scale']);t=np.array(height['translationMeters'],float)
    assert any(q['file']=='public/models/atlas.json'and q['sha256']==digest_file(SOURCE)for q in b['provenance']['sourceCatalogues'])
    reader=AtlasReader();rows=[]
    for pid in ['FJ1723','FJ1723M','FJ1679','FJ1679M','FJ1758','FJ1806','FJ1755','FJ1803','FJ1771','FJ1818','FJ3394']:
        if pid not in bp:continue
        x=reader.read_part(a,ap[pid],frame='source:BodyParts3D-4.0:packed-native-metres-y-up',catalogue_path=SOURCE);y=reader.read_part(b,bp[pid],frame='model:local-ta98:assembled-metres',catalogue_path=TARGET)
        samefaces=np.array_equal(x.faces,y.faces);samecounts=x.vertices.shape==y.vertices.shape
        exact=bool(samecounts and np.array_equal((x.vertices*s+t).astype('<f4').astype(float),y.vertices));error=float(np.linalg.norm((x.vertices*s+t)-y.vertices,axis=1).max())if samecounts else None
        normals=bool(np.array_equal(x.attributes['normals'],y.attributes['normals']))
        rows.append(dict(id=pid,originalName=ap[pid]['name'],targetName=bp[pid]['name'],sameActualSourceFaceIndices=samefaces,sameSourceControlCount=samecounts,recordedTransformFloat32ReplayExact=exact,maximumUnquantizedReplayError_m=error,packedNormalWordsLiteral=normals,sourceGeometrySha256=x.geometry_sha,targetGeometrySha256=y.geometry_sha,sourceBankSha256=x.source['bank_sha256'],targetBankSha256=y.source['bank_sha256']))
    assert rows and all(r['sameActualSourceFaceIndices']and r['recordedTransformFloat32ReplayExact']and r['packedNormalWordsLiteral']for r in rows)
    before=reader.asset(height['previousCatalogue']);beforeproof=dict(path=str(before),exists=before.exists(),expectedSha256=height['previousAtlasSha256'])
    if before.exists():beforeproof.update(actualSha256=digest_file(before),matchesRecordedDigest=digest_file(before)==height['previousAtlasSha256'])
    out=dict(schema='recorded-bp3d-height-frame-replay/v8',sourceCatalogue=str(SOURCE),sourceCatalogueSha256=digest_file(SOURCE),currentCatalogueSha256=EXPECTED,existingRecordedUniformScale=s,existingRecordedTranslation_m=t.tolist(),transformDefinition='Only replay of the already recorded male 183 cm whole-assembly heightAdjustment. Each native source control is converted with Float32(sourceFloat32*scale+translation). No fitted parameters or new registration.',heightAdjustment=height,transformImplementationSha256=digest_file(ROOT/'scripts/resize-ta98-male-height.py'),recordedPreHeightCatalogue=beforeproof,rows=rows,passed=True,newRegistrationEstimated=False,anatomicalAlignmentToZaBrainOrSkeletonCertified=False,noSourceOrLiveMetadataChanged=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',out);assert digest_file(TARGET)==EXPECTED
    print(json.dumps(dict(output=str(OUT/'proof.json'),replayedParts=len(rows),allFloat32ControlsAndNormalsLiteral=True,scale=s,translation_m=t.tolist(),elapsedSeconds=out['elapsedSeconds'])))
if __name__=='__main__':main()
