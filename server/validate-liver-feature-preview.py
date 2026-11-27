from pathlib import Path
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader

folder=ROOT/'.local-models/geometry-workbench/liver-completion/liver-surface-controls-v143'
path=folder/'with-landmark.json';atlas=json.loads(path.read_bytes())
ordinary=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_bytes())
assert atlas['parts']==ordinary['parts'] and atlas['chunks']==ordinary['chunks']
marker=atlas['liverFeatureIntegration'];parts={p['id']:p for p in atlas['parts']}
host=parts[marker['hostPartId']];reader=AtlasReader();mesh=reader.read_part(atlas,host,frame='model:local-ta98:assembled-metres')
assert marker['bankSha256']==mesh.source['bank_sha256']
plan=json.loads((folder/'construction.json').read_text())
assert digest_file(folder/'construction.json')==marker['sourcePlanSha256']
assert mesh.geometry_sha==plan['hostGeometrySha256'] and mesh.attributes_sha==plan['hostAttributesSha256']
outer_count=host['tissueDrawGroups'][0]['count']//3
checks=[]
for g in marker['termGroups']+atlas['liverLandmarkIntegration']['termGroups']:
    anchor=g['anchor'];face=anchor['face'];w=np.array(anchor['barycentric']);point=np.array(anchor['point'])
    assert 0<=face<outer_count and w.shape==(3,)and np.all(w>=0)and np.all(w<=1)and abs(w.sum()-1)<1e-14
    actual=w@mesh.triangles[face];error=float(np.linalg.norm(actual-point))
    assert error<64*np.finfo(float).eps*max(1.,float(np.abs(mesh.triangles[face]).max()))
    record=dict(term=g['term'],representation=g['representation'],anchorReplayErrorM=error,addedTissueSolids=0,anatomicalAcceptance=False)
    if g['representation']=='surface-region':
        ids=np.asarray(g['surface']['faces'],int)
        assert len(np.unique(ids))==len(ids)and np.all(ids>=0)and np.all(ids<outer_count)and face in ids
        tri=mesh.triangles[ids]
        record.update(facets=len(ids),areaM2=float(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()/2),boundsM=[tri.min((0,1)).tolist(),tri.max((0,1)).tolist()],boundaryProvisional=True)
    else:
        assert g['completeOpeningRepresented']is False
        record.update(completeOpeningRepresented=False,addedCompletedTa98Terms=0)
    checks.append(record)
atomic_json(folder/'selection-evidence.json',dict(passed=True,candidateSha256=digest_file(path),
    hostGeometrySha256=mesh.geometry_sha,bankSha256=mesh.source['bank_sha256'],allOriginalPartsAndChunksLiteral=True,
    checks=checks,surfaceAnnotations=18,pointLandmarks=1,addedSolids=0,installed=False,anatomicalAcceptance=False))
print('18 actual facet masks and one literal rim point replay; original anatomy unchanged.')
