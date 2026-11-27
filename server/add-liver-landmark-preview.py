"""Demonstrate a non-solid, source-bound point on the actual hilar rim."""
from pathlib import Path
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

folder=ROOT/'.local-models/geometry-workbench/liver-completion/liver-surface-controls-v143'
atlas=json.loads((folder/'candidate.json').read_text())
host=next(p for p in atlas['parts']if p['id']=='TA98PART:ZA:Liver:capsule')
mesh=AtlasReader().read_part(atlas,host,frame='model:local-ta98:assembled-metres')
count=host['tissueDrawGroups'][0]['count']//3
v,back=np.unique(mesh.vertices,axis=0,return_inverse=True);f=back[mesh.faces[:count]]
loops=boundary_loops(f)
assert len(loops)==3,'The expected literal source hila changed'
# An actual rim vertex anchors the ordinary label. No centroid in a hole and
# no fabricated marker solid substitutes for the full anatomical opening.
rim=np.concatenate(loops);vertex=int(rim[np.argmax(v[rim,2])]);face=int(np.flatnonzero(np.any(f==vertex,axis=1))[0])
corner=int(np.flatnonzero(f[face]==vertex)[0]);weights=np.eye(3)[corner]
reference=dict(term='A05.8.01.016',conceptId='TA98LIVER:male:landmark:A05.8.01.016',partIds=[host['id']],
    representation='point-landmark',completeOpeningRepresented=False,
    anchor=dict(partId=host['id'],face=face,barycentric=weights.tolist(),point=v[vertex].tolist()))
atlas['liverLandmarkIntegration']=dict(schema='liver-hilar-landmark/v1',hostPartId=host['id'],
    bankSha256=atlas['liverIntegration']['bankSha256'],termGroups=[reference])
atlas['concepts'].append(dict(id=reference['conceptId'],name='Porta hepatis',elements=[host['id']],
    ta98Term=reference['term'],ta98Kind='exact',partOf='TA98LIVER:male:surface:A05.8.01.012',
    description='Landmark on the actual hilar rim. The complete extent of the porta hepatis remains to be delineated.'))
atomic_json(folder/'with-landmark.json',atlas)
atomic_json(folder/'landmark-evidence.json',dict(reference=reference,hostGeometrySha256=mesh.geometry_sha,
    sourceRimLoops=[len(loop)for loop in loops],sourceBankSha256=atlas['liverIntegration']['bankSha256'],
    pointOnLiteralFacet=True,addedSolids=0,addedCompletedTa98Terms=0,installed=False))
print('Actual hilar rim point prepared; no additional anatomical solid or completed term.')
