"""First bounded artery correction using the proven literal global oracle.

The v37 search is retained as a deterministic source-profile corridor proposal,
but its failed nearest-sheet reconstruction is not rerun. Complete native wall
and finite-envelope checks still decide admissibility after construction.
"""
from pathlib import Path
import sys,json,time,types
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,digest_file,atomic_json
from geometry_tools.model import Mesh

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
ORACLE=BASE/'biliary-actual-envelope-v39'

def actual_duct(actual,out):
    proof=json.loads((ORACLE/'proof.json').read_text())
    assert proof['passed']and actual.geometry_sha==proof['actualWallGeometrySha256']and actual.source['bank_sha256']==proof['actualBankSha256']
    mesh=Mesh.load(ORACLE/'actual-biliary-outer-computational-envelope.npz')
    assert mesh.geometry_sha==proof['envelopes'][0]['geometrySha256']
    mesh.save(out/'actual-biliary-outer-computational-envelope.npz')
    atomic_json(out/'actual-biliary-envelope-proof.json',proof)
    return mesh,dict(proofPath=str(ORACLE/'proof.json'),proofSha256=digest_file(ORACLE/'proof.json'),**proof)

def actual_vein(actual,partition,out):
    stage=BASE/'cystic-complete-duct-occupancy-v40';proof=json.loads((stage/'proof.json').read_text());row=next(r for r in proof['results']if r['id']=='TA98CYSTICBLOOD:male:vein')
    assert row['passedNoUnapprovedDuctOccupation']and row['recovery']['actualGeometrySha256']==actual.geometry_sha and row['recovery']['actualBankSha256']==actual.source['bank_sha256']and row['recovery']['partitionSha256']==digest_file(partition)
    mesh=Mesh.load(stage/'fixed-installed-vein-actual-outer-computational-envelope.npz')
    assert mesh.geometry_sha==row['recovery']['envelopes'][0]['geometrySha256']
    return mesh,dict(proofPath=str(stage/'proof.json'),proofSha256=digest_file(stage/'proof.json'),**row['recovery'])

def main():
    source=(ROOT/'server/construct-cystic-artery-exterior-corridor-v37.py').read_text()
    source=source.replace("OUT=WORK/'artery-exterior-corridor-v37'","OUT=WORK/'artery-valid-envelope-corridor-v41'")
    source=source.replace("envelope,envelope_proof=em.reconstruct(blockers[0][1],OUT)","envelope,envelope_proof=actual_duct(blockers[0][1],OUT)")
    source=source.replace("fixed_vein_envelope,fixed_vein_envelope_proof=em.vessel_envelope(blockers[1][1],WORK/'vein-recipient-partition-v28/vein/partition.npz',OUT)","fixed_vein_envelope,fixed_vein_envelope_proof=actual_vein(blockers[1][1],WORK/'vein-recipient-partition-v28/vein/partition.npz',OUT)")
    source=source.replace("currentStage=3","currentStage=1").replace("stage=3,maximumStages=3","stage=1,maximumStages=3")
    source=source.replace("/v37","/v41")
    module=types.ModuleType('valid_envelope_artery_stage1');module.__file__=str(ROOT/'server/construct-cystic-artery-exterior-corridor-v37.py');module.actual_duct=actual_duct;module.actual_vein=actual_vein
    exec(compile(source,module.__file__,'exec'),module.__dict__)
    module.main()

if __name__=='__main__':main()
