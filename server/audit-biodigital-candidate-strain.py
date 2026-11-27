"""Compare complete source-matched vessel faces, not just sampled landmarks."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
B=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
ap=argparse.ArgumentParser();ap.add_argument('--baseline-stage',required=True);ap.add_argument('--candidate-stage',required=True);ap.add_argument('--output-name',default='complete-source-strain-v1');a=ap.parse_args();started=time.monotonic()
base=B/a.baseline_stage/'candidate-v1';candidate=B/a.candidate_stage/'candidate-v1';out=B/a.candidate_stage/a.output_name;out.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=json.loads((candidate/'instance-ownership.json').read_text())['instances'];oldrows=json.loads((base/'instance-ownership.json').read_text())['instances'];assert rows==oldrows
assert json.loads((candidate/'candidate-report.json').read_text())['status']=='candidate-numerical-gates-passed-anatomy-unverified'
records=[];inputs={};unchanged=[]
for meshname in ['arterial','venous','arterial-remnant-companion']:
    p0,p1=base/(meshname+'.npz'),candidate/(meshname+'.npz');z0,z1=np.load(p0),np.load(p1);v0,v1=z0['vertices'],z1['vertices'];f0,f1=z0['faces'],z1['faces'];assert np.array_equal(f0,f1)
    x0,x1=base/(meshname+'-ancestry.npz'),candidate/(meshname+'-ancestry.npz');q0,q1=np.load(x0),np.load(x1)
    assert set(q0.files)==set(q1.files)and all(np.array_equal(q0[k],q1[k])for k in q0.files)
    inputs.update({str(p):sha(p)for p in [p0,p1,x0,x1]})
    for row in rows:
        if row['candidateMesh']!=meshname+'.npz':continue
        s,n=row['candidateVertexStart'],row['candidateVertexCount'];movement=np.linalg.norm(v1[s:s+n]-v0[s:s+n],axis=1)
        if not np.any(movement>1e-9):unchanged.append(row['sourceInstanceId']);continue
        s,nf=row['candidateFaceStart'],row['candidateFaceCount'];tri0=v0[f0[s:s+nf]];tri1=v1[f1[s:s+nf]]
        ar0=np.linalg.norm(np.cross(tri0[:,1]-tri0[:,0],tri0[:,2]-tri0[:,0]),axis=1);ar1=np.linalg.norm(np.cross(tri1[:,1]-tri1[:,0],tri1[:,2]-tri1[:,0]),axis=1);assert np.all(ar0>0)
        ratios=ar1/ar0
        records.append(dict(sourceInstanceId=row['sourceInstanceId'],name=row['name'],vertices_moved=int((movement>1e-9).sum()),faces=nf,maximumMovement_m=float(movement.max()),
            surfaceTriangleAreaRatio=dict(minimum=float(ratios.min()),maximum=float(ratios.max()),quantiles=np.quantile(ratios,[.01,.05,.5,.95,.99]).tolist(),belowHalf=int((ratios<.5).sum()),aboveDouble=int((ratios>2).sum()))))
    del v0,v1,f0,f1,z0,z1,q0,q1
records.sort(key=lambda r:r['surfaceTriangleAreaRatio']['minimum'])
report=dict(schema='complete-source-matched-vessel-deformation-comparison/v1',baseline=str(base),candidate=str(candidate),inputs=inputs,sourceParts=len(rows),
    changedParts=len(records),unchangedPartIds=unchanged,parts=records,allSourceFacesCompared=True,exactAncestryAndFaceOrder=True,
    surfaceAreaIsNotLumenSectionArea=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-started)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'output':str(out),'changedParts':len(records),'mostCompressed':records[:5],'elapsedSeconds':report['elapsedSeconds']}))
