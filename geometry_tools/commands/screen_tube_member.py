"""Screen one stored vascular member against its full bound recipient context.

python -m geometry_tools.commands.screen_tube_member PREVIEW --part ID --output REPORT
Optional --recipient-replacement DIRECTORY supplies hash-bound tissue candidates.
This reports actual surface crossings, not anatomical/containment acceptance.
"""
import argparse,json,time
from pathlib import Path
from collections import Counter
import numpy as np
from geometry_tools.course_io import load_arrays,load_course_spec,sha,write_json
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('preview');parser.add_argument('--part',required=True)
    parser.add_argument('--recipient-replacement');parser.add_argument('--output',required=True)
    args=parser.parse_args();started=time.monotonic();out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new report path')
    folder=Path(args.preview).resolve();preview=json.loads((folder/'report.json').read_text())
    spec,g,source,plan,course,*_=load_course_spec(preview['specPath'])
    row=next(r for r in preview['rows']if r['partId']==args.part)
    path=folder/row['geometryFile']
    if sha(path)!=row['sha256']:raise ValueError('Stored member changed')
    data=load_arrays(path);member=Mesh(data['vertices'],data['faces'],'atlas-world-metres')
    bindings={}
    for b in course['sectionBindings']:
        if sha(b['path'])!=b['sha256']:raise ValueError('Recipient source binding changed')
        for r in b['fixedTargetBindings']:
            if r['partId']in bindings and bindings[r['partId']]!=r['sourceGeometrySha256']:raise ValueError('Conflicting recipient identities')
            bindings[r['partId']]=r['sourceGeometrySha256']
    replacement={};arrays=None;replacement_sha=None
    if args.recipient_replacement:
        rp=Path(args.recipient_replacement);r=json.loads((rp/'report.json').read_text())
        if sha(rp/'candidate.npz')!=r['candidateArraysSha256']:raise ValueError('Recipient candidate changed')
        arrays=load_arrays(rp/'candidate.npz');replacement={p['partId']:p for p in r['records']};replacement_sha=sha(rp/'report.json')
    atlas=json.loads(Path(spec['catalogue']).read_text());by={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    targets=[];owners=[];context=[];lo,hi=member.bounds
    for pid,expected in bindings.items():
        if pid==args.part:continue
        part=by[pid]
        if part.get('suppressed'):raise ValueError('Context unexpectedly contains a suppressed donor')
        m=reader.read_part(atlas,part,frame=member.frame,catalogue_path=spec['catalogue'])
        if m.geometry_sha!=expected:raise ValueError('Recipient geometry changed: '+pid)
        replaced=False
        if pid in replacement:
            r=replacement[pid]
            if r['sourceGeometrySha256']!=m.geometry_sha:raise ValueError('Replacement source mismatch')
            k=r['arrayIndex'];m=Mesh(arrays[f'positions_{k}'],arrays[f'faces_{k}'],m.frame);replaced=True
        tri=m.triangles;ids=np.flatnonzero((tri.max(1)>=lo).all(1)&(tri.min(1)<=hi).all(1))
        if len(ids):targets.append(tri[ids]);owners.extend((pid,int(i))for i in ids)
        context.append(dict(partId=pid,sourceGeometrySha256=expected,recipientReplacement=replaced,broadPhaseFaces=len(ids)))
        if time.monotonic()-started>20:raise TimeoutError('20s native complete-context screen deadline')
    target=np.concatenate(targets);mesh=Mesh(target.reshape(-1,3),np.arange(target.size//3).reshape(-1,3),member.frame)
    first,second=BVH(member),BVH(mesh);tested=0;counts=Counter();witnesses=[]
    for batch in first.candidate_batches(second,0.,batch_size=1024):
        tested+=len(batch)
        if tested>500000:raise ValueError('Finite500000 triangle-test budget reached')
        result=triangle_pairs(first.triangles[batch[:,0]],second.triangles[batch[:,1]])
        for i in np.flatnonzero(np.isin(result['relation'],['crossing','coplanar-overlap'])):
            source_face,target_face=map(int,batch[i]);pid,face=owners[target_face];counts[pid]+=1
            witnesses.append(dict(partId=pid,memberFace=source_face,targetFace=face,
                kind=str(result['relation'][i]),point_m=result['first'][i].tolist()))
        if time.monotonic()-started>20:raise TimeoutError('20s native complete-context screen deadline')
    report=dict(schema='bound-vessel-member-surface-screen/v1',partId=args.part,
        sourceArraySha256=sha(path),sourceCatalogueSha256=sha(spec['catalogue']),
        recipientReplacementReportSha256=replacement_sha,actualContext=context,
        pairTests=tested,crossingFacePairsByOwner=dict(counts),witnesses=witnesses,
        wholeBoundContextScreenComplete=True,containmentAndAnatomicalInterfacesPending=True,
        installed=False,anatomicalAcceptance=False,sourceCodeSha256=sha(__file__),
        elapsedSeconds=time.monotonic()-started)
    out.parent.mkdir(parents=True,exist_ok=True);write_json(out,report)
    print(json.dumps({k:v for k,v in report.items()if k not in ['actualContext','witnesses']}))


if __name__=='__main__':main()
