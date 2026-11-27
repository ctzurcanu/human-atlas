"""Complete actual-triangle bone crossing comparison for a saved moved-face set."""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs


def main():
    started=time.perf_counter();parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--candidate-stage',required=True);parser.add_argument('--output-name',default='all-moved-face-bone-audit-v1');parser.add_argument('--max-faces',type=int,default=32768);args=parser.parse_args()
    assert 1<=args.max_faces<=500000
    assert Path(args.candidate_stage).name==args.candidate_stage and args.candidate_stage not in ('.','..')
    assert Path(args.output_name).name==args.output_name and args.output_name not in ('.','..')
    w=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';stage=w/args.candidate_stage;out=stage/args.output_name;out.mkdir(exist_ok=False)
    skin_evidence=stage/'moved-triangle-skin-regression-v1/actual-triangle-evidence.npz'
    with np.load(skin_evidence)as z:
        tb,ta=z['before_moved_triangles'],z['after_moved_triangles'];source_owners=z['source_instance_indices'];source_face_ids=z['source_face_ids']
    assert len(tb)<=args.max_faces, f'{len(tb)} moved faces exceed the explicit {args.max_faces} face budget'
    context_path=w/'regional-interface-v1/fixed-context/full.npz';meta_path=context_path.parent/'context.json';meta=json.loads(meta_path.read_text())
    with np.load(context_path)as z:
        v,f,owners=z['vertices'],z['faces'],z['face_owner_part_index']
    all_tri=v[f];lo=np.minimum(tb.min(axis=(0,1)),ta.min(axis=(0,1)));hi=np.maximum(tb.max(axis=(0,1)),ta.max(axis=(0,1)))
    eligible_parts=[r['partIndex']for r in meta['partOwnership']if r['sourceRecord']['system']=='skeletal']
    selected=(all_tri.max(1)>=lo).all(1)&(all_tri.min(1)<=hi).all(1)&np.isin(owners,eligible_parts)
    fixed_full_ids=np.flatnonzero(selected);bone_tri=all_tri[selected]
    bb=None
    if len(bone_tri):
        bone=Mesh(bone_tri.reshape(-1,3),np.arange(len(bone_tri)*3).reshape(-1,3),'active-male-metres');bb=BVH(bone)
    results=[]
    for label,tri in [('before',tb),('after',ta)]:
        mesh=Mesh(tri.reshape(-1,3),np.arange(len(tri)*3).reshape(-1,3),'active-male-metres');bvh=BVH(mesh);pairs=[];tested=0
        for batch in (bvh.candidate_batches(bb,0.,batch_size=2048) if bb is not None else []):
            if time.perf_counter()-started>35:raise TimeoutError('Internal35second complete-audit deadline')
            result=triangle_pairs(tri[batch[:,0]],bone_tri[batch[:,1]]);tested+=len(batch)
            pairs.extend(batch[result['relation']=='crossing'].tolist())
        results.append(dict(label=label,testedActualPairs=tested,crossingPairs=pairs))
    bs,as_=[set(map(tuple,r['crossingPairs']))for r in results]
    report=dict(schema='complete-actual-moved-face-bone-crossing/v1',candidateStage=args.candidate_stage,
        sourceMovedFaces=len(tb),fixedContextFaces=len(f),broadPhaseEligibleBoneFaces=len(bone_tri),
        fixedBonePartIds=[meta['partOwnership'][int(i)]['id']for i in np.unique(owners[selected])],
        sourceSkinEvidenceSha256=hashlib.sha256(skin_evidence.read_bytes()).hexdigest(),fixedContextSha256=hashlib.sha256(context_path.read_bytes()).hexdigest(),
        beforePairs=len(bs),afterPairs=len(as_),newPairs=len(as_-bs),resolvedPairs=len(bs-as_),
        results=results,geometryModified=False,anatomicalAcceptance=False,elapsedSeconds=time.perf_counter()-started)
    np.savez_compressed(out/'actual-bone-face-evidence.npz',source_instance_indices=source_owners,source_face_ids=source_face_ids,before_triangles=tb,after_triangles=ta,fixed_full_face_ids=fixed_full_ids,actual_fixed_triangles=bone_tri)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['sourceMovedFaces','broadPhaseEligibleBoneFaces','beforePairs','afterPairs','newPairs','elapsedSeconds']}))


if __name__=='__main__':main()
