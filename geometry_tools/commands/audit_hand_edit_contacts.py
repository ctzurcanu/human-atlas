"""Compare edited source regions with every eligible fixed atlas neighbour.

Counts source triangle crossings, not tissue penetration or anatomical validity.
Crop boundaries are diagnostic only. Geometry below the declared coordinate
precision is unchanged for this comparison; the bank audit remains separate.
"""
import argparse,json,time
from collections import Counter
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs
from geometry_tools.course_io import sha,write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['before','candidate','output']:parser.add_argument('--'+key,required=True)
    parser.add_argument('--coordinate-precision-m',type=float,default=1e-7)
    args=parser.parse_args();start=time.monotonic();precision=args.coordinate_precision_m
    if not 0<=precision<=1e-6:parser.error('Explicit coordinate precision must be at most one micrometre')
    before_path=Path(args.before);candidate_path=Path(args.candidate);out=Path(args.output)
    if out.exists():raise ValueError('Use a new contact report')
    atlas=json.loads(before_path.read_text());data=json.loads(candidate_path.read_text());changed={r['partId']for r in data['rows']}
    joins={r['firstPart']:r for r in data.get('constructedJoins',[])};patches=[]
    for row in data['rows']:
        old=Mesh.load(row['sourcePath']);new=Mesh.load(row['path']);join=joins.get(row['partId'])
        # A current-bank comparison already includes a historical annular
        # reconstruction. Do not apply its older deletion map a second time.
        if np.array_equal(old.faces,new.faces):join=None
        mapping=np.r_[join['keptSourceFaceIds'],np.full(len(join['constructedFaceIds']),-1)]if join else np.arange(len(new.faces))
        has=np.flatnonzero(mapping>=0);motion=np.linalg.norm(new.triangles[has]-old.triangles[mapping[has]],axis=2).max(1)
        live=has[motion>precision];old_ids=list(mapping[live]);new_ids=list(live)
        if join:old_ids+=join['sourceRemovedFaceIds'];new_ids+=join['constructedFaceIds']
        if not old_ids:continue
        a,b=old.triangles[np.unique(old_ids)],new.triangles[np.unique(new_ids)];points=np.r_[a.reshape(-1,3),b.reshape(-1,3)]
        patches.append(dict(partId=row['partId'],before=a,after=b,bounds=np.array([points.min(0)-precision,points.max(0)+precision])))
    reader=AtlasReader();triangles=[];owners=[];sources=[]
    for part in atlas['parts']:
        if part['id']in changed or part['id']in atlas['defaultHidden']or part.get('suppressed')or not part.get('indexCount'):continue
        box=np.array(part.get('bounds',[]))
        if box.shape!=(2,3):continue
        applicable=[p['bounds']for p in patches if np.all(box[1]>=p['bounds'][0])and np.all(box[0]<=p['bounds'][1])]
        if not applicable:continue
        mesh=reader.read_part(atlas,part,frame='model:local-ta98:assembled-metres',catalogue_path=before_path);t=mesh.triangles;mask=np.zeros(len(t),bool)
        for box in applicable:mask|=np.all(t.max(1)>=box[0],axis=1)&np.all(t.min(1)<=box[1],axis=1)
        if mask.any():
            triangles.append(t[mask]);owners.extend([part['id']]*int(mask.sum()))
            sources.append(dict(partId=part['id'],sourceGeometrySha256=mesh.geometry_sha,sourceFaceIds=np.flatnonzero(mask).tolist()))
        if time.monotonic()-start>22:raise TimeoutError('22-second native recipient gate')
    rows=[]
    if triangles:
        t=np.concatenate(triangles);fixed=BVH(Mesh(t.reshape(-1,3),np.arange(t.size//3).reshape(-1,3),'model:local-ta98:assembled-metres'))
        for patch in patches:
            counts=[]
            for key in ['before','after']:
                t=patch[key];own=BVH(Mesh(t.reshape(-1,3),np.arange(t.size//3).reshape(-1,3),fixed.mesh.frame));c=Counter();tests=0
                for batch in own.candidate_batches(fixed,0.,batch_size=512):
                    tests+=len(batch)
                    if tests>500000:raise ValueError('Finite recipient-pair budget reached')
                    q=triangle_pairs(own.triangles[batch[:,0]],fixed.triangles[batch[:,1]])
                    for k in np.flatnonzero(np.isin(q['relation'],['crossing','coplanar-overlap'])):c[owners[int(batch[k,1])]]+=1
                    if time.monotonic()-start>22:raise TimeoutError('22-second native recipient gate')
                counts.append(c)
            regress={pid:[counts[0][pid],n]for pid,n in counts[1].items()if n>counts[0][pid]}
            rows.append(dict(partId=patch['partId'],before=dict(counts[0]),after=dict(counts[1]),regressions=regress))
            print(patch['partId'].split('-')[-1],sum(counts[0].values()),sum(counts[1].values()),regress,flush=True)
    write_json(out,dict(schema='hand-source-edited-region-recipient-gate/v1',catalogueSha256=sha(before_path),candidateSha256=sha(candidate_path),
        coordinatePrecision_m=precision,contextSources=sources,rows=rows,complete=True,noContactCountRegressions=not any(r['regressions']for r in rows),
        remainingWholeHandFitNotAccepted=True,sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-start))


if __name__=='__main__':main()
