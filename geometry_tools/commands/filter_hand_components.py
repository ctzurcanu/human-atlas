"""Retain only complete source-component trials with no detected regressions.

Uses actual original triangles for context crossing counts and native self-pair
proposals confirmed by Float64 predicates. This is a diagnostic filter, not an
anatomical acceptance or atlas installer.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs
from geometry_tools.open3d_backend import diagnostics


def encoded_attribute(name,source,dtype):
    """Normalize normal representation when mixing trial and packed sources."""
    if name=='normals'and np.issubdtype(dtype,np.integer)and np.issubdtype(source.dtype,np.floating):
        length=np.linalg.norm(source,axis=1)
        if np.any(length==0):raise ValueError('Zero trial normal cannot be encoded')
        unit=source/length[:,None]
        return np.rint(np.clip(unit,-1,1)*np.iinfo(dtype).max).astype(dtype)
    return source.astype(dtype,copy=False)


def crossings(mesh,tree,start):
    first=BVH(mesh);count=0;tests=0
    for batch in first.candidate_batches(tree,0.,batch_size=512):
        tests+=len(batch)
        if tests>100000:raise ValueError('Finite component/context pair budget reached')
        r=triangle_pairs(first.triangles[batch[:,0]],tree.triangles[batch[:,1]])
        count+=int(np.isin(r['relation'],['crossing','coplanar-overlap']).sum())
        if time.monotonic()-start>22:raise TimeoutError('22-second native component filter')
    return count


def self_pairs(mesh):
    d=diagnostics(mesh)
    if d['self_intersection_pairs']>20:return None,d['self_intersection_pairs']
    return set(map(tuple,d['source_face_pair_examples'])),d['self_intersection_pairs']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['base','binding','bones','skin-reference','output']:parser.add_argument('--'+name,required=True)
    parser.add_argument('--trial',action='append',required=True)
    args=parser.parse_args();start=time.monotonic();paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in ['base','binding','bones','skin-reference']}
    base,binding,bones,sr=[json.loads(p.read_text())for p in paths.values()];side=base['side'];suffix='.l'if side=='left'else'.r'
    trials=[(Path(p).resolve(),json.loads(Path(p).read_text()))for p in args.trial]
    fixed=[]
    for r in bones['rows']:
        if r['sourceObject'].endswith(suffix):
            if sha(r['path'])!=r['sha256']:raise ValueError('Fixed skeleton changed')
            fixed.append(Mesh.load(r['path']).triangles)
    t=np.concatenate(fixed);bone=BVH(Mesh(t.reshape(-1,3),np.arange(t.size//3).reshape(-1,3),'model:local-ta98:assembled-metres'))
    if sha(sr['path'])!=sr['sha256']or not sr['defaultHiddenExcluded']:raise ValueError('Visible skin reference changed')
    skin=BVH(Mesh.load(sr['path']));out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new filtered stage')
    out.mkdir();rows=[]
    for row in base['rows']:
        if sha(row['path'])!=row['sha256']:raise ValueError('Baseline changed')
        original=Mesh.load(row['path']);v=original.vertices.copy();attrs={k:a.copy()for k,a in original.attributes.items()}
        sources=[]
        for report_path,trial in trials:
            q=next(r for r in trial['rows']if r['partId']==row['partId'])
            if sha(q['path'])!=q['sha256']or q['sourceSha256']!=row['sourceSha256']:raise ValueError('Trial source changed')
            mesh=Mesh.load(q['path'])
            if not np.array_equal(mesh.faces,original.faces)or mesh.frame!=original.frame:raise ValueError('Trial changed source face ownership/frame')
            sources.append((report_path,q,mesh))
        owners=[r for r in binding['rows']if r['side']==side and r['partId']==row['partId']]
        face_owner=np.full(len(original.faces),-1,int);changes={};records=[]
        for owner in owners:
            k=owner['physicalComponent'];fids=np.array(owner['sourceFaceIds'],int);ids=np.array(owner['rawVertexIds'],int);face_owner[fids]=k
            source=Mesh(original.vertices,original.faces[fids],original.frame)
            baseline_pairs,bp=self_pairs(source)
            if baseline_pairs is None:raise ValueError('Baseline self diagnostic exceeds finite evidence size')
            baseline_counts=[crossings(source,tree,start)for tree in [bone,skin]];best=baseline_counts;selected=None;checks=[]
            for report_path,q,m in sources:
                proof=q.get('clearanceProposal',{})
                component_proof=next((p for p in proof.get('components',[])if p['sourceComponent']==k),None)
                complete=component_proof['complete']if component_proof else proof.get('complete',False)
                if not complete:checks.append(dict(trial=str(report_path),rejected='Incomplete source component'));continue
                candidate=Mesh(m.vertices,m.faces[fids],m.frame);pairs,count=self_pairs(candidate)
                if pairs is None or not pairs<=baseline_pairs:
                    checks.append(dict(trial=str(report_path),rejected='Additional intrinsic pairs',pairs=count));continue
                counts=[crossings(candidate,tree,start)for tree in [bone,skin]]
                if all(a<=b for a,b in zip(counts,best))and sum(counts)<sum(best):best=counts;selected=(report_path,q,m)
                checks.append(dict(trial=str(report_path),boneCrossingPairs=counts[0],skinCrossingPairs=counts[1],intrinsicPairs=count))
            if selected:
                report_path,q,m=selected;v[ids]=m.vertices[ids]
                for name in attrs:attrs[name][ids]=encoded_attribute(name,m.attributes[name],attrs[name].dtype)[ids]
                changes[k]=ids
            records.append(dict(sourceComponent=k,baselineBoneSkinPairs=baseline_counts,selectedBoneSkinPairs=best,
                selectedTrial=str(selected[0])if selected else None,checks=checks,changed=bool(selected)))
        source_pairs,source_count=self_pairs(original)
        if source_pairs is None:raise ValueError('Baseline has too many source self pairs')
        # Independent components can collide after individually valid changes.
        # Revert both changed owners of every new pair; never conceal a defect.
        for _ in range(len(owners)+1):
            composite=Mesh(v,original.faces,original.frame,attrs);pairs,count=self_pairs(composite)
            extra=(pairs-source_pairs)if pairs is not None else set()
            if pairs is not None and not extra:break
            reject=set(changes)if pairs is None else {int(face_owner[f])for pair in extra for f in pair}&set(changes)
            if not reject:raise ValueError('Unresolved whole-part self regression')
            for k in reject:
                ids=changes.pop(k);v[ids]=original.vertices[ids]
                for name in attrs:attrs[name][ids]=original.attributes[name][ids]
                record=next(r for r in records if r['sourceComponent']==k);record.update(changed=False,revertedForWholePartCrossing=True)
        else:raise ValueError('Whole-part self filter did not terminate')
        path=out/f'part-{len(rows):02}.npz';composite.save(path)
        rows.append(dict(row,path=str(path),sha256=sha(path),geometrySha256=composite.geometry_sha,
            parentPosePath=row['path'],parentPoseSha256=row['sha256'],sourceComponentShapesPreserved=False,
            componentChecks=records,newDetectedIntrinsicPairs=0,retainedComponentChanges=len(changes)))
    report=dict(schema='nonregressing-hand-component-diagnostic/v1',side=side,rows=rows,installed=False,installationAllowed=False,
        inputs={k:dict(path=str(p),sha256=sha(p))for k,p in paths.items()},sourceCodeSha256=sha(__file__),
        skeletonChanged=False,sourceFacesRetained=True,actualTriangleContextChecks=True,elapsedSeconds=time.monotonic()-start)
    write_json(out/'report.json',report);print(json.dumps(dict(parts=len(rows),retainedComponents=sum(r['retainedComponentChanges']for r in rows),elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
