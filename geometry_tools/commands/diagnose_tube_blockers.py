"""Identify named recipient owners against intact, source-bound tube cuts.

python -m geometry_tools.commands.diagnose_tube_blockers SPEC --output NEW_REPORT
Section overlap is diagnostic evidence; it is not a whole-surface test or proof
that a vessel/tissue interface is unintended. No anatomy is changed.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import shapely
from geometry_tools.course_io import load_course_spec,sha,write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spec');parser.add_argument('--output',required=True)
    args=parser.parse_args();started=time.monotonic();out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new report path')
    spec,g,_,_,course,*_=load_course_spec(args.spec)
    byindex={}
    for binding in course['sectionBindings']:
        path=Path(binding['path'])
        if sha(path)!=binding['sha256']:raise ValueError('Stale named context')
        data=json.loads(path.read_text())
        for row in data['rows']:
            # A source report may retain failed rows outside its bound range.
            if row['index']in binding['indices']:byindex[row['index']]=row
    rows=[]
    for station in course['rows']:
        shapes={}
        for outline in byindex[station['index']]['sourceContextOutlines']:
            if outline['kind']!='solid':continue
            p=shapely.polygons(outline['points'])
            if not shapely.is_valid(p):raise ValueError('Invalid named source outline')
            name=outline['name']
            shapes[name]=shapely.symmetric_difference(shapes[name],p)if name in shapes else p
        for profile in station['profiles']:
            edge=np.asarray(profile['sourceRawEdges']);weights=np.asarray(profile['sourceEdgeWeights'])
            points=np.sum(g['vertices'][edge]*weights[:,:,None],axis=1)
            local=np.asarray(profile['sourcePolygonXZ_m'])+profile['referenceCentreXZ_m']
            if np.max(abs(local-points[:,[0,2]]))>1e-12:raise ValueError('Source section coordinate frame mismatch')
            p=shapely.polygons(points[:,[0,2]])
            if not shapely.is_valid(p):raise ValueError('Invalid intact source contour')
            for name,shape in shapes.items():
                if shapely.intersects(p,shape):rows.append(dict(index=station['index'],y_m=station['y_m'],
                    partId=profile['partId'],contextOwner=name,
                    overlapArea_m2=float(shapely.area(shapely.intersection(p,shape)))))
        if time.monotonic()-started>15:raise TimeoutError('15s native named-blocker diagnosis deadline')
    names=sorted({r['contextOwner']for r in rows})
    summary=[]
    for name in names:
        hits=[r for r in rows if r['contextOwner']==name]
        summary.append(dict(owner=name,profileSectionHits=len(hits),
            maximumOverlapArea_m2=max(r['overlapArea_m2']for r in hits),
            stationIndices=sorted({r['index']for r in hits})))
    report=dict(schema='source-bound-named-blocker-diagnosis/v1',specSha256=sha(args.spec),
        sourceCodeSha256=sha(__file__),stationsChecked=len(course['rows']),owners=summary,
        sectionEvidence=rows,geometryChanged=False,wholeSurfaceChecksPending=True,
        anatomicalAttachmentOrCollisionInferred=False,installed=False,elapsedSeconds=time.monotonic()-started)
    out.parent.mkdir(parents=True,exist_ok=True);write_json(out,report)
    print(json.dumps({k:v for k,v in report.items()if k!='sectionEvidence'}))


if __name__=='__main__':main()
