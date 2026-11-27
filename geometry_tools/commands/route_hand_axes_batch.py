"""Independently saved finger-route stages sharing one bounded native context."""
import argparse,json,sys,time
from pathlib import Path
from geometry_tools.course_io import sha,write_json
from geometry_tools.commands.route_hand_axis import main as route_one


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['axis-inventory','course-matches','bones','skin-reference','recipient-context','output']:p.add_argument('--'+k,required=True)
    p.add_argument('--side',choices=['left','right'],required=True)
    p.add_argument('--transverse-neighbour-cells',type=int,default=2)
    args=p.parse_args();start=time.monotonic();paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in
        ['axis-inventory','course-matches','bones','skin-reference','recipient-context']}
    source=json.loads(paths['axis-inventory'].read_text());out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new independent route-batch stage')
    out.mkdir();records=[];selected=[r for r in source['rows']if '-'+args.side+'_'in r['partId']]
    for row in selected:
        remaining=21-(time.monotonic()-start)
        if not row['complete']:
            records.append(dict(partId=row['partId'],component=row['component'],complete=False,status='SourceAxisUnresolved'));continue
        if remaining<4:
            records.append(dict(partId=row['partId'],component=row['component'],complete=False,status='DeferredBeforeStartingAtBatchDeadline'));continue
        if sha(row['path'])!=row['sha256']:raise ValueError('Source axis record changed')
        dest=out/('component-'+str(row['component']).zfill(2));before=sys.argv
        try:
            sys.argv=['route-hand-axis','--source-axis-report',row['path'],'--course-match-report',str(paths['course-matches']),
                '--bones',str(paths['bones']),'--skin-reference',str(paths['skin-reference']),
                '--recipient-context',str(paths['recipient-context']),'--free-endpoints',
                '--transverse-neighbour-cells',str(args.transverse_neighbour_cells),
                '--native-route-seconds',str(min(16.,remaining-1.)), '--output',str(dest)]
            route_one();saved=dest/'report.json';result=json.loads(saved.read_text())
            record=dict(partId=row['partId'],component=row['component'],complete=result['proof']['complete'],
                status=result['proof']['status'],path=str(saved),sha256=sha(saved))
        except (ValueError,TimeoutError)as exc:
            record=dict(partId=row['partId'],component=row['component'],complete=False,status='ProposalFailed',error=str(exc))
        finally:sys.argv=before
        records.append(record)
        write_json(out/'progress.json',dict(rows=records,installed=False,elapsedSeconds=time.monotonic()-start))
    report=dict(complete=all(r['complete']for r in records),rows=records,sourceGeometryChanged=False,
        skeletonChanged=False,installed=False,installationAllowed=False,
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
        sourceCodeSha256=sha(__file__),nativeBatchDeadlineSeconds=21.,elapsedSeconds=time.monotonic()-start)
    write_json(out/'report.json',report)
    print(json.dumps(dict(completed=sum(r['complete']for r in records),total=len(records),
                         statuses=[dict(component=r['component'],status=r['status'])for r in records],elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
