"""Supervise small independent encoded-cell audits with durable results.

No deadline is extended or job retried. A failed batch stops this stage. The
final whole-bank check remains mandatory; batch completion never installs it.
"""
from pathlib import Path
import argparse,concurrent.futures,json,os,subprocess,sys,time
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser()
for name in ['source','stage','output','checkpoints','record']:ap.add_argument('--'+name,required=True)
ap.add_argument('--batch',type=int,default=64);ap.add_argument('--workers',type=int,default=4)
ap.add_argument('--batch-seconds',type=float,default=35);ap.add_argument('--seconds',type=float,default=210)
ap.add_argument('--face-budget',type=int,default=0,help='Limit the actual triangle workload per batch as well as its cell count')
args=ap.parse_args();started=time.monotonic();record=Path(args.record).resolve();record.mkdir(parents=True,exist_ok=False)
stage=Path(args.stage).resolve();meta=json.loads((stage/'pending.json').read_text())
assert 1<=args.workers<=4 and 1<=args.batch<=128 and 1<=args.batch_seconds<=45
source=Path(args.source).resolve();output=Path(args.output).resolve();checkpoints=Path(args.checkpoints).resolve()
count=len(meta['parts']);parts=list(meta['parts'].values());rows=[];env=dict(os.environ,PYTHONPATH=str(ROOT),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
batches=[];begin=0
while begin<count:
    size=0;faces=0
    while begin+size<count and size<args.batch:
        additional=parts[begin+size]['indexCount']//3
        if size and args.face_budget and faces+additional>args.face_budget:break
        size+=1;faces+=additional
    batches.append((begin,size));begin+=size
def job(batch):
    begin,size=batch
    if time.monotonic()-started+args.batch_seconds>args.seconds:raise TimeoutError('No new batch may outlive the stage deadline')
    log=record/f'batch-{begin}.log';runtime=record/f'batch-{begin}-runtime.json';t=time.monotonic()
    command=[sys.executable,'-B','-m','geometry_tools.runtime','--seconds',str(args.batch_seconds),'--report',str(runtime),'--',
        str(ROOT/'.local-models/geometry-tools-env/bin/python'),'-B','server/pack-spinal-cord-core.py',str(source),
        '--from-stage',str(stage),'--output',str(output),'--checkpoints',str(checkpoints),
        '--precheck-start',str(begin),'--precheck-count',str(size)]
    with log.open('w')as stream:code=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,check=False).returncode
    return dict(start=begin,count=size,returncode=code,elapsed_seconds=time.monotonic()-t,log=str(log),runtime=str(runtime))
iterator=iter(batches);failed=False
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers)as executor:
    pending={executor.submit(job,next(iterator)):None for _ in range(min(args.workers,len(batches)))}
    while pending:
        finished,_=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
        for future in finished:
            pending.pop(future);row=future.result();rows.append(row);failed=failed or row['returncode']!=0
            atomic_json(record/'batches.json',dict(bank_sha256=meta['bank_sha256'],source_sha256=digest_file(source/'cells.npz'),
                rows=rows,complete=False,installed=False,elapsed_seconds=time.monotonic()-started),overwrite=True)
            print(json.dumps(dict(completed_cells=sum(r['count']for r in rows if r['returncode']==0),total_cells=count,
                batch_start=row['start'],returncode=row['returncode'],elapsed_seconds=time.monotonic()-started)),flush=True)
            if not failed:
                begin=next(iterator,None)
                if begin is not None:pending[executor.submit(job,begin)]=None
assert not failed,'An encoded batch failed; inspect its saved log and do not install'
assert sum(r['count']for r in rows)==count
atomic_json(record/'batches.json',dict(bank_sha256=meta['bank_sha256'],source_sha256=digest_file(source/'cells.npz'),
    rows=rows,complete=True,installed=False,full_occupied_volume_check_pending=True,elapsed_seconds=time.monotonic()-started),overwrite=True)
