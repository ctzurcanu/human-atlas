"""Run one saved segment through separately supervised construction stages."""
from pathlib import Path
import argparse,json,subprocess,sys,time,os
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--segment',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
started=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';python=str(ROOT/'.local-models/geometry-tools-env/bin/python')
jobs=[('parcels',35,[python,'-B','server/construct-spinal-tract-parcels.py','--fields',str(base/'female-tract-fields-v1'),
    '--segment',args.segment,'--seconds','25','--output',str(out/'parcels')]),
    ('intersect',20,[python,'-B','server/intersect-spinal-tract-parcels.py','intersect','--parcels',str(out/'parcels'),'--output',str(out/'boolean')]),
    ('charts',45,[python,'-B','server/intersect-spinal-tract-parcels.py','assemble','--parcels',str(out/'parcels'),
        '--source',str(out/'boolean'),'--output',str(out/'charts')]),
    ('facets',35,[python,'-B','server/condition-spinal-tract-facets.py',str(out/'charts'),'--output',str(out/'facets')])]
records=[]
for name,seconds,command in jobs:
    t=time.monotonic();runtime=out/(name+'-runtime.json');log=out/(name+'.log')
    with log.open('w')as stream:
        code=subprocess.run([sys.executable,'-B','-m','geometry_tools.runtime','--seconds',str(seconds),'--report',str(runtime),'--',*command],
            cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,check=False).returncode
    records.append(dict(stage=name,returncode=code,elapsed_seconds=time.monotonic()-t,runtime_report=str(runtime),log=str(log)))
    atomic_json(out/'stages.json',dict(segment=args.segment,stages=records,complete=False,installed=False),overwrite=True)
    if code:
        print(json.dumps(dict(segment=args.segment,failed=name,returncode=code,log_tail=log.read_text()[-2500:])),flush=True)
        raise SystemExit(code)
    print(json.dumps(dict(segment=args.segment,stage=name,elapsed_seconds=time.monotonic()-t)),flush=True)
atomic_json(out/'stages.json',dict(segment=args.segment,stages=records,complete=True,
    source_sha256=digest_file(out/'facets/cells.npz'),full_emitted_audit_pending=True,installed=False,
    elapsed_seconds=time.monotonic()-started),overwrite=True)
print(json.dumps(dict(segment=args.segment,constructed=True,seconds=time.monotonic()-started)),flush=True)
