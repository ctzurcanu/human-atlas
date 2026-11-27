"""Content-bound local jobs, durable state, reuse and reviewed artifacts."""
from pathlib import Path
import os,json,time,platform
from psycopg.types.json import Jsonb
from .common import CACHE,ROOT,VERSION,identity,digest_file,atomic_json,json_bytes


def implementation_signature():
    paths=[p for p in Path(__file__).parent.rglob('*.py')if 'tests'not in p.parts]
    paths.extend(p for p in Path(ROOT/'geometry').glob('*.py')if not p.name.startswith('test_'))
    from importlib.metadata import version,PackageNotFoundError
    versions={}
    for package in ['numpy','scipy','manifold3d','meshpy','point-cloud-utils','open3d','clarabel']:
        try:versions[package]=version(package)
        except PackageNotFoundError:versions[package]=None
    return identity(platform.python_version(),versions,[(str(p.relative_to(ROOT)),digest_file(p))for p in sorted(paths)])


def file_bindings(value):
    found={}
    def walk(v):
        if isinstance(v,dict):
            for x in v.values():walk(x)
        elif isinstance(v,(list,tuple)):
            for x in v:walk(x)
        elif isinstance(v,str)and len(v)<4096:
            try:
                p=Path(v)
                if p.is_file():
                    absolute=str(p.resolve())
                    # A report may mention the same large geometry bank in
                    # thousands of face witnesses. Bind each actual file once;
                    # verify_bindings still rechecks the completed snapshot.
                    if absolute not in found:found[absolute]=digest_file(p)
            except OSError:pass
    walk(value);return found


def result_file_bindings(result,input_files):
    """Retain independently audited artifacts even when they are also inputs.

    A tissue layout can be both a source binding and the literal buffer replay
    evidence used by the installer. Its actual audit, declared hash and file
    must agree; ordinary source inputs remain excluded from result artifacts.
    """
    found=file_bindings(result)
    input_paths={str(Path(p).resolve())for p in input_files}
    audited={str(Path(p).resolve())for p in result.get('audited_banks',[])}
    expected={str(Path(p).resolve()):sha for p,sha in result.get('proposal',{}).get('bank_bindings',{}).items()}
    for path in audited:
        if path not in found or path not in expected or found[path]!=expected[path]:
            raise ValueError('An independently audited artifact is missing or changed: '+path)
    return {p:sha for p,sha in found.items()if p not in input_paths or p in audited}


def verify_bindings(bindings):
    wrong=[]
    for path,sha in bindings.items():
        if not Path(path).is_file()or digest_file(path)!=sha:wrong.append(path)
    if wrong:raise ValueError('Bound input/artifact changed: '+', '.join(wrong[:5]))


def stopped_runtime(parameters):
    """Terminal supervisor evidence for this exact worker, never a stale PID."""
    path=parameters.get('runtime_report')
    if not path:return None
    try:
        runtime=json.loads(Path(path).read_text())
    except (OSError,ValueError):return None
    if runtime.get('pid')!=parameters.get('owner_pid')or runtime.get('state')not in ['timed_out','cancelled','failed','completed']:
        return None
    return dict(error='Worker '+runtime['state']+' before finalizing the job; inspect '+path,
                elapsed_seconds=runtime.get('elapsed_seconds'),runtime_report=path)


def row(db,jid):
    with db.connect()as c:r=c.execute('SELECT * FROM atlas_geom.jobs WHERE id=%s',(jid,)).fetchone()
    if not r:raise ValueError('Unknown geometry job')
    stopped=stopped_runtime(r['parameters'])if r['state']=='running'else None
    if stopped:
        with db.connect()as c:
            c.execute("UPDATE atlas_geom.jobs SET state='failed',error=%s,updated_at=now()WHERE id=%s AND state='running'",(stopped['error'],jid))
            r=c.execute('SELECT * FROM atlas_geom.jobs WHERE id=%s',(jid,)).fetchone()
    return r


def artifacts(db,jid):
    with db.connect()as c:items=c.execute('SELECT * FROM atlas_geom.job_artifacts WHERE job_id=%s',(jid,)).fetchall()
    verify_bindings({r['path']:r['sha256']for r in items});return items


def execute(db,spec,*,retry=False):
    payload=spec['parameters'];kind=spec['kind'];bindings=file_bindings(payload);mesh_bindings={};references=list(spec.get('meshes',[]));model=payload.get('model',spec.get('model'))
    for key in ['mesh','first','second','moving','protected']:
        value=payload.get(key,[]);references.extend(value if isinstance(value,list)else[value]if value else[])
    references=[x for x in references if isinstance(x,str)and not Path(x).is_file()]
    if kind=='fit'and payload.get('moving'):
        from .constraints import closure
        references.extend(x['id']for x in closure(db,payload['moving'],model)['members'])
        if payload.get('mode')=='icp':
            references.append(payload['registration_source'])
            target=payload['registration_target']
            if not Path(target).is_file():references.append(db.mesh(target,payload.get('target_model',model))['id'])
    from .registry import load_registered
    from .model import AtlasReader
    reader=AtlasReader()
    for reference in sorted(set(references)):
        m,mesh=load_registered(db,reference,model,reader=reader);mesh_bindings[m['id']]={'geometry':m['geometry_sha256'],'attributes':m['attributes_sha256']}
        source=m['source']
        if source.get('kind')=='npz':bindings[source['path']]=source['sha256']
        else:
            bindings[source['bank_path']]=source['bank_sha256']
            if source.get('catalogue_path'):bindings[source['catalogue_path']]=digest_file(source['catalogue_path'])
    constraint_bindings=[]
    if kind=='fit'and mesh_bindings:
        from .constraints import list_constraints
        constraint_bindings=list_constraints(db,mesh_ids=list(mesh_bindings))
        constraint_bindings=[{k:v for k,v in r.items()if k!='created_at'}for r in constraint_bindings]
    code=implementation_signature();signature=identity(kind,payload,bindings,mesh_bindings,constraint_bindings,code);jid=identity('geometry-job/v1',signature);attempt=1
    with db.connect()as c:
        c.execute('SELECT pg_advisory_xact_lock(hashtext(%s))',(jid,));old=c.execute('SELECT * FROM atlas_geom.jobs WHERE id=%s FOR UPDATE',(jid,)).fetchone()
        if old and old['state']in ['passed','reviewed','installed']:
            artifacts(db,jid);verify_bindings(bindings)
            return dict(job_id=jid,state=old['state'],cache_hit=True,result=old['result'])
        if old and old['state']=='running':
            stopped=stopped_runtime(old['parameters'])
            if not stopped:raise ValueError('Identical job is already running; inspect it instead of launching another solve')
            c.execute("UPDATE atlas_geom.jobs SET state='failed',error=%s,updated_at=now()WHERE id=%s",(stopped['error'],jid))
            old={**old,'state':'failed'}
        if old:
            attempt=old['parameters'].get('_attempt',1)+1
            if not retry or attempt>3:raise ValueError('Failed job requires --retry; at most three explicit attempts are allowed for the same inputs')
        params=dict(payload=payload,input_files=bindings,mesh_bindings=mesh_bindings,constraint_bindings=constraint_bindings,implementation_sha256=code,_attempt=attempt,owner_pid=os.getpid(),
                    runtime_report=os.environ.get('ATLAS_RUNTIME_REPORT'),runtime_deadline_utc=os.environ.get('ATLAS_RUNTIME_DEADLINE_UTC'))
        c.execute('''INSERT INTO atlas_geom.jobs(id,kind,input_signature,tool_version,parameters,state)VALUES(%s,%s,%s,%s,%s,'running')
          ON CONFLICT(id)DO UPDATE SET state='running',parameters=EXCLUDED.parameters,error=NULL,updated_at=now()''',(jid,kind,signature,VERSION,Jsonb(params)))
    output=CACHE/'jobs'/jid;output.mkdir(parents=True,exist_ok=True);began=time.monotonic()
    try:
        if kind=='audit':
            from .cli import load_ref
            from .audit import audit
            _,mesh=load_ref(db,payload['mesh'],payload.get('model'),payload.get('frame'));result=audit(mesh,**payload.get('controls',{}))
        elif kind=='contact':
            from .cli import load_ref
            from .contact import contacts
            _,a=load_ref(db,payload['first'],payload.get('model'),payload.get('frame'));_,b=load_ref(db,payload['second'],payload.get('model'),payload.get('frame'))
            _,result=contacts(db,a,b,payload['distance_m'],**payload.get('controls',{}));result['passed']=not result['truncated']
        elif kind=='fit':
            from .fit import fit
            result=fit(db,payload,output/('candidate-'+str(attempt)))
        elif kind=='catalogue-candidate':
            from .install import prepare
            proposal=prepare(db,{**payload,'job_id':jid},output/('prepared-'+str(attempt)))
            # These are actual independently decoded and audited output banks,
            # not arbitrary files promoted to reviewed artifacts.
            result=dict(passed=proposal['passed'],proposal=proposal,
                        audited_banks=list(proposal['bank_bindings']))
        elif kind=='tissue-layers':
            from .tissue import build_layers
            result=build_layers(payload['source'],output/('tissues-'+str(attempt)),payload['names'],payload.get('frame'),payload.get('origin',[0,0,0]),payload.get('max_pairs',250000),payload.get('max_tests',2000000))
        elif kind=='tissue-pack':
            from .tissue import pack
            import numpy as np
            results=[]
            for i,assembly in enumerate(payload.get('assemblies',[payload])):
                with np.load(assembly['source'],allow_pickle=False)as z:
                    positions=z['positions'].copy();surfaces={sid:z[key].copy()for sid,key in assembly['surfaces'].items()}
                    words={sid:z[key].copy()for sid,key in assembly['corner_words'].items()}if assembly.get('corner_words')else None
                    uvs={sid:z[key].copy()for sid,key in assembly['corner_uvs'].items()}if assembly.get('corner_uvs')else None
                results.append(pack(positions,surfaces,assembly['cells'],assembly['part_ids'],assembly['frame'],output/('tissues-'+str(attempt)+'-'+str(i)),origin=assembly.get('origin',[0,0,0]),corner_words=words,corner_uvs=uvs,max_pairs=assembly.get('max_pairs',1000000),max_tests=assembly.get('max_tests',2000000)))
            result=dict(passed=all(r['passed']for r in results),assemblies=results)
        elif kind in ['tube-path','tube-sweep','tube-join']:
            from .tubes import path,sweep,join
            result={'tube-path':path,'tube-sweep':sweep,'tube-join':join}[kind](payload,output/('tube-'+str(attempt)+'.npz'))
        else:raise ValueError('Unsupported job kind; arbitrary shell/Python execution is not a job type')
        verify_bindings(bindings)
        if implementation_signature()!=code:raise ValueError('Geometry implementation changed during this job; rerun against one consistent version')
        state='passed'if result.get('passed',True)else'failed'
        result=json.loads(json_bytes(result));result.update(elapsed_seconds=time.monotonic()-began)
        report=output/('result-'+str(attempt)+'.json');atomic_json(report,result)
        produced={str(p):digest_file(p)for p in output.rglob('*')if p.is_file()};produced.update(result_file_bindings(result,bindings))
        with db.connect()as c:
            for p,sha in produced.items():c.execute('INSERT INTO atlas_geom.job_artifacts VALUES(%s,%s,%s,%s)ON CONFLICT(job_id,role,path)DO UPDATE SET sha256=EXCLUDED.sha256',(jid,'result',p,sha))
            c.execute('UPDATE atlas_geom.jobs SET state=%s,result=%s,updated_at=now()WHERE id=%s',(state,Jsonb(result),jid))
        return dict(job_id=jid,state=state,cache_hit=False,passed=state=='passed',result=result)
    except BaseException as exc:
        with db.connect()as c:c.execute("UPDATE atlas_geom.jobs SET state='failed',error=%s,updated_at=now()WHERE id=%s",(str(exc),jid))
        raise


def review(db,jid,evidence_path):
    job=row(db,jid);artifacts(db,jid);verify_bindings(job['parameters']['input_files']);evidence=json.loads(Path(evidence_path).read_text())
    if job['state']not in ['passed','reviewed']:raise ValueError('Only a numerically passed job may receive review evidence')
    if evidence.get('input_signature')!=job['input_signature']or not evidence.get('reviewer')or evidence.get('accepted_for_local_preview')is not True:raise ValueError('Review must identify the current input signature, reviewer and explicit local-preview acceptance')
    views=evidence.get('views',[])
    if not views or not any(v.get('section')for v in views):raise ValueError('Actual surface and section views must be bound to the review')
    if not all(v.get('inspected')is True for v in views):raise ValueError('Uninspected renders are not review evidence')
    verify_bindings({str(Path(v['path']).resolve()):v['sha256']for v in views})
    if 'remaining_limits'not in evidence:raise ValueError('Record remaining limits, including partial anatomy')
    with db.connect()as c:
        c.execute('INSERT INTO atlas_geom.job_artifacts VALUES(%s,%s,%s,%s)ON CONFLICT DO NOTHING',(jid,'review',str(Path(evidence_path).resolve()),digest_file(evidence_path)))
        c.execute("UPDATE atlas_geom.jobs SET state='reviewed',updated_at=now()WHERE id=%s",(jid,))
    return dict(job_id=jid,state='reviewed',anatomical_acceptance=evidence.get('anatomical_acceptance',False),remaining_limits=evidence['remaining_limits'])


def add_parser(sub):
    p=sub.add_parser('job');s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('run');q.add_argument('spec');q.add_argument('--retry',action='store_true')
    q=s.add_parser('status');q.add_argument('id')
    q=s.add_parser('review');q.add_argument('id');q.add_argument('evidence')
    q=s.add_parser('recover');q.add_argument('id');q.add_argument('--reason',required=True)
    p.set_defaults(handler=run)


def run(args,db):
    if args.action=='run':return execute(db,json.loads(Path(args.spec).read_text()),retry=args.retry)
    if args.action=='review':return review(db,args.id,args.evidence)
    job=row(db,args.id)
    if args.action=='recover':
        if job['state']!='running':raise ValueError('Recovery applies only to an abandoned running job')
        try:os.kill(job['parameters']['owner_pid'],0)
        except ProcessLookupError:pass
        else:raise ValueError('The owner process still exists; do not recover or duplicate a live job')
        with db.connect()as c:c.execute("UPDATE atlas_geom.jobs SET state='failed',error=%s,updated_at=now()WHERE id=%s",('Recovered abandoned worker: '+args.reason,args.id))
        job=row(db,args.id)
    return job
