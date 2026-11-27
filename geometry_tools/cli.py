"""Local geometry workbench. Coordinates and distances are metres."""
from pathlib import Path
import argparse, json, sys, traceback
import numpy as np
from .common import atomic_json
from .db import Database
from .model import Mesh
from .registry import import_catalogue, load_registered, search,inventory
from .spatial import cached_bvh, neighbors
from .contact import contacts


def serial(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, Path): return str(value)
    if hasattr(value,'isoformat'): return value.isoformat()
    raise TypeError(type(value).__name__)


def emit(value, report=None):
    clean=json.loads(json.dumps(value,default=serial,allow_nan=False))
    if report: atomic_json(report,clean)
    printed=clean
    if report and len(json.dumps(clean))>10000:
        printed={k:v for k,v in clean.items()if k in ['schema','model','model_id','passed','state','parts','errors','term_count','elapsed_seconds','cache_hit','cache_key']}
        printed.update(report=str(Path(report).resolve()),full_details_in_report=True)
    print(json.dumps(printed,indent=2,allow_nan=False),flush=True)
    return clean


def load_ref(db, value, model=None, frame=None):
    if Path(value).suffix=='.npz' and Path(value).is_file(): return None,Mesh.load(value,frame=frame)
    return load_registered(db,value,model)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dsn',help='Local PostgreSQL connection string; prefer ATLAS_GEOMETRY_DSN')
    p.add_argument('--report',help='New JSON evidence file (never overwrites)')
    s=p.add_subparsers(dest='command',required=True)
    s.add_parser('db-init');s.add_parser('status')
    q=s.add_parser('registry-import');q.add_argument('catalogue');q.add_argument('--model',required=True);q.add_argument('--frame');q.add_argument('--root');q.add_argument('--max-parts',type=int,help='Save at most this many records; resume the exact revision in a later bounded call. Publish only when complete')
    q=s.add_parser('register');q.add_argument('mesh');q.add_argument('--part-id',required=True);q.add_argument('--name',required=True);q.add_argument('--frame')
    q=s.add_parser('export');q.add_argument('mesh');q.add_argument('--model');q.add_argument('--output',required=True)
    q=s.add_parser('inventory');q.add_argument('root_term');q.add_argument('--model',required=True);q.add_argument('--authority');q.add_argument('--classifications');q.add_argument('--sex',choices=['male','female'])
    q=s.add_parser('search');q.add_argument('--query');q.add_argument('--similar');q.add_argument('--model');q.add_argument('--sex');q.add_argument('--side');q.add_argument('--term');q.add_argument('--system');q.add_argument('--license',dest='license_name');q.add_argument('--include-hidden',action='store_true');q.add_argument('--limit',type=int,default=20)
    q=s.add_parser('index');q.add_argument('mesh');q.add_argument('--model');q.add_argument('--frame');q.add_argument('--leaf-size',type=int,default=12)
    q=s.add_parser('neighbors');q.add_argument('mesh');q.add_argument('--model');q.add_argument('--distance',type=float,default=.002);q.add_argument('--include-hidden',action='store_true');q.add_argument('--limit',type=int,default=1000)
    q=s.add_parser('closest');q.add_argument('first');q.add_argument('second');q.add_argument('--model');q.add_argument('--frame')
    q=s.add_parser('contact');q.add_argument('first');q.add_argument('second');q.add_argument('--model');q.add_argument('--frame');q.add_argument('--distance',type=float,required=True);q.add_argument('--self',dest='selfcheck',action='store_true');q.add_argument('--max-pairs',type=int,default=250000);q.add_argument('--area-error',type=float);q.add_argument('--area-max-depth',type=int,default=8)
    # Additional tools register their parsers here, avoiding imports of Blender
    # or optional construction engines for ordinary read-only searches.
    for module in ['audit','constraints','fit','tissue','tubes','jobs','review','install','open3d_backend']:
        try:
            imported=__import__('geometry_tools.'+module,fromlist=['add_parser'])
        except ModuleNotFoundError as exc:
            if exc.name=='geometry_tools.'+module:continue
            raise
        imported.add_parser(s)
    return p


def main(argv=None):
    args=parser().parse_args(argv);db=Database(args.dsn);c=args.command
    try:
        if c=='db-init': out=db.initialize()
        elif c=='status':out=db.status()
        elif c=='registry-import':out=import_catalogue(db,args.catalogue,model_id=args.model,frame_id=args.frame,root=args.root,max_parts=args.max_parts,progress=lambda v:print(json.dumps(v),file=sys.stderr,flush=True))
        elif c=='register':
            m=Mesh.load(args.mesh,frame=args.frame);out={'id':db.register_artifact_mesh(m,args.mesh,part_id=args.part_id,name=args.name),'geometry_sha256':m.geometry_sha}
        elif c=='export':
            row,m=load_registered(db,args.mesh,args.model);m.save(args.output);out={'id':row['id'],'geometry_sha256':m.geometry_sha,'output':args.output}
        elif c=='search':out=search(db,**{k:getattr(args,k)for k in ['query','similar','model','sex','side','term','system','license_name','include_hidden','limit']})
        elif c=='inventory':out=inventory(db,args.root_term,args.model,authority=args.authority,classifications=json.loads(Path(args.classifications).read_text())if args.classifications else None,sex=args.sex)
        elif c=='index':
            row,m=load_ref(db,args.mesh,args.model,args.frame);b,hit=cached_bvh(db,m,args.leaf_size);out={**hit,'triangles':len(m.faces),'nodes':len(b.low),'geometry_sha256':m.geometry_sha}
        elif c=='neighbors':out=neighbors(db,args.mesh,model=args.model,distance_m=args.distance,include_hidden=args.include_hidden,limit=args.limit)
        elif c=='closest':
            from .contact import require_common_frame
            ra,a=load_ref(db,args.first,args.model,args.frame);rb,b=load_ref(db,args.second,args.model,args.frame);require_common_frame(a,b)
            ba,ca=cached_bvh(db,a);bb,cb=cached_bvh(db,b);out={**ba.closest_mesh(bb),'bvh_cache':[ca,cb],'frame':a.frame,'geometry_sha256':[a.geometry_sha,b.geometry_sha],'anatomical_attachment_inferred':False}
        elif c=='contact':
            ra,a=load_ref(db,args.first,args.model,args.frame);rb,b=load_ref(db,args.second,args.model,args.frame)
            arrays,out=contacts(db,a,b,args.distance,selfcheck=args.selfcheck,max_pairs=args.max_pairs,area_error_m2=args.area_error,area_max_depth=args.area_max_depth)
            out={k:v for k,v in out.items()if not k.endswith('_patch_components')}
            out.update(first=args.first,second=args.second,geometry_sha256=[a.geometry_sha,b.geometry_sha],witness_preview={k:v[:5]for k,v in arrays.items()})
        else:out=args.handler(args,db)
        emit(out,args.report)
        return 2 if out.get('passed')is False or out.get('errors') else 0
    except (ValueError,OSError,KeyError,RuntimeError) as exc:
        emit({'error':type(exc).__name__,'message':str(exc),'passed':False},args.report)
        return 2
