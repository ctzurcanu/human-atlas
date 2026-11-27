"""Hash-bound review requests for the already running Blender session."""
from pathlib import Path
import json,time
from .common import CACHE,ROOT,atomic_json,atomic_npz,identity,digest_file
from .cli import load_ref
from .model import Mesh
import numpy as np

QUEUE=CACHE/'review'


def prepare(db,spec):
    if not spec.get('meshes'):raise ValueError('Specify the actual meshes to review')
    if len(spec['meshes'])>100:raise ValueError('Review batches are bounded to 100 meshes')
    nonce=identity(spec,time.time_ns());folder=QUEUE/nonce;folder.mkdir(parents=True);rows=[];frames=set()
    visible={e['part_id']for e in spec['meshes']if isinstance(e,dict)and e.get('part_id')}
    for i,entry in enumerate(spec['meshes']):
        entry={'mesh':entry}if isinstance(entry,str)else entry
        _,mesh=load_ref(db,entry['mesh'],spec.get('model'),spec.get('frame'));frames.add(mesh.frame)
        corner_words=None;corner_source=entry.get('corner_normal_source')
        if corner_source is not None:
            cp=Path(corner_source['path']).resolve()
            if not cp.is_relative_to(ROOT)or digest_file(cp)!=corner_source['sha256']:
                raise ValueError('Corner-normal source must be local and hash-bound')
            with np.load(cp,allow_pickle=False)as z:corner_words=z[corner_source.get('array','corner_normals')].copy()
            if corner_words.shape!=(len(mesh.faces),3,3)or not np.isfinite(corner_words).all()or np.any(np.linalg.norm(corner_words,axis=2)<=0):
                raise ValueError('One finite nonzero normal word is required per actual face corner')
        removed=0
        if spec.get('assembled_tissues')and entry.get('tissue_draw_groups')is not None:
            groups=entry['tissue_draw_groups'];end=0;keep=[]
            if not entry.get('part_id'):
                raise ValueError('Assembled tissue review needs explicit part identities')
            for g in groups:
                if g['start']!=end or g['count']<=0 or g['count']%3:
                    raise ValueError('Review tissue draw groups must cover actual consecutive triangles')
                end+=g['count']
                if g.get('neighbor')in visible:removed+=g['count']//3
                else:keep.extend(range(g['start']//3,end//3))
            if end!=mesh.faces.size:
                raise ValueError('Review tissue draw groups do not cover the source mesh')
            if not keep:
                continue
            if corner_words is not None:corner_words=corner_words[np.asarray(keep,int)]
            mesh=Mesh(mesh.vertices,mesh.faces[np.array(keep,int)],mesh.frame,mesh.attributes,
                      dict(mesh.source,diagnostic_hidden_reciprocal_interfaces=removed))
        # Shared construction banks may contain a long ureter or other cells.
        # Review framing must use this displayed mesh's actual indexed extent.
        used,inv=np.unique(mesh.faces,return_inverse=True)
        mesh=Mesh(mesh.vertices[used],inv.reshape(-1,3),mesh.frame,{k:v[used]for k,v in mesh.attributes.items()},mesh.source)
        p=folder/(str(i)+'.npz')
        if corner_words is None:mesh.save(p)
        else:
            atomic_npz(p,vertices=mesh.vertices,faces=mesh.faces,frame=np.array(mesh.frame),
                source_json=np.array(json.dumps(mesh.source,sort_keys=True)),
                **{'attribute_'+k:v for k,v in mesh.attributes.items()},corner_normals=corner_words)
        row=dict(path=str(p),sha256=digest_file(p),geometry_sha256=mesh.geometry_sha,frame=mesh.frame,name=entry.get('name',entry['mesh']),color=entry.get('color',[.68,.42,.38,1.]),diagnostic_hidden_reciprocal_interfaces=removed)
        if entry.get('part_id'):row['part_id']=entry['part_id']
        if entry.get('nonanatomical_framing_only')or mesh.source.get('nonanatomical_framing_only'):
            row['nonanatomical_framing_only']=True
        rows.append(row)
    if not rows:raise ValueError('Assembled review has no exposed source surfaces')
    if len(frames)!=1:raise ValueError('Review context must use one explicitly registered frame')
    deadline=float(spec.get('wall_seconds',60))
    if not 1<=deadline<=120:raise ValueError('A review needs a 1–120 second native deadline')
    task=dict(schema='blender-review/v1',id=nonce,action='render',meshes=rows,views=spec.get('views',[{'name':'front','direction':[0,0,1]},{'name':'right','direction':[1,0,0]}]),resolution=spec.get('resolution',640),output=str(folder),frame=next(iter(frames)),wall_seconds=deadline)
    atomic_json(folder/'request.json',task);sha=digest_file(folder/'request.json');atomic_json(QUEUE/'pending'/(nonce+'.json'),dict(request=str(folder/'request.json'),sha256=sha))
    script=Path(__file__).with_name('blender_bridge.py')
    return dict(request_id=nonce,folder=str(folder),request_sha256=sha,state='queued',bootstrap_python="atlas_geometry_review_request_id="+repr(nonce)+"; exec(compile(open("+repr(str(script))+").read(),'atlas-geometry-review-bridge','exec'))",existing_Blender_only=True,main_file_save=False)


def status(request_id):
    if len(request_id)!=64 or any(c not in '0123456789abcdef'for c in request_id):raise ValueError('Expected a review request SHA identity')
    folder=QUEUE/request_id;result=folder/'result.json'
    if result.exists():
        r=json.loads(result.read_text())
        for image in r.get('images',[]):
            if digest_file(image['path'])!=image['sha256']:raise ValueError('Review image changed after rendering')
        return r
    return dict(request_id=request_id,state='queued',bridge_status=json.loads((QUEUE/'bridge.json').read_text())if(QUEUE/'bridge.json').exists()else None)


def add_parser(sub):
    p=sub.add_parser('review');s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('prepare');q.add_argument('spec')
    q=s.add_parser('status');q.add_argument('id')
    p.set_defaults(handler=lambda args,db:prepare(db,json.loads(Path(args.spec).read_text()))if args.action=='prepare'else status(args.id))
