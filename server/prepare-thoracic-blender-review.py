"""Queue literal surface/section packing views in the existing Blender."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.review import prepare
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--catalogue',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
path=ROOT/args.catalogue;a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres';out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False);m=a['finiteCardiacIntegration'];entries=[]
def read(pid):return r.read_part(a,parts[pid],frame=frame,catalogue_path=path)
heart_ids=list(dict.fromkeys(sum([m.get(k,[])for k in ['wallPartIds','valvePartIds','papillaryPartIds','chordalPartIds','rootPartIds']],[])))
for pid in heart_ids+a['pulmonarySegmentIntegration']['partIds']+a['pulmonaryBranchIntegration']['partIds']:
 p=parts[pid];mesh=read(pid);target=out/(str(len(entries))+'.npz');mesh.save(target);color=[.70,.30,.32,1.]if pid in heart_ids else[.62,.75,.83,1.]if pid in a['pulmonarySegmentIntegration']['partIds']else[.93,.73,.43,1.]
 entry=dict(mesh=str(target),part_id=pid,name=p['name'],color=color)
 if p.get('tissueDrawGroups'):entry['tissue_draw_groups']=p['tissueDrawGroups']
 entries.append(entry)
heart=read(f'TA98HEART:{args.sex}:epicardium');lo,hi=heart.bounds+np.array([[-.055]*3,[.055]*3]);hidden=set(a.get('defaultHidden',[]))
contexts=[('Native fixed rib cage',[p['id']for p in a['parts']if p['system']=='skeletal'and p['id']not in hidden and not p.get('defaultHidden')and not p.get('suppressed')and np.all(np.array(p['bounds'])[1]>=lo)and np.all(np.array(p['bounds'])[0]<=hi)],[.83,.81,.72,1.]),('Diaphragm',a['finiteDiaphragmIntegration']['partIds'],[.58,.40,.41,1.])]
for name,ids,color in contexts:
 meshes=[read(i)for i in ids];offsets=np.r_[0,np.cumsum([len(q.vertices)for q in meshes])];mesh=Mesh(np.concatenate([q.vertices for q in meshes]),np.concatenate([q.faces+offsets[i]for i,q in enumerate(meshes)]),frame);target=out/(str(len(entries))+'.npz');mesh.save(target);entries.append(dict(mesh=str(target),name=name,color=color))
centre=heart.vertices.mean(0);spec=dict(frame=frame,assembled_tissues=True,meshes=entries,resolution=1000,wall_seconds=60,views=[dict(name='surface',direction=[.35,.12,1.]),dict(name='section',direction=[0.,0.,1.],section=dict(normal=[0.,0.,1.],point=[0.,centre[1],.025],clear_outer=True)),dict(name='axial',direction=[0.,1.,0.],up=[0.,0.,-1.],section=dict(normal=[0.,1.,0.],point=[0.,centre[1],0.],clear_outer=True))]);atomic_json(out/'review-spec.json',spec);result=prepare(None,spec);atomic_json(out/'review-request.json',result);atomic_json(out/'source.json',dict(catalogue_sha256=digest_file(path),elapsed_seconds=time.monotonic()-start));print(json.dumps(result,indent=2),flush=True)
