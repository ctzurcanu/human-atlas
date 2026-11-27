"""Render actual native covering/root surfaces against the fixed skeleton."""
from pathlib import Path
import argparse,json,time
import numpy as np
from PIL import Image,ImageDraw
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader

ap=argparse.ArgumentParser();ap.add_argument('stage');ap.add_argument('--dura-candidate');args=ap.parse_args();start=time.monotonic();stage=Path(args.stage).resolve()
data=json.loads((stage/'sources.json').read_text());a=json.loads(Path(data['source_catalogue']).read_text());reader=AtlasReader();meshes=[]
for row in data['sources']:
    name=row['name'];z=np.load(args.dura_candidate if args.dura_candidate and 'dura'in name.lower()else row['mesh']);v=z['vertices'];f=z['faces']
    if 'dura'in name.lower():color,alpha=(188,170,156),80
    elif 'Cauda'in name:color,alpha=(238,193,88),255
    elif 'Lumbar'in name:color,alpha=(201,144,58),220
    elif 'Sacral'in name:color,alpha=(230,174,75),220
    else:continue
    meshes.append((name,v,f,color,alpha))
for p in a['parts']:
    if p.get('suppressed')or not p['id'].startswith('LOCAL:female:'):continue
    if 'lumbar_vertebra_'not in p['id']and not p['id'].startswith('LOCAL:female:sacrum:'):continue
    m=reader.read_part(a,p,frame='model:local-female-ta98:assembled-metres',catalogue_path=data['source_catalogue'])
    meshes.append((p['name'],m.vertices,m.faces,(188,195,196),70))
surface=np.load(stage.parent/'female-surface-v1/surface.npz');meshes.append(('Cord',surface['positions'],surface['faces'],(239,211,164),255))
im=Image.new('RGB',(1400,1000),(19,27,35));draw=ImageDraw.Draw(im,'RGBA')
for panel,(horizontal,depth,title)in enumerate([(0,2,'Posterior: native roots and dural termination'),(2,0,'Lateral: fixed lumbar spine and sacrum')]):
    all_tri=[]
    for name,v,f,color,alpha in meshes:
        tri=v[f];keep=(tri[:,:,1].max(1)>.93)&(tri[:,:,1].min(1)<1.225);tri=tri[keep]
        normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-100)
        for t,n in zip(tri,normal):
            light=.45+.55*abs(float(n@np.array([.2,.1,.97])))
            all_tri.append((t,float(t[:,depth].mean()),tuple(int(x*light)for x in color)+(alpha,)))
    for t,_,color in sorted(all_tri,key=lambda t:t[1]):
        x=t[:,horizontal]*(2700 if panel==0 else 3300)+(350 if panel==0 else 1160)
        y=900-(t[:,1]-.93)*2700
        draw.polygon(list(zip(x,y)),fill=color)
    draw.text((30+panel*700,30),title,fill=(235,240,244,255))
    draw.text((30+panel*700,60),'Actual source surfaces; no invented section caps',fill=(167,185,196,255))
    y=900-(1.1471889-.93)*2700;draw.line((80+panel*700,y,620+panel*700,y),fill=(157,185,200,130),width=1)
    draw.text((40+panel*700,y+5),'Conus / native dural end',fill=(167,185,196,255))
render_stage=Path(args.dura_candidate).resolve().parent if args.dura_candidate else stage
path=render_stage/'native-coverings-fixed-spine.png';im.save(path)
atomic_json(render_stage/'render.json',dict(image=str(path),sha256=digest_file(path),source_sha256=digest_file(stage/'sources.json'),
    actual_source_geometry=True,section_caps_added=0,skeleton_changed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(image=str(path),elapsed_seconds=time.monotonic()-start)))
