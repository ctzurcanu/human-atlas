"""Render actual indexed tract-cell cross-sections for anatomical review."""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from PIL import Image,ImageDraw,ImageChops
from geometry_tools.common import atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);ap.add_argument('--ys',nargs='+',type=float);args=ap.parse_args()
started=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
d=json.loads((source/'cells.json').read_text())
with np.load(source/'cells.npz')as archive:z={key:archive[key]for key in archive.files}
v=z['positions']
meshes={}
for name,uses in d['cells'].items():
    meshes[name]=np.concatenate([z[u['surface']]if u['side']==1 else z[u['surface']][:,[0,2,1]]for u in uses])
def contours(f,y):
    q=v[f,1];f=f[(q.min(1)<y)&(q.max(1)>y)];points=[];cache={};edges=defaultdict(list)
    for face in f:
        hits=[]
        for a,b in [(0,1),(1,2),(2,0)]:
            i,j=sorted([int(face[a]),int(face[b])]);pa,pb=v[[i,j]]
            if (pa[1]-y)*(pb[1]-y)>=0:continue
            key=(i,j)
            if key not in cache:cache[key]=len(points);points.append(pa+(pb-pa)*(y-pa[1])/(pb[1]-pa[1]))
            hits.append(cache[key])
        if len(hits)==2:
            a,b=hits;edges[a].append(b);edges[b].append(a)
    assert all(len(n)==2 for n in edges.values()),'Open indexed contour'
    todo=set(edges);result=[]
    while todo:
        first=min(todo);prev=None;current=first;ids=[]
        while current in todo:
            todo.remove(current);ids.append(current);nxt=next(p for p in edges[current]if p!=prev);prev,current=current,nxt
        assert current==first;result.append(np.asarray(points)[ids])
    return result

colors=dict(gracile='#ddba55',cuneate='#54b8bc',**{'lateral-corticospinal':'#f06d67','anterior-spinocerebellar':'#a38bc6',
    'rubrospinal':'#91be62','spino-olivary':'#e3a2bd','lateral-vestibulospinal':'#779fdf',
    'anterior-corticospinal':'#e29245','tectospinal':'#8dd8c0',
    'bulboreticulospinal':'#a2c7e6','pontoreticulospinal':'#c9a4ed','reticulospinal-fibres':'#779da1'})
regional_colors={'marginal':'#eab463','gelatinous':'#bda4df','proprius':'#668fcf',
    'lamina-V':'#ce879d','lamina-VI':'#6bbaa6','grey-commissure-anterior':'#d9ad62',
    'grey-commissure-posterior':'#d3c67c','phrenic-nucleus':'#a4c658',
    'accessory-nucleus':'#de9672','intermediomedial':'#abbd75',
    'intermediolateral':'#d48f57','clarke':'#a099d0','sacral-parasympathetic':'#ba9ace',
    'lateral-horn':'#bb8d82','lateral-intermediate':'#b99c7e',
    'anterior-proprius':'#a9d1b3','lateral-proprius':'#a9d1b3','posterior-proprius':'#a9d1b3',
    'anteromedial-motor':'#c9a5ed','posteromedial-motor':'#ebb57f','anterolateral-motor':'#f3909e',
    'posterolateral-motor':'#94c3d7','retroposterolateral-motor':'#e9cf72','central-motor':'#8acaa3',
    'posterior-spinocerebellar':'#c6d787','lissauer-posterior':'#f4a767','lissauer-lateral':'#f4a767',
    'white-commissure-anterior-anterior':'#cddfe6','white-commissure-anterior-lateral':'#cddfe6',
    'white-commissure-posterior-posterior':'#bdd6e2','white-commissure-posterior-lateral':'#bdd6e2'}
used=np.unique(np.concatenate(list(meshes.values())));low,high=v[used,1].min(),v[used,1].max()
ys=np.array(args.ys if args.ys else[low+(high-low)*w for w in [.22,.5,.78]])+1.237e-7
image=Image.new('RGB',(max(3,len(ys))*500,650),'#131b23');draw=ImageDraw.Draw(image);records=[]
for view,y in enumerate(ys):
    loops={name:contours(f,float(y))for name,f in meshes.items()};allp=np.concatenate([p for values in loops.values()for p in values])
    lo,hi=allp[:,[0,2]].min(0),allp[:,[0,2]].max(0);center=(lo+hi)/2;scale=380/max(hi-lo)
    cx,cy=250+view*500,275
    for name,polys in sorted(loops.items(),key=lambda r:r[0].startswith('tract-')):
        key=name.split('@')[0].split('~')[0];native=d.get('attributes',{}).get(name,{}).get('native_region',key)
        if key.startswith('tract-'):color=colors[key.removeprefix('tract-').rsplit('-',1)[0]]
        elif key.rsplit('-',1)[0]in regional_colors:color=regional_colors[key.rsplit('-',1)[0]]
        elif native.startswith('pia'):color='#d6a79b'
        elif native.startswith('grey'):color='#a77d7d'
        elif native.startswith('white'):color='#e5d6bc'
        else:raise ValueError('Unknown anatomical review colour: '+name)
        mask=Image.new('1',image.size,0)
        for poly in polys:
            xy=(poly[:,[0,2]]-center)*scale;xy[:,0]+=cx;xy[:,1]=cy-xy[:,1]
            one=Image.new('1',image.size,0);ImageDraw.Draw(one).polygon([tuple(p)for p in xy],fill=1);mask=ImageChops.logical_xor(mask,one)
        image.paste(color,mask=mask)
    draw.text((cx-210,35),f"{d['segment']} / actual tissue cells / y={y:.6f} m",fill='white')
    draw.text((cx-210,480),'Anterior above / left cord on right of image',fill='#b9cad6')
    records.append(dict(y_m=float(y),cell_contours={n:len(p)for n,p in loops.items()},actual_indexed_geometry=True))
for i,(name,color)in enumerate(colors.items()):
    x=30+(i%3)*490;y=535+(i//3)*28;draw.rectangle([x,y,x+14,y+14],fill=color);draw.text((x+24,y),name,fill='white')
path=out/'tract-sections.png';image.save(path)
atomic_json(out/'review.json',dict(source_sha256=digest_file(source/'cells.npz'),image_sha256=digest_file(path),
    views=records,probabilistic_template_registration=True,anatomical_acceptance=False,installed=False,
    elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(image=str(path),segment=d['segment'],views=len(records),elapsed_seconds=time.monotonic()-started)),flush=True)
