"""Render actual emitted finite coats, cavity contours and fixed native cord."""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from PIL import Image,ImageDraw,ImageChops
from geometry_tools.common import ROOT,atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('bank');args=ap.parse_args();start=time.monotonic();bank=Path(args.bank).resolve();out=bank/'review';out.mkdir(exist_ok=True)
record=json.loads((bank/'tissue.json').read_text());assert record['passed'];raw=(bank/'cells.bin').read_bytes();meshes={}
for name,p in record['parts'].items():
    v=np.ndarray((p['vertexCount'],3),'<f4',buffer=raw,offset=p['positions']).astype(float)+[0.,1.35,-.05]
    f=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices']);v,mapping=np.unique(v,axis=0,return_inverse=True);meshes[name]=(v,mapping[f])
cord=np.load(bank.parent/'female-surface-v1/surface.npz');vv,mapping=np.unique(cord['positions'],axis=0,return_inverse=True);meshes['cord']=(vv,mapping[cord['faces']])
colors={'dura':(184,167,145),'arachnoid':(185,204,216),'cord':(224,200,162)}
image=Image.new('RGB',(1600,1120),(19,27,35));draw=ImageDraw.Draw(image)
# A literal half surface reveals the open subarachnoid space without a
# fabricated cap across it. Normal shading uses the actual emitted triangles.
for panel,(h,depth,label)in enumerate([(0,2,'Posterior / right half cut away'),(2,0,'Lateral / open sac and native cord')]):
    tris=[]
    for name,(v,f)in meshes.items():
        q=v[f]
        if name!='cord':q=q[q[:,:,0].max(1)<=.0003]
        n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);n/=np.maximum(np.linalg.norm(n,axis=1)[:,None],1e-100)
        for t,nn in zip(q,n):
            light=.45+.55*abs(float(nn@np.array([.35,.12,.93])))
            tris.append((t,float(t[:,depth].mean()),tuple(int(x*light)for x in colors[name])))
    for t,_,color in sorted(tris,key=lambda r:r[1]):
        x=(t[:,h]+(.05 if h==2 else 0))*1400+220+panel*460;y=940-(t[:,1]-.97)*1400
        draw.polygon(list(zip(x,y)),fill=color)
    draw.text((30+panel*460,35),label,fill='white')

def sections(v,f,y):
    points=[];edge_cache={};segments=[]
    for tri in f:
        dd=v[tri,1]-y;hits=[]
        for j,k in [(0,1),(1,2),(2,0)]:
            if dd[j]*dd[k]<0:
                key=tuple(sorted((int(tri[j]),int(tri[k]))))
                if key not in edge_cache:
                    a,b=v[list(key)];edge_cache[key]=len(points);points.append(a+(b-a)*(y-a[1])/(b[1]-a[1]))
                hits.append(edge_cache[key])
        if len(hits)==2:segments.append(hits)
    edges=defaultdict(list)
    for a,b in segments:edges[a].append(b);edges[b].append(a)
    assert all(len(x)==2 for x in edges.values()),'Open or branched encoded section'
    todo=set(edges);loops=[]
    while todo:
        first=min(todo);prev=None;cur=first;loop=[]
        while cur in todo:
            todo.remove(cur);loop.append(cur);nxt=next(n for n in edges[cur]if n!=prev);prev,cur=cur,nxt
        assert cur==first;loops.append(np.asarray(points)[loop])
    return loops
views=[]
for k,(name,y)in enumerate([('Cervical',1.510123),('Thoracic',1.350123),('Lumbar cord',1.182123),('Caudal dural sac',1.050123)]):
    contours={n:sections(v,f,y)for n,(v,f)in meshes.items()};allp=np.concatenate([p for ls in contours.values()for p in ls]);center=(allp.min(0)+allp.max(0))/2
    cx,cy=1210,165+k*245;scale=15500
    for family,loops in contours.items():
        mask=Image.new('1',image.size,0)
        for p in loops:
            q=(p[:,[0,2]]-center[[0,2]])*scale;q[:,0]+=cx;q[:,1]=cy-q[:,1]
            one=Image.new('1',image.size,0);ImageDraw.Draw(one).polygon(list(map(tuple,q)),fill=1);mask=ImageChops.logical_xor(mask,one)
        image.paste(colors[family],mask=mask)
    draw.text((975,cy-112),name+' / actual tissue cut',fill='white');views.append(dict(name=name,point=[0,y,0],normal=[0,1,0],contours={n:len(ls)for n,ls in contours.items()}))
draw.text((30,1040),'Beige: dura. Blue: arachnoid. Native cord: gold. Dark spaces are not filled tissue.',fill=(190,205,216))
path=out/'encoded-coverings-surfaces-sections.png';image.save(path)
atomic_json(out/'review.json',dict(bank_sha256=digest_file(bank/'cells.bin'),image=str(path),image_sha256=digest_file(path),
    actual_encoded_geometry=True,views=views,diagnostic_half_surface_is_uncapped=True,skeleton_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(image=str(path),elapsed_seconds=time.monotonic()-start)),flush=True)
