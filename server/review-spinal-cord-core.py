"""Render literal encoded tissue contours with their holes retained."""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from PIL import Image,ImageDraw,ImageChops
from geometry_tools.common import atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('bank');ap.add_argument('--source');ap.add_argument('--encoded-stage',action='store_true');args=ap.parse_args();start=time.monotonic()
base=Path(args.bank).resolve();manifest=json.loads((base/('pending.json'if args.encoded_stage else'tissue.json')).read_text());raw=(base/'cells.bin').read_bytes()
out=base/'review';out.mkdir(exist_ok=True);origin=np.array([0.,1.35,-.05]);meshes={};families=defaultdict(list)
metadata=json.loads((Path(args.source)/'cells.json').read_text())if args.source else{}
family_of=lambda name:metadata.get('attributes',{}).get(name,{}).get('native_region',name.split('@')[0]).split('-')[0]
family_by_id={p['id']:family_of(name)for name,p in manifest['parts'].items()}
first=next(iter(manifest['parts'].values()));v=np.ndarray((first['vertexCount'],3),'<f4',buffer=raw,offset=first['positions']).astype(float)+origin
for name,p in manifest['parts'].items():
    assert p['positions']==first['positions']and p['vertexCount']==first['vertexCount']
    f=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices'])
    family=family_of(name)
    selected=[]
    for g in p['tissueDrawGroups']:
        neighbor=g.get('neighbor','')
        if family_by_id.get(neighbor)==family:continue
        selected.extend(range(g['start']//3,(g['start']+g['count'])//3))
    families[family].append(f[selected])
for family,faces in families.items():
    combined=np.concatenate(faces);ids=np.unique(combined);used,inv=np.unique(v[ids],axis=0,return_inverse=True)
    meshes[family]=(used,inv[np.searchsorted(ids,combined)])

def loops_at_y(v,f,y):
    fy=v[f,1];f=f[(fy.min(1)<y)&(fy.max(1)>y)]
    segments=[];points=[];cache={}
    def cross(a,b):
        key=tuple(sorted((int(a),int(b))))
        if key not in cache:
            p,q=v[list(key)];w=(y-p[1])/(q[1]-p[1]);cache[key]=len(points);points.append(p+w*(q-p))
        return cache[key]
    for tri in f:
        d=v[tri,1]-y;hits=[]
        for i,j in [(0,1),(1,2),(2,0)]:
            if d[i]*d[j]<0:hits.append(cross(tri[i],tri[j]))
        if len(hits)==2:segments.append(hits)
    edges=defaultdict(list)
    for a,b in segments:edges[a].append(b);edges[b].append(a)
    assert all(len(n)==2 for n in edges.values()),'Cut is not closed contours'
    todo=set(edges);loops=[]
    while todo:
        a=min(todo);loop=[];previous=None;current=a
        while current in todo:
            todo.remove(current);loop.append(current);peers=edges[current];following=next(n for n in peers if n!=previous);previous,current=current,following
        assert current==a;loops.append(np.asarray(points)[loop])
    return loops

colors={'white':(232,214,185),'grey':(158,122,115),'pia':(206,171,164)}
image=Image.new('RGB',(1200,750),(19,27,35));draw=ImageDraw.Draw(image);records=[]
for i,(label,y)in enumerate(zip(['Cervical','Thoracic','Lumbar enlargement','Conus'],[1.510123,1.350123,1.182123,1.153123])):
    contours={name:loops_at_y(v,f,y)for name,(v,f)in meshes.items()}
    allp=np.concatenate([p for loops in contours.values()for p in loops]);center=(allp[:,[0,2]].min(0)+allp[:,[0,2]].max(0))/2
    cx,cy=300+(i%2)*600,220+(i//2)*340;scale=32000
    for name in ['white','grey']:
        mask=Image.new('1',image.size,0)
        for poly in contours[name]:
            xy=(poly[:,[0,2]]-center)*scale;xy[:,0]+=cx;xy[:,1]=cy-xy[:,1]
            one=Image.new('1',image.size,0);ImageDraw.Draw(one).polygon([tuple(p)for p in xy],fill=1)
            mask=ImageChops.logical_xor(mask,one)
        image.paste(colors[name],mask=mask)
    draw.text((cx-180,cy-145),label+' / encoded tissue section',fill='white')
    draw.text((cx-180,cy+145),'Anterior above; central canal remains open',fill=(163,183,196))
    records.append(dict(name=label,plane_normal=[0,1,0],plane_point=[0,y,0],contours={k:len(p)for k,p in contours.items()},caps_are_tissue_only=True))
path=out/'encoded-cord-sections.png';image.save(path)
surface_image=Image.new('RGB',(1000,750),(19,27,35));sd=ImageDraw.Draw(surface_image)
surface_views=[]
for mode,cx,label in [('whole',240,'Fitted outer cord'),('grey',740,'Grey matter / thoracic interval')]:
    triangles=[];normals=[]
    for name,p in manifest['parts'].items():
        if mode=='grey'and family_of(name)!='grey':continue
        vv=v
        ff=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices'])
        ww=np.ndarray((p['vertexCount'],3),'<i2',buffer=raw,offset=p['normals'])
        signs=np.ndarray((p['vertexCount'],),'i1',buffer=raw,offset=p['normalSigns'])
        pick=[]
        for g in p['tissueDrawGroups']:
            peer=g.get('neighbor','')
            if peer and(mode=='whole'or family_by_id.get(peer)=='grey'):continue
            pick.extend(range(g['start']//3,(g['start']+g['count'])//3))
        selected=ff[pick]
        if mode=='grey':
            yy=vv[selected,1];selected=selected[(yy.min(1)<1.39)&(yy.max(1)>1.31)]
        for face in selected:
            t=vv[face]
            poly=list(t)
            if mode=='grey':
                for level,sign in [(1.31,1),(1.39,-1)]:
                    cut=[]
                    for pa,pb in zip(poly,poly[1:]+poly[:1]):
                        da,db=sign*(pa[1]-level),sign*(pb[1]-level)
                        if da>=0:cut.append(pa)
                        if da*db<0:cut.append(pa+(pb-pa)*da/(da-db))
                    poly=cut
            for j in range(1,len(poly)-1):triangles.append([poly[0],poly[j],poly[j+1]]);normals.append((ww[face].astype(float)*signs[face,None]).mean(0))
    ts=np.asarray(triangles);nn=np.asarray(normals);center=(ts.reshape(-1,3).min(0)+ts.reshape(-1,3).max(0))/2
    rotation=np.array([[.8,0.,-.6],[0.,1.,0.],[.6,0.,.8]])
    q=(ts-center)@rotation.T;nr=nn@rotation.T;nr/=np.maximum(np.linalg.norm(nr,axis=1)[:,None],1e-50)
    scale=570/np.ptp(ts.reshape(-1,3),axis=0).max();screen=q[:,:,:2]*scale;screen[:,:,0]+=cx;screen[:,:,1]=390-screen[:,:,1]
    rgb=np.array(colors['white'if mode=='whole'else'grey'])
    for k in np.argsort(q[:,:,2].mean(1)):
        light=.4+.6*max(0.,float(nr[k]@np.array([-.35,.2,.915])))
        sd.polygon([tuple(p)for p in screen[k]],fill=tuple((rgb*light).astype(int)))
    sd.text((cx-150,55),label,fill='white');surface_views.append(dict(mode=mode,triangles=len(ts),actual_stored_corner_normals=True))
surface_path=out/'encoded-cord-surfaces.png';surface_image.save(surface_path)
atomic_json(out/'surface-sections-clipped.json',dict(bank_sha256=digest_file(base/'cells.bin'),views=records,
    image_sha256=digest_file(path),surface_image_sha256=digest_file(surface_path),surface_views=surface_views,actual_float32_geometry=True,geometry_changed=False,
    anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(image=str(path),views=records,elapsed_seconds=time.monotonic()-start)),flush=True)
