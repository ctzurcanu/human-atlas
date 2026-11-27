"""Actual source triangles and literal transverse contours; no inferred caps."""
from pathlib import Path
import argparse,json,time
import numpy as np
from PIL import Image,ImageDraw
from geometry_tools.common import atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('stage');args=ap.parse_args();start=time.monotonic()
base=Path(args.stage).resolve();z=np.load(base/'surface.npz');v=z['positions'];f=z['faces'];tri=v[f]
out=base/'review';out.mkdir(exist_ok=True)
image=Image.new('RGB',(1200,750),(19,27,35));d=ImageDraw.Draw(image)
rows=[]
ys=[1.510123,1.350123,1.182123,1.153123]
for i,(name,y)in enumerate(zip(['Cervical','Thoracic','Lumbar enlargement','Conus'],ys)):
    segments=[]
    for t in tri:
        values=t[:,1]-y;hits=[]
        for a,b in [(0,1),(1,2),(2,0)]:
            if values[a]*values[b]<0:
                w=values[a]/(values[a]-values[b]);hits.append(t[a]+w*(t[b]-t[a]))
        if len(hits)==2:segments.append(hits)
    points=np.asarray(segments);xz=points[:,:,[0,2]];center=(xz.reshape(-1,2).min(0)+xz.reshape(-1,2).max(0))/2
    scale=200/.017;cx=450+(i%2)*390;cy=230+(i//2)*330
    for pair in xz:
        screen=(pair-center)*scale;screen[:,0]+=cx;screen[:,1]=cy-screen[:,1]
        d.line([tuple(q)for q in screen],fill=(233,200,91),width=3)
    d.text((cx-130,cy-120),name+' / y='+str(y)+'m',fill='white')
    d.text((cx-130,cy+120),'Anterior (+Z) above; posterior below',fill=(163,183,196))
    d.line((cx-100,cy+100,cx-100+10*scale/1000,cy+100),fill=(163,183,196),width=3)
    rows.append(dict(name=name,y_m=y,segments=len(segments),contour_segments_m=points.tolist(),width_mm=float(np.ptp(xz[:,:,0])*1000),ap_mm=float(np.ptp(xz[:,:,1])*1000)))
for axes,cx,cy,label in [((0,1,2),130,390,'Posterior'),((2,1,0),260,390,'Lateral')]:
    q=(tri-v.mean(0))[:,:,list(axes)];s=q[:,:,:2]*1350;s[:,:,0]+=cx;s[:,:,1]=cy-s[:,:,1]
    normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-30)
    for k in np.argsort(q[:,:,2].mean(1)):
        light=.45+.55*abs(float(normal[k]@np.array([.25,.1,.96])))
        color=tuple(int(t*light)for t in [230,199,96]);d.polygon([tuple(p)for p in s[k]],fill=color)
    d.text((cx-40,65),label,fill='white')
image.save(out/'native-cord-surface-sections.png')
atomic_json(out/'surface-sections.json',dict(surface_sha256=digest_file(base/'surface.npz'),sections=rows,
    image_sha256=digest_file(out/'native-cord-surface-sections.png'),caps_created=0,geometry_changed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(image=str(out/'native-cord-surface-sections.png'),sections=[{k:r[k]for k in ['name','width_mm','ap_mm','segments']}for r in rows],elapsed_seconds=time.monotonic()-start)))
