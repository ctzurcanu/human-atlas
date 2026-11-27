"""Bounded surface/uncapped-section review of an encoded bladder point."""
from pathlib import Path
import argparse,json,time
import numpy as np
from PIL import Image,ImageDraw
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader

ap=argparse.ArgumentParser();ap.add_argument('catalogue');ap.add_argument('--output',required=True);args=ap.parse_args();started=time.monotonic()
path=Path(args.catalogue).resolve();out=Path(args.output).resolve();out.mkdir(exist_ok=True)
a=json.loads(path.read_text());m=a['finiteUrinaryIntegration'];sex=m['sex'];parts={p['id']:p for p in a['parts']};reader=AtlasReader()
g=next(g for g in m['termGroups']if g['term']=='A08.3.01.002');point=np.asarray(g['pointReference']['point'])
ids=next(g['partIds']for g in m['termGroups']if g['term']=='A08.3.01.001');tri=[];normal=[];color=[]
for pid in ids:
 p=parts[pid];mesh=reader.read_part(a,p,frame=f'model:local-{sex}-ta98:assembled-metres',catalogue_path=path)
 fs=np.array([f for q in p['tissueDrawGroups']if q.get('neighbor')not in ids for f in range(q['start']//3,(q['start']+q['count'])//3)],int)
 if not len(fs):continue
 words=mesh.attributes['normals'].astype(float)*mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices)))[:,None]
 tri.extend(mesh.triangles[fs]);normal.extend(words[mesh.faces[fs]])
 rgb=[218,140,125]if 'mucosa'in pid else[165,76,72]if any(s in pid for s in ['detrusor','trigone','longitudinal','circular'])else[215,175,154]
 color.extend([rgb]*len(fs))
tri=np.asarray(tri);normal=np.asarray(normal);color=np.asarray(color)
center=(tri.reshape(-1,3).min(0)+tri.reshape(-1,3).max(0))/2;span=np.ptp(tri.reshape(-1,3),axis=0).max()*1.15
records=[]
for name,axes,section in [('surface',(0,1,2),False),('sagittal',(2,1,0),True)]:
 ts=[];ns=[];cs=[]
 for t,n,c in zip(tri,normal,color):
  poly=list(zip(t,n))
  if section:
   poly=[];d=t[:,0]-point[0];inside=d<=0
   for i,j in [(0,1),(1,2),(2,0)]:
    if inside[i]:poly.append((t[i],n[i]))
    if inside[i]!=inside[j]:
     w=d[i]/(d[i]-d[j]);poly.append((t[i]+w*(t[j]-t[i]),n[i]+w*(n[j]-n[i])))
  for i in range(1,len(poly)-1):
   pick=[0,i,i+1];ts.append([poly[k][0]for k in pick]);ns.append([poly[k][1]for k in pick]);cs.append(c)
 ts=np.asarray(ts);ns=np.asarray(ns)[:,:,list(axes)];cs=np.asarray(cs);v=(ts-center)[:,:,list(axes)]
 screen=v[:,:,:2]*570/span;screen[:,:,0]+=320;screen[:,:,1]=345-screen[:,:,1];depth=v[:,:,2]
 zbuf=np.full((640,640),-np.inf);rgb=np.full((640,640,3),[19,27,35],np.uint8);light=np.array([-.3,.4,1.]);light/=np.linalg.norm(light)
 for k,t in enumerate(screen):
  lo=np.maximum(np.floor(t.min(0)).astype(int),[0,30]);hi=np.minimum(np.ceil(t.max(0)).astype(int),[639,639])
  if np.any(hi<lo):continue
  den=(t[1,1]-t[2,1])*(t[0,0]-t[2,0])+(t[2,0]-t[1,0])*(t[0,1]-t[2,1])
  if abs(den)<1e-18:continue
  x,y=np.meshgrid(np.arange(lo[0],hi[0]+1)+.5,np.arange(lo[1],hi[1]+1)+.5)
  w0=((t[1,1]-t[2,1])*(x-t[2,0])+(t[2,0]-t[1,0])*(y-t[2,1]))/den;w1=((t[2,1]-t[0,1])*(x-t[2,0])+(t[0,0]-t[2,0])*(y-t[2,1]))/den;w2=1-w0-w1
  zz=w0*depth[k,0]+w1*depth[k,1]+w2*depth[k,2];window=zbuf[lo[1]:hi[1]+1,lo[0]:hi[0]+1];mask=(w0>=0)&(w1>=0)&(w2>=0)&(zz>window)
  if not mask.any():continue
  nn=w0[:,:,None]*ns[k,0]+w1[:,:,None]*ns[k,1]+w2[:,:,None]*ns[k,2];nn/=np.maximum(np.linalg.norm(nn,axis=2)[:,:,None],1e-20)
  shade=.3+.7*np.maximum(np.einsum('ijk,k->ij',nn,light),0);pixels=np.clip(shade[:,:,None]*cs[k],0,255).astype(np.uint8)
  window[mask]=zz[mask];rgb[lo[1]:hi[1]+1,lo[0]:hi[0]+1][mask]=pixels[mask]
 image=Image.fromarray(rgb);draw=ImageDraw.Draw(image);anchor=(point-center)[list(axes)][:2]*570/span;sx,sy=float(anchor[0]+320),float(345-anchor[1])
 draw.line((sx,sy,480,70),fill=(148,212,222),width=2);draw.ellipse((sx-4,sy-4,sx+4,sy+4),fill=(148,212,222),outline=(25,50,60),width=2);draw.rounded_rectangle((472,55,630,85),radius=5,fill=(35,49,61),outline=(85,109,125));draw.text((484,65),'Apex of bladder',fill=(226,236,241))
 draw.text((12,12),sex+' / '+name+(' / diagnostic uncapped cut'if section else' / encoded normals'),fill='white')
 image_path=out/(name+'.png');image.save(image_path);records.append(dict(name=name,path=str(image_path),sha256=digest_file(image_path),section=section,triangles=len(ts)))
atomic_json(out/'render.json',dict(catalogue_sha256=digest_file(path),views=records,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(views=records,elapsed_seconds=time.monotonic()-started)),flush=True)
