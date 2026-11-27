"""Render literal hand triangles against fixed bones and visible skin locally.

Flat diagnostic shading and a translucent skin overlay; this does not replace
the atlas material/packed-normal or anatomical review. No interactive app moves.
"""
import argparse
import json
import struct
import time
import zlib
from pathlib import Path
import numpy as np
from geometry_tools.course_io import sha,write_json
from geometry_tools.model import Mesh


def png(path,pixels):
    h,w,_=pixels.shape
    def chunk(kind,data): return struct.pack('!I',len(data))+kind+data+struct.pack('!I',zlib.crc32(kind+data)&0xffffffff)
    raw=b''.join(b'\0'+row.tobytes() for row in pixels)
    path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('!IIBBBBB',w,h,8,2,0,0,0))+chunk(b'IDAT',zlib.compress(raw))+chunk(b'IEND',b''))


def raster(triangles,colour,matrix,centre,scale,size):
    t=(triangles-centre)@matrix.T; xy=t[:,:,:2]*scale+size/2
    xy[:,:,1]=size-xy[:,:,1]
    image=np.zeros((size,size,3),float);depth=np.full((size,size),-np.inf)
    normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);length=np.linalg.norm(normal,axis=1)
    light=np.array([.3,.5,1]);light/=np.linalg.norm(light)
    shade=.55+.45*np.abs(normal@light/np.maximum(length,1e-20))
    for q,z,c in zip(xy,t[:,:,2],shade):
        low=np.maximum(np.floor(q.min(0)).astype(int),0);high=np.minimum(np.ceil(q.max(0)).astype(int),size-1)
        if np.any(low>high):continue
        a,b,cross=q;den=(b[1]-cross[1])*(a[0]-cross[0])+(cross[0]-b[0])*(a[1]-cross[1])
        if abs(den)<1e-12:continue
        yy,xx=np.mgrid[low[1]:high[1]+1,low[0]:high[0]+1];xx=xx+.5;yy=yy+.5
        u=((b[1]-cross[1])*(xx-cross[0])+(cross[0]-b[0])*(yy-cross[1]))/den
        v=((cross[1]-a[1])*(xx-cross[0])+(a[0]-cross[0])*(yy-cross[1]))/den;w=1-u-v
        zz=u*z[0]+v*z[1]+w*z[2];section=np.s_[low[1]:high[1]+1,low[0]:high[0]+1]
        mask=(u>=0)&(v>=0)&(w>=0)&(zz>depth[section]);depth[section][mask]=zz[mask]
        image[section][mask]=np.array(colour)*c
    return image,depth


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['left','right','bones','skin-reference','output']:parser.add_argument('--'+name,required=True)
    args=parser.parse_args();start=time.monotonic();out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new render stage')
    out.mkdir();bone_report=json.loads(Path(args.bones).read_text());sr=json.loads(Path(args.skin_reference).read_text())
    if sha(sr['path'])!=sr['sha256']:raise ValueError('Skin reference changed')
    skin=Mesh.load(sr['path']);panels=[];images=[]
    for side in ['left','right']:
        report_path=Path(getattr(args,side));data=json.loads(report_path.read_text());suffix='.l'if side=='left'else'.r'
        pieces=[]
        for r in data['rows']:
            if sha(r['path'])!=r['sha256']:raise ValueError('Candidate changed')
            pieces.append(Mesh.load(r['path']).triangles)
        arteries=np.concatenate(pieces)
        bones=np.concatenate([Mesh.load(r['path']).triangles for r in bone_report['rows']if r['sourceObject'].endswith(suffix)])
        bounds=np.array([bones.min((0,1)),bones.max((0,1))]);centre=bounds.mean(0)
        skin_tri=skin.triangles
        overlap=np.all(skin_tri.max(1)>=bounds[0]-.015,axis=1)&np.all(skin_tri.min(1)<=bounds[1]+.015,axis=1)
        skin_tri=skin_tri[overlap]
        for view,angle in [('front',0.),('back',np.pi),('oblique',.65)]:
            matrix=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
            projected=(bones-centre)@matrix.T;scale=450/np.ptp(projected[:,:,:2].reshape(-1,2),axis=0).max()
            bone_image,bd=raster(bones,[220,220,207],matrix,centre,scale,512)
            artery_image,ad=raster(arteries,[241,78,77],matrix,centre,scale,512)
            skin_image,sd=raster(skin_tri,[174,147,125],matrix,centre,scale,512)
            image=np.full((512,512,3),[18,26,33],float);mask=bd> -np.inf;image[mask]=bone_image[mask]
            mask=ad>bd;image[mask]=artery_image[mask]
            mask=(sd>np.maximum(ad,bd));image[mask]=image[mask]*.8+skin_image[mask]*.2
            image=np.clip(image,0,255).astype(np.uint8);path=out/f'{side}-{view}.png';png(path,image);panels.append(image)
            images.append(dict(side=side,view=view,path=str(path),sha256=sha(path),candidateSha256=sha(report_path)))
            if time.monotonic()-start>18:raise TimeoutError('18-second native hand rendering stage')
    sheet=np.concatenate([np.concatenate(panels[:3],axis=1),np.concatenate(panels[3:],axis=1)],axis=0)
    png(out/'contact-sheet.png',sheet);write_json(out/'report.json',dict(images=images,contactSheet=str(out/'contact-sheet.png'),
        contactSheetSha256=sha(out/'contact-sheet.png'),literalTriangles=True,skeletonChanged=False,
        skinOverlayOpacity=.2,diagnosticFlatShading=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start))
    print(json.dumps(dict(images=len(images),elapsedSeconds=time.monotonic()-start)))


if __name__=='__main__':main()
