"""Render actual proposedbank ownership pieces and literal source context."""
from pathlib import Path
import sys,json,numpy as np
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
B=ROOT/".local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1";O=B/"native-terminal-ownership-v12";A=json.loads((O/"proposed-atlas.json").read_text());R=AtlasReader();meshes={};colors={};frame="male-native-world-metres"
for side in["l","r"]:
 for suffix,k,c in[("terminal-ICA","stem",(218,144,24)),("A12.2.07.047","m1",(196,44,66))]:
  pid=f"TA98PART:ZA:Middle cerebral artery (M1-segment).{side}:{suffix}";p=next(p for p in A["parts"]if p["id"]==pid);meshes[f"{k}-{side}"]=R.read_part(A,p,frame=frame);colors[f"{k}-{side}"]=c
 for stem,pid,c in[("ica",f"TA98PART:ZA:Internal carotid artery.{side}:A12.2.06.015",(48,109,178)),("aca",f"ZA:Anterior cerebral artery.{side}",(50,141,85))]:meshes[f"{stem}-{side}"]=R.read_part(A,next(p for p in A["parts"]if p["id"]==pid),frame=frame);colors[f"{stem}-{side}"]=c
im=Image.new("RGB",(1800,1250),"white");d=ImageDraw.Draw(im);rootlow=np.array([.005,1.711,.016]);roothigh=np.array([.020,1.728,.033]);wholelow=np.array([-.044,1.711,.014]);wholehigh=np.array([.044,1.738,.034])
configs=[(rootlow,roothigh,np.array([1.,0,0]),np.array([0.,1,0]),"Left terminal source ownership: coronal",True),(wholelow,wholehigh,np.array([1.,0,0]),np.array([0.,0,1]),"Both native ownership regions: superior",False),(rootlow,roothigh,np.array([0.,0,1]),np.array([0.,1,0]),"Left terminal source ownership: sagittal",True)]
for panel,(low,high,u,v,title,leftonly)in enumerate(configs):
 ctr=(low+high)/2;axis=np.cross(u,v);scale=min(540/((high-low)@abs(u)),510/((high-low)@abs(v)));off=np.array([panel*600+300,320]);triangles=[]
 def proj(q):return np.c_[(q-ctr)@u,-(q-ctr)@v]*scale+off
 for name,m in meshes.items():
  if leftonly and name.endswith("-r"):continue
  t=m.triangles;keep=np.all(t.max(1)>=low,axis=1)&np.all(t.min(1)<=high,axis=1)
  for tri in t[keep]:
   n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);shade=.65+.35*abs(n@axis)/max(np.linalg.norm(n),1e-30);triangles.append((float(tri.mean(0)@axis),proj(tri),tuple(int(q*shade)for q in colors[name])))
 for _,q,c in sorted(triangles,key=lambda x:x[0]):d.polygon([tuple(x)for x in q],fill=c)
 d.text((panel*600+15,20),title,fill=(0,0,0))
# Actual triangle/plane intersections: no solids invented at openings.
ss=np.load(B/"native-carotid-terminal-cohort-v7b/mca-l-stations.npz");cc=ss["centres"];aa=np.load(B/"aca-l-source-stations.npz")["centres"];n17=cc[18]-cc[16];n17/=np.linalg.norm(n17);na=aa[1]-aa[0];na/=np.linalg.norm(na)
for panel,(point,normal,title)in enumerate([(cc[17],n17,"True cut at shared open ring17"),(aa[0],na,"True cut at unchanged ACA junction"),(cc[0],cc[1]-cc[0],"True cut at unchanged ICA/MCA seam")]):
 normal=normal/np.linalg.norm(normal);u=np.cross(normal,np.array([0.,0,1.]));u/=np.linalg.norm(u);v=np.cross(normal,u);off=np.array([panel*600+300,950]);scale=40000.
 for name,m in meshes.items():
  if name.endswith("-r"):continue
  tri=m.triangles;dist=(tri-point)@normal;keep=(dist.min(1)<=0)&(dist.max(1)>=0)
  for t,dd in zip(tri[keep],dist[keep]):
   q=[]
   for i,j in[(0,1),(1,2),(2,0)]:
    if dd[i]*dd[j]<0:q.append(t[i]+dd[i]/(dd[i]-dd[j])*(t[j]-t[i]))
   if len(q)==2:
    q=np.array(q);xy=np.c_[(q-point)@u,-(q-point)@v]*scale+off
    if np.max(abs(xy-off))<315:d.line([tuple(p)for p in xy],fill=colors[name],width=3)
 d.text((panel*600+15,710),title,fill=(0,0,0))
d.text((15,630),"Gold: terminal cerebral-ICA stem. Red: M1/sphenoid region. Blue: unchanged ICA .015. Green: unchanged ACA.",fill=(0,0,0));d.text((15,655),"Actual decoded bank: every source position, triangle corner, normal word and UV retained. Ring17 is provisional ownership.",fill=(0,0,0));d.text((15,1215),"No geometric fit or junction repair: inherited ICA/MCA seam and ACA-root defects remain unchanged. No cap faces were added.",fill=(0,0,0));p=O/"actual-ownership-and-cuts.png";im.save(p);print(p)
