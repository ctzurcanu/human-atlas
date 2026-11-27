"""Draw literal triangle/plane intersection lines for packing diagnosis."""
from pathlib import Path
import argparse,json,sys,time,html
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--catalogue');ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
path=ROOT/(args.catalogue or f'.local-models/ta98-review/bowel-wall-local-atlas/{args.sex}/atlas.json');a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres';out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
def read(pid):return reader.read_part(a,parts[pid],frame=frame,catalogue_path=path)
heart=read(f'TA98HEART:{args.sex}:epicardium');hidden=set(a.get('defaultHidden',[]));lo,hi=heart.bounds+np.array([[-.05]*3,[.05]*3])
groups=[('Chest wall','#828b98',[p['id']for p in a['parts']if p['system']in ['skeletal','muscular']and p['id']not in hidden and not p.get('suppressed')and not p.get('defaultHidden')and np.all(np.array(p['bounds'])[1]>=lo)and np.all(np.array(p['bounds'])[0]<=hi)]),('Lungs','#31adc8',a['pulmonarySegmentIntegration']['partIds']),('Airways','#f4c354',a['pulmonaryBranchIntegration']['partIds']),('Heart','#f16d77',[f'TA98HEART:{args.sex}:epicardium'])]
triangles=[(name,color,np.concatenate([read(i).triangles for i in ids]))for name,color,ids in groups]
levels=np.linspace(heart.bounds[0,1]+.015,heart.bounds[1,1]-.015,4);svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="850" viewBox="0 0 1200 850"><rect width="1200" height="850" fill="#131b23"/>']
def text(x,y,s,size=15,color='#dce4eb'):svg.append(f'<text x="{x}" y="{y}" fill="{color}" font-family="sans-serif" font-size="{size}">{html.escape(str(s))}</text>')
text(24,30,args.sex.capitalize()+' — actual axial triangle sections',23)
for n,(name,color,_)in enumerate(triangles):text(25+n*260,57,name,16,color)
for panel,level in enumerate(levels):
 ox=20+(panel%2)*600;oy=90+(panel//2)*380;text(ox+8,oy+15,f'Y = {level*1000:.1f} mm')
 def point(q):return np.c_[ox+280+q[:,0]*1750,oy+205-q[:,2]*1750]
 for name,color,tri in triangles:
  d=tri[:,:,1]-level;keep=(d.min(1)<0)&(d.max(1)>0);t=tri[keep];d=d[keep];segs=[]
  for e,(i,j)in enumerate([(0,1),(1,2),(2,0)]):
   hit=d[:,i]*d[:,j]<0;rows=np.where(hit)[0];v=t[hit,i]+(t[hit,j]-t[hit,i])*(-d[hit,i]/(d[hit,j]-d[hit,i]))[:,None];segs.extend(zip(rows.tolist(),v[:,[0,2]].tolist()))
  endpoints={}
  for row,p in segs:endpoints.setdefault(row,[]).append(p)
  paths=[]
  for ends in endpoints.values():
   if len(ends)!=2:continue
   v=point(np.c_[np.array(ends)[:,0],np.zeros(2),np.array(ends)[:,1]]);paths.append(f'M{v[0,0]:.2f},{v[0,1]:.2f}L{v[1,0]:.2f},{v[1,1]:.2f}')
  svg.append(f'<path d="{" ".join(paths)}" stroke="{color}" stroke-width="{1.2 if name=="Heart" else .7}" fill="none"/>')
 text(ox+8,oy+345,'Right ← x → Left     Anterior ↑',12)
svg.append('</svg>');(out/'sections.svg').write_text(''.join(svg));atomic_json(out/'source.json',dict(catalogue_sha256=digest_file(path),levels_m=levels.tolist(),method='literal triangle/plane crossings; no replacement or smoothing',elapsed_seconds=time.monotonic()-start));print(out/'sections.svg')
