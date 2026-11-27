"""Actual triangle-plane sections of three immutable native donor assets."""
from pathlib import Path
import collections, hashlib, json, sys, time
import numpy as np
from PIL import Image,ImageDraw,ImageFont

ROOT=Path(__file__).resolve().parents[1]

def inside(point,poly):
 x,y=point;px,py=poly.T;qx,qy=np.roll(poly,1,axis=0).T;mask=(py>y)!=(qy>y)
 return bool(np.count_nonzero(x<(qx[mask]-px[mask])*(y-py[mask])/(qy[mask]-py[mask])+px[mask])%2)

def draw_image(name,v,f,sections,path):
 im=Image.new('RGB',(1800,1040),'white');draw=ImageDraw.Draw(im);font=ImageFont.truetype('/System/Library/Fonts/Helvetica.ttc',18);small=ImageFont.truetype('/System/Library/Fonts/Helvetica.ttc',14)
 draw.text((25,12),name+' — actual source triangles and uncapped section curves',fill='black',font=font);draw.text((25,39),'Native source units; physical scale undeclared. No mesh edits, caps or schematic geometry.',fill='black',font=small)
 def panel(number,title):
  x=(number%3)*600+25;y=(number//3)*480+95;draw.text((x,y-22),title,fill='black',font=small);return x,y,550,415
 def project(points,box):
  lo=points.min(0);hi=points.max(0);scale=min(box[2]/max(hi[0]-lo[0],1e-12),box[3]/max(hi[1]-lo[1],1e-12));center=(lo+hi)/2;return np.c_[box[0]+box[2]/2+(points[:,0]-center[0])*scale,box[1]+box[3]/2-(points[:,1]-center[1])*scale]
 for number,(angle,title)in enumerate([(0,'Front: actual triangle surface'),(.7,'Oblique: actual triangle surface')]):
  box=panel(number,title);projected=np.c_[v[:,0]*np.cos(angle)+v[:,2]*np.sin(angle),v[:,1]];depth=v[:,2]*np.cos(angle)-v[:,0]*np.sin(angle);pp=project(projected,box);t=v[f];norm=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);shade=np.abs(norm@np.array([.4,.2,.8]))/np.maximum(np.linalg.norm(norm,axis=1),1e-15);shade=np.clip(shade,0,1)
  for face in np.argsort(depth[f].mean(1)):
   color=(int(135+95*shade[face]),int(45+80*shade[face]),int(45+80*shade[face]));draw.polygon([tuple(p)for p in pp[f[face]]],fill=color)
 for number,s in enumerate(sections,2):
  box=panel(number,f"{s['normal_axis']}={s['plane_coordinate_native']:.6g}; {len(s['loops'])} curves; depths "+','.join(str(l['containment_depth'])for l in s['loops']));pts=np.asarray(s['points_native']);axes=['XYZ'.index(x)for x in s['projection_axes']];pp=project(pts[:,axes],box)
  for segment in s['segments']:draw.line([tuple(p)for p in pp[segment]],fill='#111111',width=2)
  low=pts[:,axes].min(0);high=pts[:,axes].max(0);draw.text((box[0],box[1]+box[3]+5),f"Axes {s['projection_axes']}; native bounds {np.round(low,4).tolist()} to {np.round(high,4).tolist()}",fill='black',font=small)
 im.save(path)

def cut(v,f,axis,value):
 # Avoid vertex-coincident planes explicitly, without moving mesh controls.
 while np.any(v[:,axis]==value):value=np.nextafter(value,np.inf)
 all_edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);edges,inv=np.unique(np.sort(all_edges,axis=1),axis=0,return_inverse=True);inv=inv.reshape(3,-1).T
 d=v[:,axis]-value;cross=(d[edges[:,0]]<0)!=(d[edges[:,1]]<0);ids=np.flatnonzero(cross);lookup=np.full(len(edges),-1,int);lookup[ids]=np.arange(len(ids));e=edges[ids];t=d[e[:,0]]/(d[e[:,0]]-d[e[:,1]]);points=v[e[:,0]]+t[:,None]*(v[e[:,1]]-v[e[:,0]])
 segments=[];faceids=[]
 for face,es in enumerate(inv):
  hits=lookup[es];hits=hits[hits>=0]
  if len(hits)==2:segments.append(hits.tolist());faceids.append(face)
 segments=np.asarray(segments,int);adj=collections.defaultdict(list)
 for i,(x,y)in enumerate(segments):adj[int(x)].append(int(y));adj[int(y)].append(int(x))
 seen=set();loops=[];bad=[];projection=[x for x in range(3)if x!=axis]
 for start in sorted(adj):
  if start in seen:continue
  stack=[start];component=[]
  while stack:
   x=stack.pop()
   if x in seen:continue
   seen.add(x);component.append(x);stack.extend(adj[x])
  if any(len(adj[x])!=2 for x in component):bad.append(component);continue
  order=[start];previous=None;current=start
  while True:
   nxt=next(x for x in adj[current]if x!=previous)
   if nxt==start:break
   order.append(nxt);previous,current=current,nxt
  poly=points[order][:,projection];area=abs(float(np.sum(poly[:,0]*np.roll(poly[:,1],-1)-poly[:,1]*np.roll(poly[:,0],-1))/2));loops.append(dict(point_ids=order,absolute_area_native2=area))
 for i,loop in enumerate(loops):
  witness=points[loop['point_ids'][0],projection];loop['containment_depth']=sum(inside(witness,points[other['point_ids']][:,projection])for j,other in enumerate(loops)if i!=j)
 return dict(normal_axis='XYZ'[axis],plane_coordinate_native=float(value),projection_axes=['XYZ'[x]for x in projection],points_native=points.tolist(),segments=segments.tolist(),segment_source_face_ids=faceids,loops=loops,noncycle_components=bad,even_odd_cross_section_area_native2=sum((-1)**l['containment_depth']*l['absolute_area_native2']for l in loops),method='Exact triangle-edge plane intersections; original triangles retained; no cap geometry or proximity merge')

def main():
 started=time.monotonic();src=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1/source';out=src.parent/'source-topology-v1/major-trunk-sections';out.mkdir(parents=True,exist_ok=True)
 assets=json.loads((src/'assets.json').read_text())['assets'];instances=json.loads((src/'instances.json').read_text())['instances'];positions=np.memmap(src/'positions.f32',dtype='<f4',mode='r');indices=np.memmap(src/'indices.u32',dtype='<u4',mode='r');report=[]
 for name in ['Thoracic aorta','Abdominal aorta','Superior vena cava']:
  instance=next(i for i in instances if i['name']==name);asset=next(a for a in assets if a['geometryId']==instance['geometryId']);pr,ir=asset['positions'],asset['indices'];v=np.asarray(positions[pr['byteOffset']//4:(pr['byteOffset']+pr['byteLength'])//4]).reshape(-1,3).astype(float);f=np.asarray(indices[ir['byteOffset']//4:(ir['byteOffset']+ir['byteLength'])//4]).reshape(-1,3);matrix=np.asarray(instance['worldMatrix']).reshape((4,4),order='F');v=v@matrix[:3,:3].T+matrix[:3,3]
  # Canonical exact weld only lets shared triangle edges share section points.
  v,inv=np.unique(v,axis=0,return_inverse=True);f=inv[f];bounds=np.array([v.min(0),v.max(0)]);sections=[cut(v,f,1,float(bounds[0,1]+fraction*(bounds[1,1]-bounds[0,1])))for fraction in [.3,.5,.7]]+[cut(v,f,0,float(bounds[:,0].mean()))]
  slug=name.lower().replace(' ','-');png=out/(slug+'.png');draw_image(name,v,f,sections,png);record=dict(source_id=instance['id'],geometryId=asset['geometryId'],source_name=name,coordinate_frame='original-donor-world-native-units',source_transform_applied_once=True,source_faces_unchanged=True,sections=sections,png=str(png),png_sha256=hashlib.sha256(png.read_bytes()).hexdigest());(out/(slug+'.json')).write_text(json.dumps(record,separators=(',',':'))+'\n');report.append(dict(name=name,png=str(png),cuts=[dict(axis=s['normal_axis'],position=s['plane_coordinate_native'],closed_curves=len(s['loops']),depths=[l['containment_depth']for l in s['loops']],areas=[l['absolute_area_native2']for l in s['loops']],noncycle_components=len(s['noncycle_components']))for s in sections]))
 summary=dict(assets=report,source_unchanged=True,selfIntersections='UNTESTED',anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started);(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
