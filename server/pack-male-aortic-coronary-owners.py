"""Emit three disjoint owners with literal native far corners and shared annuli."""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh
from geometry_tools.common import ROOT,atomic_json
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim
from geometry_tools.tissue import pack,normal_words

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-coronary-ostia-v2';O=D/'owned-bank-v1';z=np.load(D/'surfaces-and-ancestry.npz');v=z['vertices'];outer=z['outer_faces'];inner=z['inner_faces'];outer_component=z['outer_source_component'];inner_component=z['inner_source_component'];outer_face=z['outer_source_face'];inner_face=z['inner_source_face'];native_context=ROOT/'.local-models/geometry-workbench/heart/male-native-aortic-source-defects-v2';raw=[Mesh.load(ROOT/'.local-models/geometry-workbench/heart/male-aortic-section-transport-v1/candidate-outer.npz'),Mesh.load(native_context/'Left-coronary-artery.npz'),Mesh.load(native_context/'Right-coronary-artery.npz')];rightbranch=max(recover_station_rings(raw[2].vertices,raw[2].faces),key=lambda b:max(raw[2].vertices[b['rings'][0],1].mean(),raw[2].vertices[b['rings'][-1],1].mean()));root_ids=np.concatenate(rightbranch['rings'][:13]);right_mask=np.all(np.isin(raw[2].faces,root_ids),axis=1);raw[2]=Mesh(raw[2].vertices,raw[2].faces[right_mask],raw[2].frame,raw[2].attributes);sources=[];stations=[];rings=[]
for i,m in enumerate(raw):
 used,iv=np.unique(m.faces,return_inverse=True);source=Mesh(m.vertices[used],iv.reshape(-1,3),m.frame,{k:a[used]for k,a in m.attributes.items()});b=recover_station_rings(source.vertices,source.faces)[0];r=list(b['rings']);
 if i in [1,2]and source.vertices[r[0],1].mean()<source.vertices[r[-1],1].mean():r=r[::-1]
 st=np.empty(len(source.vertices),int)
 for j,rr in enumerate(r):st[rr]=j
 sources.append(source);stations.append(st);rings.append(r)
thresholds={1:14,2:3};names={0:'aorta',1:'left-coronary',2:'right-coronary'};surfaces={};cells={n:[]for n in names.values()};corner_uvs={};corner_words={};interpolation=[]
def ownership(component,face):
 result=np.zeros(len(component),int)
 for i,cut in thresholds.items():
  ids=np.flatnonzero((component==i)&(face>=0)&(face<len(sources[i].faces)));kept=ids[stations[i][sources[i].faces[face[ids]]].min(1)>=cut];result[kept]=i
 return result
outer_owner=ownership(outer_component,outer_face);inner_owner=ownership(inner_component,inner_face)
for i,name in names.items():
 for sheet,faces,owners in [('outer',outer,outer_owner),('inner',inner,inner_owner)]:
  key=name+':'+sheet;surfaces[key]=faces[owners==i];cells[name].append(dict(surface=key,side=1,neighbor=None))
for port in json.loads((D/'construction-proof.json').read_text())['ports']:
 i=port['component'];key=f'external:{i}:{port["station"]}';surfaces[key]=z[f'annulus_{[(0,0),(0,24),(1,24),(2,12)].index((i,port["station"]))}_faces'];cells[names[i]].append(dict(surface=key,side=1,neighbor=None))
for i,cut in thresholds.items():
 expected=sources[i].vertices[rings[i][cut]];center=expected.mean(0);outer_tail=surfaces[names[i]+':outer'];inner_tail=surfaces[names[i]+':inner'];ol=min(boundary_loops(outer_tail),key=lambda r:np.linalg.norm(v[r].mean(0)-center));il=min(boundary_loops(inner_tail),key=lambda r:np.linalg.norm(v[r].mean(0)-center));target=v[ol]-(v[ol]-center)/np.linalg.norm(v[ol]-center,axis=1)[:,None]*.0003;order=np.argmin(np.linalg.norm(v[il][None]-target[:,None],axis=2),axis=1);assert len(ol)==len(il)==len(set(order))==12;assert np.linalg.norm(v[il[order]]-target,axis=1).max()<1e-11;key='shared:'+names[i];surfaces[key]=annular_rim(ol,il[order]);cells[names[i]].append(dict(surface=key,side=1,neighbor='aorta'));cells['aorta'].append(dict(surface=key,side=-1,neighbor=names[i]))
words=normal_words(v.astype('<f4').astype(float),surfaces);literal_source_corners=0;interpolated_source_corners=0;corrected_corners=0
for i,name in names.items():
 mask=outer_owner==i;ff=outer[mask];cc=outer_component[mask];anc=outer_face[mask];key=name+':outer';uv=np.zeros((len(ff),3,2),'<f4');normal=words[key].copy()
 for j,(face,component,source_face)in enumerate(zip(ff,cc,anc)):
  source=sources[int(component)]
  if source_face>=len(source.faces):continue
  tri=source.vertices[source.faces[source_face]];fn=np.cross(v[face[1]]-v[face[0]],v[face[2]]-v[face[0]]);e=tri[1:]-tri[0];gram=e@e.T
  for k,point in enumerate(v[face]):
   equal=np.flatnonzero(np.all(tri==point,axis=1))
   if len(equal):weights=np.eye(3)[equal[0]];literal_source_corners+=1
   else:
    xy=np.linalg.solve(gram,e@(point-tri[0]));weights=np.r_[1-xy.sum(),xy];assert np.min(weights)>-1e-8 and np.max(weights)<1+1e-8;interpolated_source_corners+=1
   uv[j,k]=(weights@source.attributes['uvs'][source.faces[source_face]]).astype('<f4');n=weights@source.attributes['normals'][source.faces[source_face]].astype(float)
   if n@fn<=0:corrected_corners+=1;continue
   normal[j,k]=np.rint(n).astype('<i2')if len(equal)else np.rint(n/np.linalg.norm(n)*32767).astype('<i2')
 corner_uvs[key]=uv;words[key]=normal
 for sheet in ['inner']:
  sid=name+':'+sheet;corner_uvs[sid]=np.zeros((len(surfaces[sid]),3,2),'<f4')
for sid,ff in surfaces.items():corner_uvs.setdefault(sid,np.zeros((len(ff),3,2),'<f4'))
ids={'aorta':'TA98ARTERIAL:male:ascending-aorta:finite-wall','left-coronary':'TA98ARTERIAL:male:left-coronary:root-wall','right-coronary':'TA98ARTERIAL:male:right-coronary:root-wall'};result=pack(v,surfaces,cells,ids,raw[0].frame,O,corner_words=words,corner_uvs=corner_uvs,max_pairs=250000,max_tests=2000000);assert result['passed'];report=dict(schema='male-aortic-coronary-owned-bank/v1',source_construction=str(D/'construction-proof.json'),bank_sha256=result['bank_sha256'],bank_bytes=result['total_bytes'],finite_cells=3,reciprocal_owned_annuli=2,all_cells_and_pair_occupied_volumes_passed=result['passed'],native_source_uv_corners_literal_or_barycentric=True,literal_source_outer_corners=literal_source_corners,interpolated_clipped_source_corners=interpolated_source_corners,contradictory_outer_normal_corners_reconciled=corrected_corners,all_far_native_vessel_positions_literal=True,heart_geometry_and_all113controls_unchanged=True,heart_and_arch_joins_pending=True,installable=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T);atomic_json(D/'owned-bank-stage.json',report);print(json.dumps(report))
