#!/usr/bin/env python3
"""Materialized source of the passed inline surface-owner diagnostic.

No geometry is emitted on import. This file preserves that diagnostic's
algorithm; the earlier receipt predates this source-file materialization.
"""
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file

def run(base):
 B=Path(base);plan=B/'surface-owner-diagnostic-plan-v1.json';p=json.loads(plan.read_text());source=B/p['source'];z=np.load(source);v=z['vertices'];f=z['faces'];tags=z['source_original_ids'];sourcefaces=z['source_face_ids'];idmap=json.loads(str(z['source_id_map_json']));assert idmap=={'1':'nativeouter','3':'nativefluid','5':'joinedDPAIPAouter','7':'joinedDPAIPAinner','9':'greaterouter','11':'greaterinner','13':'caudalouter','15':'caudalinner'}
 portplan=B/'dorsal-inferior-port-plan-v1.json';pp=json.loads(portplan.read_text());a=np.array(pp['source_witness']['first']);b=np.array(pp['source_witness']['second']);normal=b-a;normal/=np.linalg.norm(normal);point=(a+b)/2;offset=float(point@normal);bank=SharedCutBank(v);dist=(v-point)@normal;faces=[];orig=[];owners=[];sheet=[];newtags=[]
 for fi,t in enumerate(f):
  tag=int(tags[fi]);which='outer'if tag in [1,5,9,13]else'inner';owner='parent'if tag in[1,3]else('greater'if tag in[9,11]else'caudal')
  if tag not in [5,7]:faces.append(t.tolist());orig.append(fi);owners.append(owner);sheet.append(which);newtags.append(tag);continue
  for side,own in [(-1,'dorsal'),(1,'inferior')]:
   poly=[]
   for aa,bb in zip(t,np.roll(t,-1)):
    da,db=dist[aa],dist[bb]
    if side*da>=-1e-12:poly.append(int(aa))
    if da*db<0 and abs(da)>1e-12 and abs(db)>1e-12:poly.append(bank.intersection(aa,bb,normal,offset,'DPA-IPA-declared-midplane'))
   clean=[]
   for j in poly:
    if not clean or clean[-1]!=j:clean.append(j)
   if len(clean)>1 and clean[0]==clean[-1]:clean.pop()
   if clean:
    start=clean.index(min(clean));clean=clean[start:]+clean[:start]
   for j in range(1,len(clean)-1):
    tt=[clean[0],clean[j],clean[j+1]]
    if len(set(tt))==3:faces.append(tt);orig.append(fi);owners.append(own);sheet.append(which);newtags.append(tag)
 v2=np.array(bank.positions);f2=np.array(faces);orig=np.array(orig);owners=np.array(owners);sheet=np.array(sheet);newtags=np.array(newtags);m=Mesh(v2,f2,str(z['frame']));O=B/'surface-owner-diagnostic-v1';O.mkdir(exist_ok=False);t0=v[f];ar0=np.linalg.norm(np.cross(t0[:,1]-t0[:,0],t0[:,2]-t0[:,0]),axis=1)/2;t1=v2[f2];ar1=np.linalg.norm(np.cross(t1[:,1]-t1[:,0],t1[:,2]-t1[:,0]),axis=1)/2;err=np.bincount(orig,weights=ar1,minlength=len(f))-ar0;unselected=np.flatnonzero(~np.isin(tags,[5,7]));literal=all(np.array_equal(v2[f2[np.flatnonzero(orig==i)[0]]],v[f[i]])and np.sum(orig==i)==1 for i in unselected);n=np.cross(t0[:,1]-t0[:,0],t0[:,2]-t0[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];planeerror=abs(np.einsum('fki,fi->fk',t1-t0[orig,0,None,:],n[orig])).max();edges=defaultdict(list)
 for fi,t in enumerate(f2):
  for aa,bb in zip(t,np.roll(t,-1)):edges[tuple(sorted([int(aa),int(bb)]))].append((fi,int(aa),int(bb)))
 seams=defaultdict(list);amb=[]
 for edge,uses in edges.items():
  if len(uses)!=2:continue
  first,second=uses;oa,ob=owners[first[0]],owners[second[0]]
  if oa==ob:continue
  pair=tuple(sorted([str(oa),str(ob)]));sa,sb=sheet[first[0]],sheet[second[0]]
  if sa!=sb:amb.append({'edge':edge,'owners':pair,'sheets':[str(sa),str(sb)],'faces':[first[0],second[0]]});continue
  use=first if oa==pair[0]else second;seams[(pair,str(sa))].append((use[1],use[2]))
 loops={};loop_fail=[]
 for (pair,which),ee in seams.items():
  follow=defaultdict(list);incoming=Counter()
  for aa,bb in ee:follow[aa].append(bb);incoming[bb]+=1
  cycles=[]
  if any(len(bs)!=1 or incoming[aa]!=1 for aa,bs in follow.items())or set(follow)!=set(incoming):loop_fail.append({'pair':pair,'sheet':which,'reason':'branched or open directed boundary','edges':ee});continue
  remaining=set(follow)
  while remaining:
   start=min(remaining);cur=start;cycle=[]
   while cur in remaining:remaining.remove(cur);cycle.append(cur);cur=follow[cur][0]
   if cur!=start:loop_fail.append({'pair':pair,'sheet':which,'reason':'does not return to start'})
   cycles.append(cycle)
  loops['|'.join(pair)+'|'+which]={'cycles':cycles,'coordinates_metres':[v2[c]for c in cycles],'edge_count':len(ee)}
 expected={tuple(sorted(pair))for pair in p['expected_pairs']};actual={pair for pair,which in seams};loopcounts={key:len(val['cycles'])for key,val in loops.items()};gate=audit(m,require_closed=True,expected_components=1,expected_genus=[3]);okay=gate['passed']and literal and np.array_equal(v2[:len(v)],v)and np.max(abs(err))<1e-14 and planeerror<1e-13 and not amb and not loop_fail and actual==expected and all('|'.join(pair)+'|'+s in loops and len(loops['|'.join(pair)+'|'+s]['cycles'])==1 for pair in expected for s in ['outer','inner'])
 atomic_npz(O/'conforming-whole-boundary.npz',vertices=v2,faces=f2,frame=z['frame'],whole_source_face_ancestry=orig,source_original_ids=newtags,source_physical_face_ids=sourcefaces[orig],owners=owners,sheets=sheet,source_original_vertices=v,source_original_faces=f,midplane_point=point,midplane_normal=normal)
 result={'schema':'pancreatic-conforming-surface-owner-diagnostic/v1','passed':okay,'plan_sha256':digest_file(plan),'portplan_sha256':digest_file(portplan),'source_wall_sha256':digest_file(source),'source_id_map':idmap,'original_prefix_literal':bool(np.array_equal(v2[:len(v)],v)),'all_non_DPAIPA_source_facets_literal':literal,'max_source_face_area_replay_error_m2':float(np.max(abs(err))),'max_new_triangle_distance_from_original_parent_plane_m':float(planeerror),'additional_edge_intersection_controls':len(v2)-len(v),'source_triangles':len(f),'conforming_triangles':len(f2),'owner_face_counts':dict(Counter(owners.tolist())),'gate':gate,'expected_owner_pairs':sorted(expected),'actual_owner_pairs':sorted(actual),'directed_seam_loops':loops,'loop_counts':loopcounts,'ambiguous_sheet_edges':amb,'loop_failures':loop_fail,'internal_caps_created':False,'anatomical_acceptance':False,'installed':False};atomic_json(O/'proof.json',result);return result

if __name__=='__main__':
 run(ROOT/'.local-models/geometry-workbench/circulatory-completion/pancreatic-public-cohort-v1/finite-network-construction-v1')
