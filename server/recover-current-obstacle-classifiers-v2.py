#!/usr/bin/env python3
"""Projected numerical closures of literal current obstacle surfaces only."""
import argparse,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import manifold3d as md
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

def directed_cycles(faces):
 edges=defaultdict(list)
 for t in faces:
  for a,b in zip(t,np.roll(t,-1)):edges[tuple(sorted((int(a),int(b))))].append((int(a),int(b)))
 if any(len(uses)>2 for uses in edges.values()):raise ValueError('Nonmanifold source edges')
 boundary=[uses[0]for uses in edges.values()if len(uses)==1];follow=defaultdict(list);incoming=Counter()
 for a,b in boundary:follow[a].append(b);incoming[b]+=1
 if not boundary or any(len(bs)!=1 or incoming[a]!=1 for a,bs in follow.items())or set(follow)!=set(incoming):raise ValueError('Boundary is not disjoint coherent directed cycles')
 loops=[];remaining=set(follow)
 while remaining:
  start=min(remaining);cur=start;loop=[]
  while cur in remaining:remaining.remove(cur);loop.append(cur);cur=follow[cur][0]
  if cur!=start or len(loop)<3:raise ValueError('Invalid directed source cycle')
  loops.append(np.array(loop,int))
 return loops

def cap(v,source_loop):
 loop=source_loop[::-1];points=v[loop];center=points.mean(0);basis=np.linalg.svd(points-center,full_matrices=False)[2][:2].T;xy=(points-center)@basis;nxy=np.roll(xy,-1,axis=0);area=float((xy[:,0]*nxy[:,1]-xy[:,1]*nxy[:,0]).sum()/2)
 if area<0:basis[:,1]*=-1;xy=(points-center)@basis;area=-area
 local=np.asarray(md.triangulate([xy],epsilon=0.,allow_convex=False),dtype=int)
 if not len(local):raise ValueError('Projected source polygon triangulation is empty')
 ff=loop[local];lookup={int(i):j for j,i in enumerate(loop)};restored=[]
 # Restore only original source controls omitted on polygon boundary edges.
 while True:
  edgeuses=defaultdict(list)
  for fi,t in enumerate(ff):
   for a,b in zip(t,np.roll(t,-1)):edgeuses[tuple(sorted([int(a),int(b)]))].append((fi,int(a),int(b)))
  changed=False
  for uses in edgeuses.values():
   if len(uses)!=1:continue
   fi,a,b=uses[0]
   if a not in lookup or b not in lookup:raise ValueError('Unexpected boundary control')
   ia,ib=lookup[a],lookup[b];chain=[];i=ia
   while i!=ib:
    chain.append(int(loop[i]));i=(i+1)%len(loop)
    if len(chain)>len(loop):raise ValueError('Invalid polygon chain')
   chain.append(b)
   if len(chain)==2:continue
   pa,pb=xy[ia],xy[ib];ab=pb-pa;den=float(ab@ab)
   if den==0:raise ValueError('Coincident projected source boundary endpoints')
   idx=np.array([lookup[x]for x in chain]);fraction=((xy[idx]-pa)@ab)/den;distance=np.linalg.norm(xy[idx]-(pa+fraction[:,None]*ab),axis=1)
   if distance.max()>2e-12 or fraction.min()<-1e-10 or fraction.max()>1+1e-10:raise ValueError('Tool omitted a materially noncollinear source boundary chain')
   t=ff[fi];third=int(next(x for x in t if x not in [a,b]));new=np.array([[u,w,third]for u,w in zip(chain[:-1],chain[1:])]);ff=np.r_[ff[:fi],ff[fi+1:],new];restored.append({'omitted_source_controls':chain[1:-1],'max_projected_edge_deviation_m':float(distance.max())});changed=True;break
  if not changed:break
 actual=Counter((int(a),int(b))for t in ff for a,b in zip(t,np.roll(t,-1)));net=Counter({e:n-actual[e[::-1]]for e,n in actual.items()if n>actual[e[::-1]]});expected=Counter((int(a),int(b))for a,b in zip(loop,np.roll(loop,-1)))
 if net!=expected:raise ValueError('Closure does not reproduce exact opposed source boundary chain')
 projected=(v-center)@basis;tri=projected[ff];aa=tri[:,1]-tri[:,0];bb=tri[:,2]-tri[:,0];areas=(aa[:,0]*bb[:,1]-aa[:,1]*bb[:,0])/2
 if np.any(areas<=0):raise ValueError('Projected closure has reversed or zero-area triangles')
 if abs(float(areas.sum())-area)>max(1e-16,area*1e-9):raise ValueError('Projected closure area does not replay source polygon')
 return ff,{'source_controls':len(loop),'closure_triangles':len(ff),'all_source_boundary_controls_used':set(loop).issubset(set(ff.ravel())),'exact_opposed_boundary_chain':True,'projected_area_m2':area,'projected_triangle_area_sum_m2':float(areas.sum()),'source_nonplanarity_preserved_m':float(abs((points-center)@np.cross(basis[:,0],basis[:,1])).max()),'restored_collinear_source_chains':restored}

def main():
 p=argparse.ArgumentParser();p.add_argument('--organ',choices=['stomach','vein'],required=True);p.add_argument('--output-name');args=p.parse_args();B=ROOT/'.local-models/geometry-workbench/circulatory-completion/pancreatic-public-cohort-v1/current-critical-obstacle-recovery-v1';P=B/'projected-closure-plan-v2.json';plan=json.loads(P.read_text());source=Path(plan['sources'][args.organ]['path']);assert digest_file(source)==plan['sources'][args.organ]['sha256'];z=np.load(source);v=z['vertices'];f=z['faces'];uv,iv=np.unique(v,axis=0,return_inverse=True);f0=iv[f];O=B/(args.output_name or(args.organ+'-projected-closure-v2'));O.mkdir(exist_ok=False);proof={'schema':'literal-source-numerical-obstacle-closure/v2','plan_sha256':digest_file(P),'source_sha256':digest_file(source),'organ':args.organ,'source_positions_moved':False,'anatomical_caps_created':False,'physiological_wall_claimed':False,'passed':False,'closures':[]}
 try:
  loops=directed_cycles(f0);closures=[]
  for loop in loops:
   ff,e=cap(uv,loop);closures.append(ff);proof['closures'].append(e);atomic_json(O/'checkpoint.json',proof,overwrite=(O/'checkpoint.json').exists())
  ff=np.r_[f0,*closures];m=Mesh(uv,ff,str(z['frame']));q=audit(m,require_closed=True,expected_components=1,expected_genus=[0]);literal=Mesh(uv,f0,m.frame).geometry_sha==Mesh(v,f,m.frame).geometry_sha;proof.update(gate=q,source_physical_facets_literal=literal,source_physical_faces=len(f0),computational_faces=len(ff)-len(f0),passed=bool(q['passed']and literal and q['topology']['oriented_triangle_integral_m3']>0));normal=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-300);atomic_npz(O/'closed-classifier.npz',vertices=uv,faces=ff,frame=z['frame'],source_physical_faces=f0,computational_cap_faces=np.concatenate(closures),geometric_face_unit_normals=normal)
 except Exception as e:proof.update(error_type=type(e).__name__,error=str(e))
 atomic_json(O/'proof.json',proof);print(args.organ,'PASS',proof['passed'],proof.get('error'),proof.get('gate',{}).get('failures'),flush=True)
if __name__=='__main__':main()
