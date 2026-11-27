"""Explicit finite vessel walls from paired exterior/lumen triangle banks.

Each open rim is bridged by a tissue annulus. No triangle spans the blood-space
disk. Open exterior classifiers, if needed for volume subtraction, belong to
the caller and must never be mistaken for these anatomical wall surfaces.
"""
from __future__ import annotations
import collections,numpy as np

def oriented_boundary_loops(faces):
 f=np.asarray(faces,dtype=np.uint32);counts=collections.Counter(tuple(sorted((int(a),int(b)))) for t in f for a,b in zip(t,np.roll(t,-1)));directed=[(int(a),int(b)) for t in f for a,b in zip(t,np.roll(t,-1)) if counts[tuple(sorted((int(a),int(b))))]==1];successor={a:b for a,b in directed};assert len(successor)==len(directed)
 unseen=set(successor);out=[]
 while unseen:
  start=min(unseen);ring=[];j=start
  while j not in ring:ring.append(j);unseen.remove(j);j=successor[j]
  assert j==start;out.append(np.array(ring,np.uint32))
 return out

def finite_wall(outer,inner,faces):
 outer=np.asarray(outer,np.float32);inner=np.asarray(inner,np.float32);faces=np.asarray(faces,np.uint32);assert outer.shape==inner.shape;count=len(outer);caps=[];loops=oriented_boundary_loops(faces)
 for ring in loops:
  cap=[]
  for a,b in zip(ring,np.roll(ring,-1)):
   a,b=int(a),int(b);cap.extend([[b,a,a+count],[b,a+count,b+count]])
  caps.append(np.array(cap,np.uint32))
 annuli=np.concatenate(caps) if caps else np.empty((0,3),np.uint32);inner_faces=faces[:,::-1]+count
 return {'PairedPositions':np.r_[outer,inner],'OuterFaces':faces,'InnerFaces':inner_faces,'MouthCapFaces':annuli,'FiniteFaces':np.r_[faces,inner_faces,annuli]},loops,caps
