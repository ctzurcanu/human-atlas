"""Exact finite-pose two-member packing with an explicit native deadline.

Inputs have already cleared the fixed context and any preceding member.
The objective minimizes the greater of two caller-supplied nonnegative costs.
A finite-pose optimum is not continuous or three-dimensional acceptance.
"""
import time
import numpy as np
import shapely


def minimum_cost_clear_pair(first,second,first_cost,second_cost,*,margin_m=.0001,wall_seconds=2.):
 first,second=np.asarray(first,object),np.asarray(second,object)
 a,b=np.asarray(first_cost,float),np.asarray(second_cost,float)
 if first.ndim!=1 or second.ndim!=1 or a.shape!=first.shape or b.shape!=second.shape or not np.isfinite(np.r_[a,b,margin_m,wall_seconds]).all()or min(margin_m,wall_seconds)<=0 or np.any(a<0)or np.any(b<0):
  raise ValueError('Finite pose shapes, nonnegative costs, positive margin and deadline required')
 if not np.all(shapely.is_valid(first))or not np.all(shapely.is_valid(second))or np.any(shapely.area(first)<=0)or np.any(shapely.area(second)<=0):raise ValueError('Qualified positive-area polygons required')
 start=time.monotonic();deadline=start+wall_seconds;best=None;bound=np.inf;examined=0;complete=True
 first_order=np.argsort(a,kind='stable');second_order=np.argsort(b,kind='stable')
 for i in first_order:
  if a[i]>=bound:break
  if time.monotonic()>deadline:complete=False;break
  candidates=second_order[b[second_order]<bound]
  if not len(candidates):break
  clear=shapely.disjoint(first[i],second[candidates])&(shapely.distance(first[i],second[candidates])>=margin_m);examined+=len(candidates)
  possible=candidates[clear]
  if len(possible):
   j=int(possible[0]);score=float(max(a[i],b[j]))
   if score<bound:bound=score;best=(int(i),j)
 # Even the last native vector predicate may consume the remaining budget.
 if time.monotonic()>deadline:complete=False
 return dict(feasible=best is not None,complete=complete,indices=None if best is None else list(best),minimumMaximumCost=None if best is None else bound,examinedPairs=examined,elapsedSeconds=time.monotonic()-start)
