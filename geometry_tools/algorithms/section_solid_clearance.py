"""Whole polygon footprint clearance in a literal closed source section.

Unlike testing only sampled vessel vertices this also rejects edge crossings
and a forbidden structure enclosed by a vessel outline. No contour is closed
or repaired here. The caller must qualify anatomical occupancy separately.
"""
import numpy as np
import shapely

def translated_section_clearance(polygons, offsets, obstacles, *, margin_m):
 offsets=np.asarray(offsets,float)
 if offsets.ndim!=2 or offsets.shape[1]!=2 or not np.isfinite(offsets).all()or not np.isfinite(margin_m)or margin_m<0:raise ValueError('Finite planar translations and nonnegative clearance required')
 def polygon(points):
  p=np.asarray(points,float)
  if p.ndim!=2 or p.shape[1]!=2 or len(p)<3 or not np.isfinite(p).all():raise ValueError('Actual finite closed-section polygons required')
  out=shapely.polygons(p)
  if not shapely.is_valid(out)or shapely.area(out)<=0:raise ValueError('Invalid source section; no automatic closure or repair permitted')
  return p,out
 source=[polygon(p)for p in polygons];target=[polygon(p)[1]for p in obstacles]
 if not source:raise ValueError('A complete source bundle is required')
 if not target:return np.ones(len(offsets),bool)
 forbidden=shapely.union_all(target);clear=np.ones(len(offsets),bool)
 for p,_ in source:
  placed=shapely.polygons(p[None]+offsets[:,None]);clear&=shapely.disjoint(placed,forbidden)&(shapely.distance(placed,forbidden)>=margin_m)
 return clear
