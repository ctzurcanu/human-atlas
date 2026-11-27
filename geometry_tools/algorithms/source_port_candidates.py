"""Enumerate source port candidates without making motion a dependency gate.

Candidate proximity is only a search filter. The caller must compare complete
source contours and record unresolved or unequal mouths explicitly.
"""
import numpy as np
from scipy.spatial import cKDTree

def source_port_candidates(ports, *, include_unmoved=True, radius_m=.005,
                           include_same_part=False):
 if not np.isfinite(radius_m)or radius_m<=0:raise ValueError('Positive source search radius required')
 if not ports:return
 centres=np.array([np.asarray(p['source'],float).mean(0)for p in ports])
 if centres.shape!=(len(ports),3)or not np.isfinite(centres).all():raise ValueError('Finite source port polygons required')
 for i,j in sorted(cKDTree(centres).query_pairs(radius_m)):
  a,b=ports[i],ports[j]
  if a['system']!=b['system']:continue
  if a['partId']==b['partId']:
   if not include_same_part:continue
   if 'component'not in a or 'component'not in b:raise ValueError('Same-part port search requires explicit topological component ownership')
   if a['component']==b['component']:continue
  if not include_unmoved and not(a['moving']or b['moving']):continue
  yield i,j
