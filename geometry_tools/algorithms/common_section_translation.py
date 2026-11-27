"""Finite shared-translation tests retaining an entire source section bundle.

Every contour receives the same translation, preserving their relative anatomy.
This tests the declared samples only; no continuous-domain or 3D fit is proved.
"""
import time
import numpy as np
import shapely


def common_section_translations(polygons, translations, context, *,
                                clearance_m=.0001, maximum_motion_m=.025,
                                wall_seconds=5., preserve_source_mutual_gaps=False):
    started=time.monotonic();offsets=np.asarray(translations,float)
    if offsets.ndim!=2 or offsets.shape[1]!=2 or len(offsets)>50000 or not np.isfinite(offsets).all():
        raise ValueError('At most50000 finite explicit planar translations required')
    if not 0<wall_seconds<=10 or clearance_m<0 or maximum_motion_m<=0 or 'container'not in context:
        raise ValueError('Explicit enclosure and finite physical/job limits required')
    clear=np.linalg.norm(offsets,axis=1)<=maximum_motion_m
    source=[shapely.polygons(np.asarray(p,float))for p in polygons]
    if any(not shapely.is_valid(p)for p in source):raise ValueError('Valid complete source contours required')
    # Rigid common translation cannot fix a source mutual collision.
    mutual_margin=0. if preserve_source_mutual_gaps else clearance_m
    source_gaps=[float(shapely.distance(a,b))for i,a in enumerate(source)for b in source[i+1:]]
    mutual=all(shapely.disjoint(a,b)and shapely.distance(a,b)>=mutual_margin
               for i,a in enumerate(source)for b in source[i+1:])
    clear&=mutual
    for points in polygons:
        ids=np.flatnonzero(clear)
        if not len(ids):break
        placed=shapely.polygons(np.asarray(points)[None,:,:]+offsets[ids,None,:])
        fits=shapely.contains(context['container'],placed)
        fits&=shapely.distance(shapely.boundary(context['container']),placed)>=clearance_m
        for key,shape in context.items():
            if key=='container':continue
            fits&=shapely.disjoint(placed,shape)
            fits&=shapely.distance(placed,shape)>=clearance_m
        clear[ids]=fits
        if time.monotonic()-started>wall_seconds:raise TimeoutError('native shared-translation test deadline')
    return np.flatnonzero(clear),dict(samplesTested=len(offsets),validSamples=int(clear.sum()),
        sourceMutualClearance=bool(mutual),relativeSourceGeometryLiteral=True,
        sourceMutualGaps_m=source_gaps,sourceMutualGapsPreservedRatherThanEnlarged=preserve_source_mutual_gaps,
        finiteSamplesOnly=True,whole3DFitNotProven=True,elapsedSeconds=time.monotonic()-started)
