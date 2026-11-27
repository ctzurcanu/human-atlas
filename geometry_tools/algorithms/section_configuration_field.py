"""Bounded clearance proposals using literal non-circular section outlines."""
import time
import numpy as np


class SectionConfigurationField:
    def __init__(self, offsets, clearance, *, maximum_surface_queries=8000000, wall_seconds=16.,
                 ball_clearance=None):
        self.offsets=np.asarray(offsets,float)
        if self.offsets.ndim!=3 or self.offsets.shape[2]!=3 or self.offsets.shape[1]<3 or not np.isfinite(self.offsets).all():
            raise ValueError('Finite actual section outlines at explicit guide stations required')
        if not 0<wall_seconds<=20 or not 1<=maximum_surface_queries<=12000000:
            raise ValueError('Explicit finite section query budgets required')
        self.clearance=clearance;self.maximum=maximum_surface_queries;self.wall=wall_seconds
        self.ball_clearance=ball_clearance
        self.started=time.monotonic();self.queries=0;self.boundingBallAccepted=0

    def __call__(self, points, stations):
        p=np.asarray(points,float);s=np.asarray(stations,float)
        if p.shape!=(len(s),3) or not np.isfinite(p).all() or not np.isfinite(s).all() or np.any(s<0)or np.any(s>len(self.offsets)-1):
            raise ValueError('Finite configuration positions and source station indices required')
        # Nearest original template stays literal. Linear interpolation of
        # offsets can silently narrow an outline when neighbouring frames turn.
        ids=np.rint(s).astype(int);ring=self.offsets[ids]
        result=np.empty(len(p));active=np.arange(len(p))
        if self.ball_clearance is not None:
            # A caller-supplied signed distance lower bound can certify an
            # entire containing ball. This only skips already-safe outline
            # queries; tight configurations still test every literal point.
            radius=np.linalg.norm(ring,axis=2).max(1)
            if self.queries+len(p)>self.maximum:raise ValueError('Finite section surface-query budget reached')
            bound=np.asarray(self.ball_clearance(p),float)-radius
            if bound.shape!=(len(p),):raise ValueError('One bounding-ball lower bound per state required')
            self.queries+=len(p);safe=np.isfinite(bound)&(bound>=0)
            result[safe]=bound[safe];self.boundingBallAccepted+=int(safe.sum())
            active=np.flatnonzero(~safe)
        p=p[active];ring=ring[active]
        # The centre is a declared course/lumen sample. Clear walls cannot
        # justify routing a vascular lumen around an intervening bone sliver.
        samples=np.concatenate([np.zeros((len(ring),1,3)),ring,(ring+np.roll(ring,-1,axis=1))*.5],axis=1)
        count=samples.shape[1];needed=len(p)*count
        if self.queries+needed>self.maximum:raise ValueError('Finite section surface-query budget reached')
        margins=[]
        for first in range(0,len(p),1024):
            if time.monotonic()-self.started>self.wall:raise TimeoutError('Native section configuration deadline')
            q=(p[first:first+1024,None]+samples[first:first+1024]).reshape(-1,3)
            gap=np.asarray(self.clearance(q),float)
            if gap.shape!=(len(q),):raise ValueError('One actual context margin per surface query required')
            margins.append(gap.reshape(-1,count).min(1));self.queries+=len(q)
        result[active]=np.concatenate(margins)if margins else np.empty(0)
        return result
