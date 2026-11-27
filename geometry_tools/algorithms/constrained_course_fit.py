"""Finite sequential geometric-course fitting with immutable hard constraints.

Every phase enforces the same physical halfspaces, protected cells and motion
domain. Curvature slack is diagnostic only. A result qualifies only at zero
slack with an independent current-chord curvature replay and strict status.
Surface packing, tube walls and anatomical acceptance remain separate checks.
"""
import time
import numpy as np
from .bounded_cell_fairing import fair_rigid_cells


def fit_course(centres,edges,moving_cells,fixed_cells,targets,*,profile_cells,
               halfspaces=(),plane_count=0,curvature_edges=None,
               curvature_limits_per_m=None,reference_centres=None,
               maximum_motion_m=.025,maximum_steps=4,step_seconds=3.,
               wall_seconds=18.,on_step=None,start_strict=False,surface_guards=None):
    start=time.monotonic();c=np.asarray(centres,float)
    if not isinstance(maximum_steps,int)or not 1<=maximum_steps<=8 or not 0<step_seconds<=8 or not 0<wall_seconds<=60:
        raise ValueError('Finite declared phase count, native step and total deadlines required')
    reference=c.copy()if reference_centres is None else np.array(reference_centres,float,copy=True)
    if reference.shape!=c.shape or not np.isfinite(reference).all():raise ValueError('Finite source-bound starting reference required')
    steps=[];best=None;best_index=None;best_metric=(np.inf,np.inf);qualified=False;strict=bool(start_strict);reason='Geometric quality not reached within phase budget'
    for index in range(maximum_steps):
        remaining=wall_seconds-(time.monotonic()-start)
        if remaining<.25:reason='Native total fitting deadline reached';break
        delta,proof=fair_rigid_cells(c,edges,moving_cells,fixed_cells,targets,
            profile_cells=profile_cells,halfspaces=halfspaces,plane_count=plane_count,
            maximum_motion_m=maximum_motion_m,time_limit_seconds=min(step_seconds,remaining),
            curvature_edges=curvature_edges,curvature_reference_centres=reference,
            curvature_limits_per_m=curvature_limits_per_m,elastic_curvature_phase=not strict,surface_guards=surface_guards)
        step=dict(index=index,proof=proof,strictCurvatureStep=strict,elapsedSeconds=time.monotonic()-start)
        qualified=bool(delta is not None and proof.get('feasible')and proof.get('geometricCurvatureQualified'))
        expired=time.monotonic()-start>wall_seconds
        if expired:qualified=False
        step['qualifiedCourse']=qualified;steps.append(step)
        if on_step is not None:on_step(step,delta)
        if delta is None:
            if strict:strict=False;continue
            reason='No independently verified hard-constraint proposal';break
        metric=(proof['maximumGeometricCurvatureViolation_per_m'],proof['curvaturePhaseSlack_per_m'])
        if metric<best_metric:best=np.array(delta,copy=True);best_metric=metric;best_index=index
        if expired:reason='Native total fitting deadline reached';break
        if qualified:best=np.array(delta,copy=True);best_index=index;reason='Hard proposal and actual discrete course curvature qualified';break
        reference=c+delta
        strict=proof['curvaturePhaseSlack_per_m']<=.01
    report=dict(qualifiedCourse=qualified,complete=qualified,reason=reason,
                steps=steps,selectedStep=best_index,maximumSteps=maximum_steps,nativeStepSeconds=step_seconds,
                nativeWallSeconds=wall_seconds,elapsedSeconds=time.monotonic()-start,
                anatomicalAcceptance=False,installed=False,wholeSurfaceAuditRequired=True)
    return best,report
