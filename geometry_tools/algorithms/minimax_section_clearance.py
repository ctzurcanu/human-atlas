"""Bound the worst clearance deficit with one explicit epigraph scalar."""
import numpy as np
from scipy import sparse
from .corotated_soft_section_clearance import solve_soft_sections

def solve_minimax_sections(stiffness,hard,hard_rhs,soft,soft_rhs,*,clearance_rows,clearance_rhs,lexicographic=False,**kwargs):
    S=sparse.csc_matrix(stiffness);A=sparse.csc_matrix(hard);B=sparse.csc_matrix(soft)
    C=sparse.csc_matrix(clearance_rows);h=np.asarray(clearance_rhs,float)
    n=A.shape[1];count=S.shape[1]
    if n!=count*3 or C.shape!=(len(h),n):raise ValueError('Explicit section and clearance dimensions required')
    # t is the first coordinate of an auxiliary vector. Its two unused
    # coordinates are minimized to zero by the existing positive ridge.
    t=sparse.csc_matrix((np.ones(len(h)),(np.arange(len(h)),np.zeros(len(h),int))),shape=(len(h),3))
    lower=sparse.csc_matrix(([1.],([0],[n])),shape=(1,n+3))
    augmented_hard=sparse.vstack([sparse.hstack([A,sparse.csc_matrix((A.shape[0],3))]),
        sparse.hstack([C,t]),lower],format='csc')
    augmented_soft=sparse.vstack([sparse.hstack([B,sparse.csc_matrix((B.shape[0],3))]),-lower],format='csc')
    first_proof=None
    if lexicographic:
        phase=dict(kwargs);limit=phase.pop('time_limit_seconds',5.)
        phase['time_limit_seconds']=limit*.25
        # This objective is measured in millimetres. A 1e-10 objective gap
        # is 1e-13 metres; it is separate from the unchanged feasibility gate.
        phase['objective_gap_tolerance']=1e-10
        objective=np.zeros((count+1,3));objective[count,0]=1.
        # Remove the flat nullspace without meaningful smoothing tradeoff.
        # Motion balls give an explicit upper bound on objective perturbation.
        radius=kwargs['maximum_motion_m']/.001
        ridge=min(1e-10,2e-9/((count+1)*radius**2*.001))
        peak_bias=.5*ridge*(count+1)*radius**2*.001
        first,first_proof=solve_soft_sections(sparse.csc_matrix((0,count+1)),augmented_hard,
            np.r_[hard_rhs,h,0.],sparse.csc_matrix((1,n+3)),np.zeros(1),
            regularizer_weight=ridge,linear_section_objective=objective,**phase)
        first_proof['maximumRegularizerPeakBias_m']=peak_bias
        if first is None:
            return None,dict(complete=True,feasible=False,lexicographicClearance=True,
                reason='First-priority native minimax stage did not return a verified feasible point',firstPriority=first_proof)
        optimum=max(float(first[count,0]),0.)
        # Second priority cannot exchange more than 10 nm of clearance for
        # a better smoothing objective. All original hard rows remain literal.
        augmented_hard=sparse.vstack([augmented_hard,-lower],format='csc')
        hard_rhs=np.r_[hard_rhs,h,0.,-(optimum+1e-8)]
        kwargs['time_limit_seconds']=limit*.75
    else:hard_rhs=np.r_[hard_rhs,h,0.]
    value,proof=solve_soft_sections(sparse.hstack([S,sparse.csc_matrix((S.shape[0],1))]),
        augmented_hard,hard_rhs,augmented_soft,np.r_[soft_rhs,0.],**kwargs)
    proof=dict(proof,clearanceObjective='Squared worst deficit plus all exposed-control deficits and original movement/smoothness',
        epigraphRows=len(h),oneStepPercentageQuota=False)
    if lexicographic:
        proof.update(lexicographicClearance=True,firstPriority=first_proof,
            firstPriorityWorstDeficit_m=optimum,maximumSecondPriorityClearanceTradeoff_m=1e-8,
            firstPriorityGlobalOptimalityClaimed=False,firstPriorityAttainedPrimalClearanceIndependentlyChecked=True)
    if value is None:return None,proof
    d=value[:count];peak=float(value[count,0])
    residual=float(np.maximum(h-C@d.ravel()-peak,0).max(initial=0))
    proof.update(worstClearanceDeficit_m=peak,maximumEpigraphResidual_m=residual)
    if peak < -1e-10 or residual>1e-10:return None,dict(proof,feasible=False)
    return d,proof
