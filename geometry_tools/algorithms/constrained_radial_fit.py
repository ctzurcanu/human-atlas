"""Finite two-stage LP for one symmetric, positive-Jacobian radial map.

Only a scalar C2 cubic spline is constructed.  This module does not move an
anatomical mesh, infer tissue identity, fill missing rays or accept anatomy.
All optimization lengths use millimetres to avoid poorly scaled metre LPs.
"""
import time
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import linprog


class RadialFitInfeasible(ValueError):
    """The original supplied host/shape bounds have no feasible solution."""


def _solve(objective, matrix, upper, bounds, deadline, stage):
    left=deadline-time.monotonic()
    if left<=0:
        raise TimeoutError('Common radial fit wall deadline before '+stage)
    started=time.monotonic()
    result=linprog(objective,A_ub=matrix,b_ub=upper,bounds=bounds,method='highs',
        options=dict(time_limit=min(1.,left),presolve=False,primal_feasibility_tolerance=1e-9,
                     dual_feasibility_tolerance=1e-9))
    evidence=dict(stage=stage,status=int(result.status),message=result.message,
        iterations=int(result.nit),elapsedSeconds=time.monotonic()-started,
        nativeTimeLimitSeconds=min(1.,left),presolve=False)
    if result.status==2:
        raise RadialFitInfeasible(stage+': '+result.message)
    if result.status==1 or time.monotonic()>=deadline:
        raise TimeoutError(stage+': native LP or common wall deadline')
    if not result.success:
        raise ValueError(stage+': '+result.message)
    return result,evidence


def fit_radial_shell(radii_m,lower_m,upper_m,*,pad_m=.0005,
        knot_spacing_m=.0005,derivative_bounds=(.75,1.5),
        ratio_bounds=(.8,1.05),lexicographic_slack_m=1e-6,
        wall_seconds=3.,max_knots=64,max_observations=100000):
    """Fit g(r) with explicit host bounds and analytic spline coefficient gates.

    Stage one minimizes max coefficient displacement from identity.  Stage
    two keeps that optimum plus <=one micrometre and minimizes L1 variation
    of the quadratic derivative coefficients.  No candidate or bound retry.
    """
    started=time.monotonic();deadline=started+wall_seconds
    r,low,high=[np.asarray(x,float)for x in [radii_m,lower_m,upper_m]]
    if (r.ndim!=1 or low.shape!=r.shape or high.shape!=r.shape or not len(r)
            or not np.isfinite(np.r_[r,low,high]).all() or np.any(r<=0)
            or np.any(low>high)):
        raise ValueError('Finite corresponding positive radii and ordered host bounds required')
    if not(0<wall_seconds<=3 and 0<pad_m<r.min() and knot_spacing_m>0
           and 0<=lexicographic_slack_m<=1e-6 and 8<=max_knots<=64
           and 1<=max_observations<=100000):
        raise ValueError('Explicit bounded LP resources and positive spline domain required')
    if len(r)>max_observations:
        raise ValueError('Finite radial observation budget exceeded')
    dmin,dmax=map(float,derivative_bounds);smin,smax=map(float,ratio_bounds)
    if not(0<dmin<=dmax and 0<smin<=smax and np.isfinite([dmin,dmax,smin,smax]).all()):
        raise ValueError('Positive finite radial and tangential eigenvalue bounds required')
    # Cubic clamping; every interior knot is simple, hence C2 throughout.
    domain=np.array([r.min()-pad_m,r.max()+pad_m])*1000
    intervals=max(1,int(np.ceil(np.ptp(domain)/(knot_spacing_m*1000))))
    grid=np.linspace(*domain,intervals+1)
    knots=np.r_[np.repeat(grid[0],4),grid[1:-1],np.repeat(grid[-1],4)]
    if len(knots)>max_knots:
        raise ValueError('Finite cubic radial knot budget exceeded')
    n=len(knots)-4;greville=np.array([knots[i+1:i+4].mean()for i in range(n)])
    basis=BSpline.design_matrix(r*1000,knots,3,extrapolate=False).toarray()
    derivative=np.zeros((n-1,n))
    factor=3/(knots[np.arange(n-1)+4]-knots[np.arange(n-1)+1])
    derivative[np.arange(n-1),np.arange(n-1)]=-factor
    derivative[np.arange(n-1),np.arange(n-1)+1]=factor
    matrix=np.concatenate([basis,-basis,derivative,-derivative])
    upper=np.r_[high*1000,-low*1000,np.full(n-1,dmax),np.full(n-1,-dmin)]
    ratio_limits=np.c_[smin*greville,smax*greville]
    # |c_i - greville_i| <= t, coefficient-wise L-infinity identity distance.
    first_matrix=np.vstack([np.c_[matrix,np.zeros(len(matrix))],
        np.c_[np.eye(n),-np.ones(n)],np.c_[-np.eye(n),-np.ones(n)]])
    first_upper=np.r_[upper,greville,-greville]
    first,e1=_solve(np.r_[np.zeros(n),1.],first_matrix,first_upper,
                    [tuple(x)for x in ratio_limits]+[(0,None)],deadline,'minimum-identity-displacement')
    optimum=float(first.x[-1]);limit=optimum+lexicographic_slack_m*1000
    second_limits=np.c_[np.maximum(ratio_limits[:,0],greville-limit),
                        np.minimum(ratio_limits[:,1],greville+limit)]
    variation=np.diff(derivative,axis=0);nv=len(variation)
    second_matrix=np.vstack([np.c_[matrix,np.zeros((len(matrix),nv))],
        np.c_[variation,-np.eye(nv)],np.c_[-variation,-np.eye(nv)]])
    second_upper=np.r_[upper,np.zeros(nv*2)]
    second,e2=_solve(np.r_[np.zeros(n),np.ones(nv)],second_matrix,second_upper,
        [tuple(x)for x in second_limits]+[(0,None)]*nv,deadline,'minimum-derivative-variation')
    coefficients=second.x[:n];values=basis@coefficients;dc=derivative@coefficients
    residual=max(float(np.maximum(low*1000-values,0).max()),
                 float(np.maximum(values-high*1000,0).max()))
    derivative_violation=max(float(np.maximum(dmin-dc,0).max()),float(np.maximum(dc-dmax,0).max()))
    ratios=coefficients/greville
    ratio_violation=max(float(np.maximum(smin-ratios,0).max()),float(np.maximum(ratios-smax,0).max()))
    if residual>2e-8 or max(derivative_violation,ratio_violation)>2e-8:
        raise ValueError('LP result exceeds its independently replayed numerical constraint residual bound')
    assert time.monotonic()<deadline,'Three-second common radial fit deadline'
    return dict(complete=True,feasible=True,degree=3,knots_mm=knots.tolist(),coefficients_mm=coefficients.tolist(),
        greville_mm=greville.tolist(),derivativeCoefficients=dc.tolist(),domain_m=(domain/1000).tolist(),
        observations=len(r),coefficientCount=n,knotCount=len(knots),derivativeBounds=[dmin,dmax],ratioBounds=[smin,smax],
        identityCoefficientDisplacementOptimum_m=optimum/1000,
        lexicographicSlack_m=lexicographic_slack_m,
        actualMaximumIdentityCoefficientDisplacement_m=float(abs(coefficients-greville).max()/1000),
        derivativeCoefficientL1Variation=float(abs(np.diff(dc)).sum()),
        maximumHostConstraintViolation_m=residual/1000,
        maximumDerivativeCoefficientViolation=derivative_violation,maximumRatioCoefficientViolation=ratio_violation,
        analyticProof=dict(basisNonnegativeAndPartitionOfUnity=True,grevilleLinearPrecision=True,
            simpleInteriorKnotsGiveC2=True,derivativeQuadraticBasisNonnegative=True,
            tangentialEigenvalueLower=float(ratios.min()),tangentialEigenvalueUpper=float(ratios.max()),
            radialEigenvalueLower=float(dc.min()),radialEigenvalueUpper=float(dc.max()),
            jacobianDeterminantLower=float(ratios.min()**2*dc.min()),
            scope='Only the stated positive-radius spline interval; no extrapolation, blending or anatomical geometry acceptance'),
        nativeStages=[e1,e2],wallSeconds=wall_seconds,elapsedSeconds=time.monotonic()-started,
        geometryModified=False,anatomicalAcceptance=False)


def evaluate_radial_fit(fit,radii_m,*,derivative=0):
    """Evaluate the saved scalar map/derivative without extrapolating its domain."""
    if not fit.get('complete')or not fit.get('feasible')or derivative not in [0,1,2]:
        raise ValueError('A complete radial fit and derivative order 0, 1 or 2 required')
    r=np.asarray(radii_m,float);domain=np.asarray(fit['domain_m'])
    if not np.isfinite(r).all()or np.any(r<domain[0])or np.any(r>domain[1]):
        raise ValueError('Radial evaluation outside the certified spline interval')
    spline=BSpline(fit['knots_mm'],fit['coefficients_mm'],fit['degree'],extrapolate=False)
    # g has millimetre output and millimetre input; convert only dimensions.
    return spline(r*1000,nu=derivative)*(1000**derivative/1000)
