"""Explicit axial fitting field with a fixed anchor curve and positive stretch.

Each transverse coordinate has independent positive/negative-side scales.
The smooth absolute-value integral preserves the anchor exactly and bounds
the transverse derivative between those scales. Longitudinal coordinates
stay literal. Hence the continuous Jacobian is positive for positive scales,
even when the anchor and scales vary along the axis. This does not establish
discrete triangle quality, tube section area or anatomical acceptance.
"""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicHermiteSpline


@dataclass
class AnchoredAxialCage:
    longitudinal_knots: object
    negative_scales: object
    positive_scales: object
    anchor_knots: object
    anchor_values: object
    transition_width_m: float = .002
    longitudinal_axis: int = 1
    minimum_stretch: float = .5

    def __post_init__(self):
        self.y=np.asarray(self.longitudinal_knots,float)
        self.negative=np.asarray(self.negative_scales,float)
        self.positive=np.asarray(self.positive_scales,float)
        self.ay=np.asarray(self.anchor_knots,float)
        self.av=np.asarray(self.anchor_values,float)
        if self.longitudinal_axis not in [0,1,2]:raise ValueError('An explicit coordinate axis is required')
        self.transverse=[i for i in range(3)if i!=self.longitudinal_axis]
        if any(a.ndim!=1 or len(a)<2 or not np.isfinite(a).all()or np.any(np.diff(a)<=0)for a in [self.y,self.ay]):
            raise ValueError('Finite strictly increasing knot arrays are required')
        if self.negative.shape!=(len(self.y),2)or self.positive.shape!=self.negative.shape or self.av.shape!=(len(self.ay),2):
            raise ValueError('Two transverse values are required at every corresponding knot')
        if not all(np.isfinite(a).all()for a in [self.negative,self.positive,self.av]):raise ValueError('Nonfinite cage values')
        if not np.isfinite(self.minimum_stretch)or self.minimum_stretch<=0:raise ValueError('Positive stretch bound required')
        if min(self.negative.min(),self.positive.min())<self.minimum_stretch:raise ValueError('Declared positive stretch bound violated')
        if not np.isfinite(self.transition_width_m)or self.transition_width_m<=0:raise ValueError('Positive finite transition width required')
        def curve(x,values):
            p=PchipInterpolator(x,values,axis=0,extrapolate=False);slopes=p.derivative()(x)
            # Zero endpoint slopes retain the monotone PCHIP bounds and make
            # the constant extension C1. Literal endpoint values are restored
            # separately, avoiding evaluation roundoff in fixed regions.
            slopes[[0,-1]]=0
            return CubicHermiteSpline(x,values,slopes,axis=0,extrapolate=False)
        self.n=curve(self.y,self.negative)
        self.p=curve(self.y,self.positive)
        self.a=curve(self.ay,self.av)

    @staticmethod
    def _values(interpolator,y,ordinates):
        clamped=np.clip(y,interpolator.x[0],interpolator.x[-1])
        values=interpolator(clamped);derivative=interpolator.derivative()(clamped)
        low=y<=interpolator.x[0];high=y>=interpolator.x[-1]
        values[low]=ordinates[0];values[high]=ordinates[-1];derivative[low|high]=0
        return values,derivative

    def _field(self,positions):
        v=np.asarray(positions,float)
        if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all():raise ValueError('Finite physical positions with shape [N,3] required')
        y=v[:,self.longitudinal_axis];anchor,ap=self._values(self.a,y,self.av)
        negative,np_=self._values(self.n,y,self.negative);positive,pp=self._values(self.p,y,self.positive)
        mean=(positive+negative)/2;difference=(positive-negative)/2
        mean_p=(pp+np_)/2;difference_p=(pp-np_)/2
        q=v[:,self.transverse]-anchor;r=np.hypot(q,self.transition_width_m)
        # q**2/(r+epsilon) avoids cancellation near the fixed curve.
        smooth_abs=q*q/(r+self.transition_width_m)
        mapped=anchor+mean*q+difference*smooth_abs
        derivative=mean+difference*q/r
        longitudinal=(1-derivative)*ap+mean_p*q+difference_p*smooth_abs
        identity=(positive==1)&(negative==1)
        mapped[identity]=v[:,self.transverse][identity]
        derivative[identity]=1;longitudinal[identity]=0
        return v,mapped,derivative,longitudinal

    def map(self,positions):
        v,mapped,_,_=self._field(positions);result=v.copy();result[:,self.transverse]=mapped
        return result

    def jacobians(self,positions):
        v,_,derivative,longitudinal=self._field(positions);j=np.tile(np.eye(3),(len(v),1,1))
        for k,axis in enumerate(self.transverse):
            j[:,axis,axis]=derivative[:,k];j[:,axis,self.longitudinal_axis]=longitudinal[:,k]
        return j

    def evidence(self):
        lower=np.minimum(self.negative.min(0),self.positive.min(0))
        upper=np.maximum(self.negative.max(0),self.positive.max(0))
        return dict(continuous_positive_jacobian=True,transverse_derivative_lower_bounds=lower.tolist(),
            transverse_derivative_upper_bounds=upper.tolist(),continuous_determinant_lower_bound=float(np.prod(lower)),
            longitudinal_coordinate_literal=True,anchor_curve_fixed=True,anatomical_acceptance=False,
            discrete_mesh_quality_and_tube_sections_require_independent_checks=True)
