"""Caller-bounded geometric co-move pairs, explicitly not source joins."""
import numpy as np
from geometry_tools.source_port_dependencies import check_source_mouth


def check_offset_pair(first, second, a, b, declaration):
    """Validate complete boundaries without changing the source-join gate."""
    try:
        gap=float(declaration['max_source_gap_m'])
        area=float(declaration['max_source_relative_area_error'])
    except (KeyError,TypeError,ValueError) as exc:
        raise ValueError('Offset pairs require explicit max_source_gap_m and max_source_relative_area_error') from exc
    if not np.isfinite([gap,area]).all() or gap<0 or area<0:
        raise ValueError('Offset-pair source bounds must be finite and nonnegative')
    # This existing call still uses its original strict join thresholds.
    # Only its literal complete-boundary geometry is used for co-movement.
    p=check_source_mouth(first,second,a,b)
    if (p['maximumGap_m']>gap or p['relativeAreaError']>area or p['normalCosine']>=-.99):
        raise ValueError('Offset-pair source geometry exceeds caller bounds or opposed-normal requirement')
    return dict(geometricOnly=True,unjoined=True,sourceJoinPassed=False,
                completeSourceBoundaries=True,maximumSourceGap_m=p['maximumGap_m'],
                sourceRelativeAreaError=p['relativeAreaError'],normalCosine=p['normalCosine'],
                max_source_gap_m=gap,max_source_relative_area_error=area,
                originalOffsetVectors_m=(second.vertices[b]-first.vertices[a]).tolist(),
                anatomicalAttachmentEstablished=False,wholeLumenPatencyEstablished=False)
