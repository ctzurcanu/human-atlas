"""Literal non-regression gate for a partial superficial-vessel repair."""
import numpy as np


def skin_transition_guard(old_distance, old_known, new_distance, new_known,
                          old_pairs, new_pairs, *, exterior_tolerance_m=1e-7):
    od,ok,nd,nk=[np.asarray(x)for x in [old_distance,old_known,new_distance,new_known]]
    if not(od.shape==ok.shape==nd.shape==nk.shape)or od.ndim!=1:
        raise ValueError('Matching per-control source classifications required')
    if ok.dtype.kind!='b'or nk.dtype.kind!='b'or not np.isfinite(od).all()or not np.isfinite(nd).all():
        raise ValueError('Finite distances and literal Boolean classifications required')
    old=set(map(tuple,old_pairs));new=set(map(tuple,new_pairs));ex_old=ok&(od>exterior_tolerance_m);ex_new=nk&(nd>exterior_tolerance_m)
    lost=int((ok&~nk).sum());added_controls=np.flatnonzero(ex_new&~ex_old);worsened=np.flatnonzero(ex_new&(nd>od+exterior_tolerance_m));added_pairs=sorted(new-old)
    improvement=bool(ex_new.sum()<ex_old.sum()or len(new)<len(old))
    passed=bool(not lost and not len(added_controls)and not len(worsened)and not added_pairs and improvement)
    return dict(passed=passed,strictImprovement=improvement,oldKnownBecomesUnknown=lost,
                beforeExteriorControls=int(ex_old.sum()),afterExteriorControls=int(ex_new.sum()),
                addedExteriorControls=added_controls.tolist(),worsenedExteriorControls=worsened.tolist(),
                beforeSkinTrianglePairs=len(old),afterSkinTrianglePairs=len(new),addedSkinTrianglePairs=added_pairs,
                anatomicalAcceptance=False)
