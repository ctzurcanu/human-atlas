"""Separate source-profile discrepancy from allowed physical section change.

The original reconciler and its historical hash-bound callers stay intact.
This wrapper applies a declared input correspondence budget, then enforces
the independently declared physical area-change limit on every output rim.
"""
from .declared_mouth_batch import reconcile_declared_mouths as _reconcile


def reconcile_declared_mouths(meshes,pairs,*,maximum_gap_m,maximum_motion_m,
                              maximum_relative_area_change,
                              maximum_source_area_discrepancy,
                              maximum_facet_rotation_degrees):
    if not 0<=maximum_relative_area_change<=maximum_source_area_discrepancy<1:
        raise ValueError('Explicit separate source and physical section-area budgets required')
    result,proof=_reconcile(meshes,pairs,maximum_gap_m=maximum_gap_m,
        maximum_motion_m=maximum_motion_m,
        maximum_relative_area_change=maximum_source_area_discrepancy,
        maximum_facet_rotation_degrees=maximum_facet_rotation_degrees)
    if any(max(p['relativeSectionAreaChanges'])>maximum_relative_area_change for p in proof['pairs']):
        raise ValueError('Output section-area change exceeds its unchanged physical budget')
    proof.update(maximumSourceAreaDiscrepancy=maximum_source_area_discrepancy,
                 maximumPhysicalSectionAreaChange=maximum_relative_area_change,
                 separateInputAndOutputAreaBudgets=True)
    return result,proof
