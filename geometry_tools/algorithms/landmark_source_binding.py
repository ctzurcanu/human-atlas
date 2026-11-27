"""Reject pose reuse when landmark guides belong to another source geometry."""
import numpy as np


def verify_landmark_source(meshes, reference, *, coordinate_tolerance_m=0.):
    """Verify identities, frames, literal face ancestry and source coordinates.

    Frame names alone do not establish a common pose. A caller may explicitly
    allow a small known export-rounding tolerance; it is never inferred from
    anatomical proximity. This check does not prove landmark correspondence.
    """
    if not np.isfinite(coordinate_tolerance_m) or not 0 <= coordinate_tolerance_m <= 1e-5:
        raise ValueError('Source export tolerance must be explicit and at most 10 micrometres')
    rows = []
    for pid, mesh in meshes.items():
        if pid not in reference:
            raise ValueError('Landmark source reference is missing '+pid)
        source = reference[pid]
        if mesh.frame != source.frame:
            raise ValueError('Landmark source frame mismatch: '+pid)
        if mesh.vertices.shape != source.vertices.shape or not np.array_equal(mesh.faces, source.faces):
            raise ValueError('Landmark source face ancestry mismatch: '+pid)
        if not np.isfinite(mesh.vertices).all() or not np.isfinite(source.vertices).all():
            raise ValueError('Non-finite landmark source: '+pid)
        error = float(np.linalg.norm(mesh.vertices-source.vertices, axis=1).max(initial=0))
        if error > coordinate_tolerance_m:
            raise ValueError(f'Landmark source pose mismatch for {pid}: {error*1000:.6f} mm. '
                             'Use the original source or derive guides for this actual pose.')
        rows.append(dict(partId=pid, maximumCoordinateError_m=error))
    return dict(complete=True, coordinateTolerance_m=coordinate_tolerance_m,
                meshes=rows, framesAloneAccepted=False, anatomicalCorrespondenceProven=False)
