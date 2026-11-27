"""Finite attachments between explicitly supplied physical surface patches.

The caller supplies both tissue contacts, a smooth course and section sizes.
This constructs no anatomical identity or attachment from proximity. Source
patch coordinates and triangulations remain literal; independent whole-cell,
host-contact and occupied-volume checks are still required.
"""
import numpy as np

from .annular_pipe_envelopes import boundary_loops, face_components
from .tube_sections import rotation_minimizing_frames
from .tube_rim_bridge import bridge_boundaries


def bridge_patches(first, second, centers, radii, *, first_side=-1,
                   second_side=-1, section_controls=12, projected_phase_alignment=True):
    """Bridge two explicit one-boundary patches with elliptical sections.

    Sides specify cap orientation relative to the supplied host face winding.
    A sheet reference can have a different side from a closed host volume;
    callers must record that distinction rather than assert a solid interface.
    Returned surface names distinguish the two literal caps and the free wall.
    Projected angular correspondence is optional for irregular source patches;
    topology and independent whole-surface intersection checks remain required.
    """
    if first.frame != second.frame:
        raise ValueError('Attachment patches require one explicitly common frame')
    if first_side not in (-1, 1) or second_side not in (-1, 1):
        raise ValueError('Both physical cap orientations must be explicit signs')
    centers, radii = np.asarray(centers, float), np.asarray(radii, float)
    if centers.ndim != 2 or centers.shape[1] != 3 or len(centers) < 4:
        raise ValueError('At least four finite course stations are required')
    if radii.shape != (len(centers), 2) or not np.isfinite(centers).all() \
            or not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError('Every station needs two finite positive section radii')
    if not isinstance(section_controls, int) or not 6 <= section_controls <= 64:
        raise ValueError('Section control count must be between six and sixty-four')
    if np.any(np.linalg.norm(np.diff(centers, axis=0), axis=1) == 0):
        raise ValueError('Consecutive course stations must be distinct')
    points, caps, loops, start = [], [], [], 0
    for mesh, side in [(first, first_side), (second, second_side)]:
        # Literal coordinate aliases may be shading/UV seams. This exact
        # identification changes no physical control or anatomical boundary.
        p, inverse = np.unique(mesh.triangles.reshape(-1, 3), axis=0, return_inverse=True)
        f = inverse.reshape(-1, 3)
        if side == -1:
            f = f[:, [0, 2, 1]]
        boundary = boundary_loops(f)
        if len(boundary) != 1 or face_components(f)[0] != 1:
            raise ValueError('Each contact must be one connected patch with one boundary')
        points.append(p); caps.append(f + start); loops.append(boundary[0] + start)
        start += len(p)
    for i, (patch, face) in enumerate(zip(points, caps)):
        local = face - (0 if i == 0 else len(points[0]))
        tri = patch[local]
        area = np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)
        if np.any(area == 0):
            raise ValueError('A contact patch contains a degenerate source triangle')
        center = np.einsum('i,ij->j', area, tri.mean(1)) / area.sum()
        if np.linalg.norm(center-centers[0 if i == 0 else -1]) > 1e-9:
            raise ValueError('Course ends must be the actual area centroids of the supplied patches')
        normal = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]).sum(0)
        direction = centers[1]-centers[0] if i == 0 else centers[-1]-centers[-2]
        if (normal @ direction) * (-1 if i == 0 else 1) <= 0:
            raise ValueError('A cap faces away from its declared physical end')
    frames = rotation_minimizing_frames(centers)
    n, b = frames[:, :, 0], frames[:, :, 1]
    angle = 2*np.pi*np.arange(section_controls)/section_controls
    rings = []
    for i in range(1, len(centers)-1):
        p = centers[i] + radii[i, 0]*np.cos(angle)[:, None]*n[i] \
            + radii[i, 1]*np.sin(angle)[:, None]*b[i]
        rings.append(np.arange(start, start+section_controls)); start += section_controls
        points.append(p)
    positions = np.concatenate(points)
    first_strip, first_evidence = bridge_boundaries(positions, loops[0], rings[0],align_projected_phase=projected_phase_alignment)
    last_strip, last_evidence = bridge_boundaries(positions, rings[-1][::-1], loops[1],align_projected_phase=projected_phase_alignment)
    strips = [first_strip]
    for a, b in zip(rings, rings[1:]):
        aj, bj = np.roll(a, -1), np.roll(b, -1)
        strips.extend([np.column_stack([a, aj, bj]), np.column_stack([a, bj, b])])
    strips.append(last_strip)
    surfaces = dict(first_contact=caps[0], free_wall=np.concatenate(strips), second_contact=caps[1])
    return positions, surfaces, dict(
        first_contact_triangles=len(caps[0]), second_contact_triangles=len(caps[1]),
        source_contact_coordinates_and_faces_literal=True,
        first_cap_side=first_side, second_cap_side=second_side,
        first_strip=first_evidence, last_strip=last_evidence,
        course_stations=len(centers), section_controls=section_controls,
        section_dimensions_are_explicit_caller_inputs=True,
        anatomy_inferred=False, whole_geometry_verified=False)
