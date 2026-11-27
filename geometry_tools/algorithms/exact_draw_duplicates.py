"""Remove only literally identical oriented triangles with identical fields.

No proximity welding, vertex removal, normal reconstruction or opposite-sheet
merging. Differently shaded or textured coincident faces remain explicit.
"""
import numpy as np


def exact_draw_duplicates(mesh):
    _, physical = np.unique(mesh.vertices, axis=0, return_inverse=True)
    faces = physical[mesh.faces]
    _, groups, counts = np.unique(np.sort(faces, axis=1), axis=0,
                                  return_inverse=True, return_counts=True)
    keep = np.ones(len(faces), bool); removed = []; conflicts = []
    for group in np.flatnonzero(counts > 1):
        ids = np.flatnonzero(groups == group); representatives = []
        for face_id in ids:
            matched = None
            for prior in representatives:
                for shift in range(3):
                    order = np.roll(np.arange(3), shift)
                    if not np.array_equal(faces[face_id, order], faces[prior]):
                        continue
                    a, b = mesh.faces[face_id, order], mesh.faces[prior]
                    if all(np.array_equal(field[a], field[b]) for field in mesh.attributes.values()):
                        matched = (prior, order); break
                if matched is not None: break
            if matched is None:
                representatives.append(int(face_id))
            else:
                keep[face_id] = False
                removed.append(dict(removedFace=int(face_id), retainedSourceFace=matched[0],
                                    sourceCornerOrder=matched[1].tolist()))
        if len(representatives) > 1:
            conflicts.append(dict(sourceFaces=ids.tolist(), retainedDistinctFields=representatives))
    kept = np.flatnonzero(keep)
    return kept, dict(complete=True, retainedSourceFaceIds=kept.tolist(),
                     duplicateCopiesRemoved=removed, retainedAppearanceOrWindingConflicts=conflicts,
                     rawCoordinatesAndAttributesLiteral=True, sourceSurfaceSupportUnchanged=True,
                     newTrianglesCreated=0, vertexWeldingPerformed=False,
                     anatomicalAcceptance=False)
