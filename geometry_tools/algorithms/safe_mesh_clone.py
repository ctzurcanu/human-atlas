"""Fresh same-topology mesh poses without inherited derived-property caches."""
from collections.abc import Mapping
from copy import deepcopy

import numpy as np

from geometry_tools.model import Mesh


def clone_with_geometry(mesh, vertices=None, *, attribute_updates=None):
    """Copy an explicit source bank into a fresh ``Mesh`` with a new pose.

    Omitted attributes remain literal. Frame, faces and source ancestry are
    preserved, so changing the vertex-bank shape requires a separate explicit
    reconstruction. Mesh's constructor copies and freezes all stored arrays;
    no cached hashes, bounds, triangles or other instance state are copied.
    """
    if not isinstance(mesh, Mesh):
        raise TypeError('An existing Mesh is required')
    points = mesh.vertices if vertices is None else np.asarray(vertices)
    if points.shape != mesh.vertices.shape:
        raise ValueError('A same-topology pose must retain the source vertex-bank shape')
    if attribute_updates is not None and not isinstance(attribute_updates, Mapping):
        raise TypeError('Attribute updates must be an explicit mapping')
    attributes = dict(mesh.attributes)
    for name, values in (attribute_updates or {}).items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError('Attribute names must be nonempty strings')
        attributes[name] = values
    return Mesh(points, mesh.faces, mesh.frame, attributes, deepcopy(mesh.source))
