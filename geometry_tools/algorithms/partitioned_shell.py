"""Shared finite covering and internal tissue cells from explicit interfaces.

This module supplies no anatomy, thickness, partitions or automatic interior
caps. Callers must provide both reviewed envelopes, source face ownership and
every actual internal interface. Each interface is stored once and referenced
with opposite orientation by the adjoining cells. No anatomical source or
catalogue is written here.
"""
from dataclasses import dataclass
from collections import defaultdict
import numpy as np


@dataclass(frozen=True)
class SharedSurface:
    id: str
    faces: np.ndarray
    kind: str


@dataclass(frozen=True)
class SurfaceUse:
    surface: str
    side: int
    neighbor: str | None


@dataclass(frozen=True)
class PartitionedShell:
    positions: np.ndarray
    surfaces: tuple[SharedSurface, ...]
    cells: dict[str, tuple[SurfaceUse, ...]]
    outer_vertex_count: int

    def cell_faces(self, cell):
        by_id = {surface.id: surface for surface in self.surfaces}
        return np.concatenate([by_id[use.surface].faces if use.side == 1
                               else by_id[use.surface].faces[:, [0, 2, 1]]
                               for use in self.cells[cell]])

    def exterior_faces(self):
        return np.concatenate([surface.faces for surface in self.surfaces
                               if surface.kind == 'outer'])


def _immutable(value, dtype):
    result = np.array(value, dtype=dtype, copy=True)
    result.flags.writeable = False
    return result


def build_partitioned_shell(outer, inner, faces, owners, regions,
                            interior_points=(), interior_interfaces=()):
    """Construct cell references, never infer the anatomical interior split.

    Interface faces index the combined bank: outer vertices, corresponding
    inner vertices, then explicit interior points. Their winding points out
    of region `a` and into region `b`. Completeness/topology/containment must
    be independently checked before anatomical use. The returned cell names
    are `capsule:<region>` and `core:<region>`; these are geometric roles.
    """
    outer = np.asarray(outer)
    inner = np.asarray(inner)
    faces = np.asarray(faces)
    owners = np.asarray(owners)
    regions = tuple(regions)
    if (outer.ndim != 2 or outer.shape[1] != 3 or inner.shape != outer.shape
            or outer.dtype not in (np.dtype('float32'), np.dtype('float64'))
            or inner.dtype != outer.dtype
            or not np.isfinite(outer).all() or not np.isfinite(inner).all()):
        raise ValueError('Explicit finite corresponding envelopes required')
    if (faces.ndim != 2 or faces.shape[1] != 3 or not len(faces)
            or not np.issubdtype(faces.dtype, np.integer)
            or faces.min() < 0 or faces.max() >= len(outer)):
        raise ValueError('Invalid carrier faces')
    if (owners.shape != (len(faces),)
            or not np.issubdtype(owners.dtype, np.integer)
            or len(regions) < 1 or len(set(regions)) != len(regions)
            or not all(isinstance(r, str) and r for r in regions)
            or set(owners.tolist()) != set(range(len(regions)))):
        raise ValueError('Every source face must have one valid nonempty owner')
    extra = np.asarray(interior_points, dtype=outer.dtype).reshape(-1, 3)
    if not np.isfinite(extra).all():
        raise ValueError('Nonfinite internal interface points')
    bank = _immutable(np.concatenate((outer, inner, extra)), outer.dtype)
    count = len(outer)
    surfaces = []
    cells = {role+':'+r: [] for r in regions for role in ('capsule', 'core')}

    def surface(identifier, triangles, kind):
        triangles = _immutable(triangles, np.int64).reshape(-1, 3)
        if (not len(triangles) or triangles.min() < 0
                or triangles.max() >= len(bank)):
            raise ValueError('Empty or invalid shared interface')
        surfaces.append(SharedSurface(identifier, triangles, kind))

    def use(cell, identifier, side, neighbor=None):
        cells[cell].append(SurfaceUse(identifier, side, neighbor))

    # Original carrier face ownership changes no source coordinate/triangle.
    for index, region in enumerate(regions):
        owned = faces[owners == index]
        outer_id, inner_id = 'outer:'+region, 'inner:'+region
        surface(outer_id, owned, 'outer')
        surface(inner_id, owned+count, 'capsule-core')
        capsule, core = 'capsule:'+region, 'core:'+region
        use(capsule, outer_id, 1)
        use(capsule, inner_id, -1, core)
        use(core, inner_id, 1, capsule)

    edges = defaultdict(list)
    for owner, triangle in zip(owners, faces):
        for a,b in zip(triangle, np.roll(triangle, -1)):
            edges[tuple(sorted((int(a), int(b))))].append((int(owner), int(a), int(b)))
    rims = defaultdict(list)
    for adjacent in edges.values():
        if len(adjacent) != 2:
            raise ValueError('Whole carrier must be closed and manifold')
        left, right = adjacent
        if left[1] != right[2] or left[2] != right[1]:
            raise ValueError('Carrier winding must be consistent')
        if left[0] == right[0]:
            continue
        if left[0] > right[0]:
            left, right = right, left
        owner,a,b = left
        # Rim closes only the thickness between the paired envelopes.
        rims[(owner, right[0])].extend(((b,a,a+count), (b,a+count,b+count)))
    for (a,b), triangles in sorted(rims.items()):
        identifier = 'rim:'+regions[a]+':'+regions[b]
        surface(identifier, triangles, 'capsule-capsule')
        left, right = 'capsule:'+regions[a], 'capsule:'+regions[b]
        use(left, identifier, 1, right)
        use(right, identifier, -1, left)

    seen = set()
    for interface in interior_interfaces:
        a,b = interface['a'], interface['b']
        if (a not in regions or b not in regions or a == b
                or frozenset((a,b)) in seen):
            raise ValueError('Unique explicitly owned internal interface required')
        seen.add(frozenset((a,b)))
        identifier = 'internal:'+a+':'+b
        triangles = np.asarray(interface['faces'])
        if not np.issubdtype(triangles.dtype, np.integer):
            raise ValueError('Internal interface indices must be integers')
        if triangles.size and triangles.min() < count:
            raise ValueError('Internal interface must not reach the exterior')
        surface(identifier, triangles, 'core-core')
        left, right = 'core:'+a, 'core:'+b
        use(left, identifier, 1, right)
        use(right, identifier, -1, left)

    return PartitionedShell(bank, tuple(surfaces),
                            {key: tuple(value) for key,value in cells.items()}, count)


def validate_reciprocal_ownership(assembly):
    """Prove literal shared-interface identity and exact opposite orientation."""
    uses = defaultdict(list)
    for cell, refs in assembly.cells.items():
        for ref in refs:
            uses[ref.surface].append((cell, ref))
    for surface in assembly.surfaces:
        refs = uses[surface.id]
        expected = 1 if surface.kind == 'outer' else 2
        if len(refs) != expected:
            raise ValueError('Incorrect number of tissue owners')
        if expected == 1:
            if refs[0][1].side != 1 or refs[0][1].neighbor is not None:
                raise ValueError('Invalid exterior ownership')
            continue
        (a,aa),(b,bb) = refs
        if aa.side != -bb.side or aa.neighbor != b or bb.neighbor != a:
            raise ValueError('Nonreciprocal tissue boundary')
    return True
