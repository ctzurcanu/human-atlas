"""Pack explicitly owned finite tissues with one shared attribute bank.

This module invents no tissue, closure, material or anatomical identity. Every
surface and outward corner field must already have been constructed and
reviewed. Reciprocal cells reuse the same indices/words with opposite signs.
"""
from collections import defaultdict
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PackedFiniteCells:
    data: bytes
    positions: np.ndarray
    normal_words: np.ndarray
    position_offset: int
    normal_offset: int
    cells: dict
    surfaces: dict


def pack_finite_cells(positions, surfaces, corner_words, cells, part_ids, geometry_origin=(0.,0.,0.), *, corner_uvs=None,share_orientation_buffers=True):
    origin=np.asarray(geometry_origin,dtype=np.float64)
    if origin.shape!=(3,) or not np.isfinite(origin).all():raise ValueError('A finite three-coordinate origin is required')
    positions = np.asarray(np.asarray(positions,dtype=np.float64)-origin, dtype='<f4')
    if positions.ndim != 2 or positions.shape[1] != 3 or not np.isfinite(positions).all():
        raise ValueError('A finite shared control bank is required')
    if set(surfaces) != set(corner_words) or set(cells) != set(part_ids) or len(set(part_ids.values())) != len(part_ids):
        raise ValueError('Complete unique surface fields and cell identities are required')
    if corner_uvs is not None and set(corner_uvs)!=set(surfaces):
        raise ValueError('Every surface in a textured chart needs explicit corner UVs')
    owners = defaultdict(list)
    for cell, references in cells.items():
        if not references: raise ValueError('Empty finite cell')
        for ref in references:
            if ref['surface'] not in surfaces or ref['side'] not in [-1, 1]:
                raise ValueError('Unknown surface or invalid orientation')
            owners[ref['surface']].append((cell, ref))
    for sid in surfaces:
        uses = owners[sid]
        if len(uses) == 1:
            if uses[0][1].get('neighbor') is not None: raise ValueError('Unreciprocated interface')
        elif len(uses) == 2:
            (a, ra), (b, rb) = uses
            if a == b or ra['side'] != -rb['side'] or ra.get('neighbor') != b or rb.get('neighbor') != a:
                raise ValueError('A shared surface needs two distinct opposed owners')
        else: raise ValueError('Unused or multiply duplicated surface')

    vertices, words, lookup, mapped, uv_rows = [], [], {}, {}, []
    for sid, triangles in surfaces.items():
        triangles = np.asarray(triangles)
        field = np.asarray(corner_words[sid])
        uv=np.asarray(corner_uvs[sid],dtype='<f4')if corner_uvs is not None else None
        if triangles.ndim != 2 or triangles.shape[1] != 3 or not len(triangles) or not np.issubdtype(triangles.dtype, np.integer):
            raise ValueError('Each surface requires actual indexed triangles')
        if triangles.min() < 0 or triangles.max() >= len(positions) or field.shape != (len(triangles), 3, 3) or field.dtype != np.dtype('<i2'):
            raise ValueError('Invalid surface index or packed corner field')
        if uv is not None and (uv.shape!=(len(triangles),3,2)or not np.isfinite(uv).all()):
            raise ValueError('Finite corner UV pairs must cover every actual surface triangle')
        xyz = positions[triangles].astype(float)
        fn = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
        if not np.all(np.einsum('tki,ti->tk', field.astype(float), fn) > 0):
            raise ValueError('Every emitted corner must be nonzero and outward after float32 packing')
        output = np.empty_like(triangles, dtype='<u4')
        for ti, triangle in enumerate(triangles):
            for k, control in enumerate(triangle):
                word = field[ti, k]
                # A chart boundary remains explicit, preventing signed-normal
                # conflicts when one cell owns several differently oriented
                # interfaces at the same geometric control.
                key = (sid, int(control), tuple(map(int, word)))
                if uv is not None:key=(*key,uv[ti,k].tobytes())
                if key not in lookup:
                    lookup[key] = len(vertices); vertices.append(positions[control]); words.append(word)
                    if uv is not None:uv_rows.append(uv[ti,k])
                output[ti, k] = lookup[key]
        mapped[sid] = output
    vertices = np.asarray(vertices, dtype='<f4'); words = np.asarray(words, dtype='<i2')
    buffer = bytearray()
    def append(value):
        buffer.extend(b'\0'*((-len(buffer))%4)); offset = len(buffer); buffer.extend(value.tobytes()); return offset
    po, no = append(vertices), append(words)
    uo=append(np.asarray(uv_rows,dtype='<f4'))if corner_uvs is not None else None
    # Only neighbours sharing a chart require different orientation buffers.
    # Nonadjacent cells can reuse one buffer: each owns disjoint render
    # vertices because the vertex key includes the chart identity. This keeps
    # regional peeling from multiplying a whole-bank byte array per cell.
    adjacency={cell:set()for cell in cells}
    for uses in owners.values():
        if len(uses)==2:
            a,b=uses[0][0],uses[1][0];adjacency[a].add(b);adjacency[b].add(a)
    colors={}if share_orientation_buffers else{cell:i for i,cell in enumerate(cells)}
    while len(colors)<len(cells):
        remaining=set(cells)-set(colors)
        cell=max(remaining,key=lambda c:(len({colors[n]for n in adjacency[c]if n in colors}),len(adjacency[c]),c))
        forbidden={colors[n]for n in adjacency[cell]if n in colors};color=0
        while color in forbidden:color+=1
        colors[cell]=color
    shared_signs={color:np.zeros(len(vertices),dtype='i1')for color in set(colors.values())}
    packed = {}
    for cell, references in cells.items():
        indices, groups, signs, start = [], [], shared_signs[colors[cell]], 0
        for ref in references:
            triangles = mapped[ref['surface']]
            if ref['side'] == -1: triangles = triangles[:, [0, 2, 1]]
            used = np.unique(triangles)
            if not np.all((signs[used] == 0)|(signs[used] == ref['side'])):
                raise ValueError('A render chart has conflicting ownership signs')
            signs[used] = ref['side']; indices.append(triangles.ravel())
            neighbor = ref.get('neighbor')
            groups.append(dict(start=start, count=int(triangles.size), **({'neighbor':part_ids[neighbor]} if neighbor else {})))
            start += int(triangles.size)
        indices = np.concatenate(indices).astype('<u4')
        io=append(indices);xyz = vertices[np.unique(indices)]
        packed[cell] = dict(id=part_ids[cell], positions=po, normals=no, indices=io,
                            vertexCount=len(vertices), indexCount=len(indices), triangles=len(indices)//3,
                            bounds=[(xyz.min(0).astype(float)+origin).tolist(), (xyz.max(0).astype(float)+origin).tolist()], geometryOrigin=origin.tolist(),
                            tissueDrawGroups=groups)
        if uo is not None:packed[cell]['uvs']=uo
    offsets={}
    for color,signs in sorted(shared_signs.items()):
        signs[signs==0]=1;offsets[color]=append(signs)
    for cell,record in packed.items():record['normalSigns']=offsets[colors[cell]]
    return PackedFiniteCells(bytes(buffer), vertices, words, po, no, packed, mapped)
