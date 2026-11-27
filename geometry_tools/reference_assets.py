"""Read-only CID -> supplied assets.json -> exact GLB-node resolution.

Only requested accessors are read. This helper neither registers nor exports
geometry. Source coordinates/units are not assumed to match an installed atlas.
"""
from pathlib import Path
import hashlib
import json
import struct
import numpy as np


class ReferenceAssets:
    def __init__(self, assets, glb, descriptions=None):
        self.assets_path, self.glb_path = Path(assets), Path(glb)
        raw = self.assets_path.read_bytes()
        self.assets = json.loads(raw)
        self.assets_sha256 = hashlib.sha256(raw).hexdigest()
        self.names = {} if descriptions is None else {int(r['cid']): r.get('name') for r in json.loads(Path(descriptions).read_bytes())['items']}
        with self.glb_path.open('rb') as stream:
            header = stream.read(12)
            magic, version, length = struct.unpack('<III', header)
            assert magic == 0x46546c67 and version == 2 and length == self.glb_path.stat().st_size
            count, kind = struct.unpack('<II', stream.read(8))
            assert kind == 0x4e4f534a
            source_json = stream.read(count)
            self.gltf = json.loads(source_json)
            binary_length, kind = struct.unpack('<II', stream.read(8))
            assert kind == 0x004e4942
            self.binary_offset, self.binary_length = stream.tell(), binary_length
        self.json_sha256 = hashlib.sha256(source_json).hexdigest()
        self.read_ranges = {}
        self.parents = {}
        for i, node in enumerate(self.gltf['nodes']):
            for child in node.get('children', []):
                assert child not in self.parents, 'Multiple-parent nodes require a scene-specific path'
                self.parents[child] = i

    def resolve(self, cids):
        wanted = set(map(int, cids))
        rows = []
        for i, structure in enumerate(self.assets['structures']):
            if structure.get('cid') not in wanted:
                continue
            for attachment in structure.get('sni', []):
                index, expected = int(attachment['gltf']), attachment['id']
                node = self.gltf['nodes'][index]
                assert node.get('name') == expected, ('CID/node join mismatch', structure['cid'], index, expected, node.get('name'))
                rows.append(dict(cid=structure['cid'], name=self.names.get(structure['cid']), structure_index=i,
                                 structure_id=structure['id'], attr=structure.get('attr'), twin=structure.get('twin'),
                                 node_index=index, node_id=expected, mesh_index=node.get('mesh'),
                                 node_transform=self.world_matrix(index).tolist()))
        return rows

    def world_matrix(self, index):
        node = self.gltf['nodes'][index]
        if 'matrix' in node:
            matrix = np.asarray(node['matrix'], float).reshape(4, 4).T
        else:
            x, y, z, w = node.get('rotation', [0, 0, 0, 1])
            rotation = np.array([[1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
                                 [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
                                 [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]])
            matrix = np.eye(4)
            matrix[:3, :3] = rotation @ np.diag(node.get('scale', [1, 1, 1]))
            matrix[:3, 3] = node.get('translation', [0, 0, 0])
        return self.world_matrix(self.parents[index]) @ matrix if index in self.parents else matrix

    def accessor(self, index):
        accessor = self.gltf['accessors'][index]
        assert 'sparse' not in accessor, 'Sparse accessor not implemented'
        view = self.gltf['bufferViews'][accessor['bufferView']]
        assert view.get('buffer', 0) == 0
        dtype = np.dtype({5120:'i1', 5121:'u1', 5122:'<i2', 5123:'<u2', 5125:'<u4', 5126:'<f4'}[accessor['componentType']])
        width = {'SCALAR':1, 'VEC2':2, 'VEC3':3, 'VEC4':4, 'MAT4':16}[accessor['type']]
        stride = view.get('byteStride', width*dtype.itemsize)
        count = accessor['count']
        offset = self.binary_offset+view.get('byteOffset', 0)+accessor.get('byteOffset', 0)
        length = (count-1)*stride+width*dtype.itemsize if count else 0
        with self.glb_path.open('rb') as stream:
            stream.seek(offset)
            raw = stream.read(length)
        assert len(raw) == length
        self.read_ranges[index] = dict(accessor=index, byte_offset=offset, byte_length=length, sha256=hashlib.sha256(raw).hexdigest())
        result = np.ndarray((count, width), dtype=dtype, buffer=raw, strides=(stride, dtype.itemsize)).copy()
        assert not accessor.get('normalized', False), 'Normalized geometry accessor requires explicit handling'
        return result

    def node_mesh(self, row):
        vertices, faces, offset = [], [], 0
        node = self.gltf['nodes'][row['node_index']]
        assert 'skin' not in node, 'Skinned geometry requires evaluated bind transforms'
        for primitive in self.gltf['meshes'][row['mesh_index']]['primitives']:
            assert not primitive.get('extensions') and not primitive.get('targets')
            points = self.accessor(primitive['attributes']['POSITION']).astype(float)
            indices = self.accessor(primitive['indices']).ravel().astype(np.int64) if 'indices' in primitive else np.arange(len(points))
            mode = primitive.get('mode', 4)
            if mode == 4:
                triangles = indices.reshape(-1, 3)
            elif mode == 5:
                triangles = np.array([indices[j:j+3] if j%2 == 0 else indices[j:j+3][[1,0,2]] for j in range(len(indices)-2)])
            elif mode == 6:
                triangles = np.array([[indices[0], indices[j], indices[j+1]] for j in range(1, len(indices)-1)])
            else:
                raise ValueError('Requested mesh is not a triangle primitive')
            triangles = triangles[(triangles[:,0] != triangles[:,1]) & (triangles[:,1] != triangles[:,2]) & (triangles[:,2] != triangles[:,0])]
            matrix = np.asarray(row['node_transform'])
            vertices.append(points @ matrix[:3,:3].T+matrix[:3,3])
            faces.append(triangles+offset)
            offset += len(points)
        return np.concatenate(vertices), np.concatenate(faces)

    def provenance(self):
        return dict(assets_path=str(self.assets_path), assets_sha256=self.assets_sha256,
                    glb_path=str(self.glb_path), glb_bytes=self.glb_path.stat().st_size,
                    glb_json_sha256=self.json_sha256, binary_offset=self.binary_offset, binary_bytes=self.binary_length,
                    requested_accessor_ranges=list(self.read_ranges.values()), complete_bank_sha256=None,
                    geometry_exported=False, registration_assumed=False)
