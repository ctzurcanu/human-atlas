"""Add original University of Denver Visible Human Female lower-limb STLs.

Requires NumPy. Run after import-primary-female.mjs. The source ZIP
stays in .local-models/source and is never copied to public/models.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import zipfile

import numpy as np

SOURCE_URL = 'https://digitalcommons.du.edu/visiblehuman/1/'
SOURCE_SHA256 = '9886eda040f6087bbb65182530f8bd262be456f366b155025c12b3dc56da774e'
STL_DTYPE = np.dtype([('normal', '<f4', 3), ('vertices', '<f4', (3, 3)), ('attribute', '<u2')])
STL_NAME = re.compile(r'VHF_(Left|Right)_Muscle_([A-Za-z0-9]+)_smooth\.stl$')


def stl_vertices(archive, member):
    data = archive.read(member)
    count = int.from_bytes(data[80:84], 'little')
    if len(data) == 84 + 50 * count:
        return np.frombuffer(data, dtype=STL_DTYPE, count=count, offset=84)['vertices'].copy()
    if data.startswith(b'solid'):
        rows = re.findall(rb'^\s*vertex\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)', data, re.M)
        if not rows or len(rows) % 3:
            raise ValueError(f'Invalid ASCII STL vertices: {member}')
        return np.asarray(rows, dtype='<f4').reshape(-1, 3, 3)
    raise ValueError(f'Invalid STL length: {member}')


def hra_vertices(atlas, directory, name):
    part = next(p for p in atlas['parts'] if p['id'] == name)
    chunk = directory / Path(atlas['chunks'][part['chunk']]['url']).name
    with chunk.open('rb') as stream:
        stream.seek(part['positions'])
        return np.frombuffer(stream.read(part['vertexCount'] * 12), dtype='<f4').reshape(-1, 3).copy()


def nearest_distances(source, target):
    distances = []
    for start in range(0, len(source), 256):
        block = source[start:start + 256]
        distances.extend(np.sqrt(np.min(np.sum((block[:, None, :] - target[None, :, :]) ** 2, axis=2), axis=1)))
    return np.asarray(distances)


def fit_side(side, archive, atlas, directory):
    records = []
    members = archive.namelist()
    for bone in ('Femur', 'Tibia', 'Fibula'):
        suffix = f'VHF_{side}_Bone_{bone}_smooth.stl'
        member = next(n for n in members if n.endswith(suffix))
        source = np.unique(stl_vertices(archive, member).reshape(-1, 3), axis=0)
        target = hra_vertices(atlas, directory, f'HRA:VH_F_{bone.lower()}_{side[0]}')
        records.append((bone, source, target))

    # Axis correspondence follows the source landmarks: source Z is superior,
    # source Y is anterior, and source X is reflected to atlas left/right.
    # The atlas and STL bones differ in breadth, so solve each side separately.
    source_bounds = np.array([[s.min(0), s.max(0)] for _, s, _ in records])
    target_bounds = np.array([[t.min(0), t.max(0)] for _, _, t in records])
    source_width = source_bounds[:, 1] - source_bounds[:, 0]
    target_width = target_bounds[:, 1] - target_bounds[:, 0]
    x_scale = -np.sum(source_width[:, 0] * target_width[:, 0]) / np.sum(source_width[:, 0] ** 2)
    fit = np.zeros((3, 3), dtype=np.float64)
    offset = np.zeros(3, dtype=np.float64)
    fit[0, 0] = x_scale
    for source_axis, target_axis in ((2, 1), (1, 2)):
        x = source_bounds[:, :, source_axis].reshape(-1)
        y = target_bounds[:, :, target_axis].reshape(-1)
        fit[target_axis, source_axis], offset[target_axis] = np.linalg.lstsq(
            np.column_stack((x, np.ones(len(x)))), y, rcond=None
        )[0]
    source_mid = source_bounds.mean(axis=1)
    target_mid = target_bounds.mean(axis=1)
    fit[0, 2], offset[0] = np.linalg.lstsq(
        np.column_stack((source_mid[:, 2], np.ones(3))),
        target_mid[:, 0] - x_scale * source_mid[:, 0], rcond=None
    )[0]

    errors = []
    for bone, source, target in records:
        aligned = source[::max(1, len(source) // 5000)] @ fit.T + offset
        distances = nearest_distances(aligned, target)
        center_error = float(np.linalg.norm((source_mid[len(errors)] @ fit.T + offset) - target_mid[len(errors)]))
        errors.append({'bone': bone, 'centerMm': round(center_error * 1000, 2),
                       'medianSurfaceMm': round(float(np.median(distances)) * 1000, 2),
                       'p90SurfaceMm': round(float(np.percentile(distances, 90)) * 1000, 2)})
    if max(error['centerMm'] for error in errors) > 16:
        raise ValueError(f'{side} bone alignment exceeds 16 mm: {errors}')
    return fit, offset, errors


def display_name(raw, side):
    words = re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', raw)
    return f'{words} ({side.lower()})'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--zip', dest='source', default='.local-models/source/Final 3D STL Models-stl-female.zip')
    parser.add_argument('--atlas', default='public/models/atlas-hra-female.json')
    args = parser.parse_args()
    source = Path(args.source)
    directory = Path(args.atlas).resolve().parent
    atlas = json.loads(Path(args.atlas).read_text())
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError(f'University of Denver source hash differs: {digest}')
    if any(p['id'].startswith('VHF:') for p in atlas['parts']):
        raise ValueError('VHF muscles are already imported; rebuild the HRA base first')
    if len(atlas['parts']) != 956:
        raise ValueError('Expected 948 body meshes from v1.10 plus eight v1.5 pelvic meshes, after the placenta split')

    with zipfile.ZipFile(source) as archive:
        alignments = {side: fit_side(side, archive, atlas, directory) for side in ('Left', 'Right')}
        original_members = sorted(n for n in archive.namelist() if STL_NAME.search(n))
        if len(original_members) != 76:
            raise ValueError(f'Expected 76 original female muscle STLs, found {len(original_members)}')
        # HRA already supplies both rectus femoris muscles. Keep its versions
        # so the same anatomy is not drawn twice with overlapping surfaces.
        member_names = [n for n in original_members if not n.endswith('_Muscle_RectusFemoris_smooth.stl')]
        if len(member_names) != 74:
            raise ValueError('Expected two existing HRA rectus femoris meshes')

        segments = []
        length = 0
        added_triangles = 0
        def append(array):
            nonlocal length
            padding = (-length) % 4
            if padding:
                segments.append(bytes(padding))
                length += padding
            position = length
            data = array.tobytes()
            segments.append(data)
            length += len(data)
            return position

        def flush():
            nonlocal segments, length
            if not length:
                return
            name = f'female-muscle-{len(atlas["chunks"])}.bin'
            raw = b''.join(segments)
            compressed = gzip.compress(raw, compresslevel=9, mtime=0)
            (directory / name).write_bytes(raw)
            (directory / (name + '.gz')).write_bytes(compressed)
            atlas['chunks'].append({'url': '/models/' + name, 'bytes': len(raw),
                                    'gzip': '/models/' + name + '.gz', 'gzipBytes': len(compressed)})
            segments = []
            length = 0

        atlas['materials']['VHF:muscular'] = {'color': [165, 84, 79], 'sourceColor': [180, 94, 88]}
        for member in member_names:
            match = STL_NAME.search(member)
            side, raw_name = match.groups()
            fit, offset, _ = alignments[side]
            triangles = stl_vertices(archive, member)
            vertices, inverse = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
            positions = (vertices.astype(np.float64) @ fit.T + offset).astype('<f4')
            indices = inverse.reshape(-1, 3).astype('<u4')
            if np.linalg.det(fit) < 0:
                indices[:, [1, 2]] = indices[:, [2, 1]]
            vectors = positions[indices[:, 1]] - positions[indices[:, 0]]
            cross = np.cross(vectors, positions[indices[:, 2]] - positions[indices[:, 0]])
            normals = np.zeros_like(positions, dtype=np.float64)
            for corner in range(3):
                np.add.at(normals, indices[:, corner], cross)
            norms = np.linalg.norm(normals, axis=1)
            normals /= np.maximum(norms[:, None], 1e-20)
            packed_normals = np.round(np.clip(normals, -1, 1) * 32767).astype('<i2')
            if length > 6_000_000:
                flush()
            id_ = f'VHF:{side}:{raw_name}'
            name = display_name(raw_name, side)
            part = {'id': id_, 'conceptId': id_, 'name': name, 'system': 'muscular',
                    'chunk': len(atlas['chunks']), 'positions': append(positions),
                    'normals': append(packed_normals), 'indices': append(indices),
                    'vertexCount': len(positions), 'indexCount': indices.size,
                    'bounds': [positions.min(axis=0).tolist(), positions.max(axis=0).tolist()],
                    'material': 'VHF:muscular', 'sourceId': Path(member).name,
                    'license': 'CC BY 4.0',
                    'provenance': {'label': 'University of Denver Visible Human Female original STL',
                                   'url': SOURCE_URL, 'detail': f'Original archive member {member}; CC BY 4.0.'}}
            atlas['parts'].append(part)
            atlas['concepts'].append({'id': id_, 'name': name, 'elements': [id_]})
            added_triangles += indices.size
        flush()

    atlas['triangles'] += added_triangles // 3
    atlas['version'] += ' + original Visible Human Female lower-limb muscles'
    atlas['scope'] = 'Original HRA references with separately aligned original VHF lower-limb muscles; upper-body muscles remain partial.'
    atlas['provenance']['muscleSourceUrl'] = SOURCE_URL
    atlas['provenance']['muscleSourceSha256'] = digest
    atlas['provenance']['muscleSourceArchive'] = source.name
    atlas['provenance']['muscleCount'] = 74
    atlas['provenance']['muscleSourceTriangles'] = added_triangles // 3
    atlas['provenance']['muscleArchiveCount'] = 76
    atlas['provenance']['muscleOmittedAsHraDuplicates'] = ['VHF_Left_Muscle_RectusFemoris_smooth.stl', 'VHF_Right_Muscle_RectusFemoris_smooth.stl']
    atlas['provenance']['muscleAlignment'] = {
        side: {'matrix': fit.tolist(), 'offset': offset.tolist(), 'boneFit': errors}
        for side, (fit, offset, errors) in alignments.items()
    }
    atlas['provenance']['notes'] += ' Original University of Denver STL muscles are separately aligned to HRA femur, tibia, and fibula landmarks; source anatomy differs, and fit error is documented.'
    Path(args.atlas).write_text(json.dumps(atlas, separators=(',', ':')))
    print(json.dumps({'parts': len(atlas['parts']), 'concepts': len(atlas['concepts']),
                      'triangles': atlas['triangles'], 'alignment': atlas['provenance']['muscleAlignment']}, indent=2))


if __name__ == '__main__':
    main()
