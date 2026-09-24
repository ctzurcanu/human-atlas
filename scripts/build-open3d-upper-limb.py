"""Append original Open 3D Model upper-limb geometry to a direct Z-Anatomy atlas.

Run in Blender background mode with scripts disabled. The publisher GLB is
decoded by Blender; no viewer exports, aliases, or catalogues are inputs.
"""
import argparse
import gzip
import json
import math
import re
import struct
import sys
from collections import defaultdict
from pathlib import Path

import bpy

parser = argparse.ArgumentParser()
parser.add_argument('--glb', type=Path, required=True)
parser.add_argument('--atlas', type=Path, required=True)
parser.add_argument('--zip-sha256', required=True)
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
atlas = json.loads(args.atlas.read_text())
assert atlas['source'] == 'Z-Anatomy'
assert not any(p['id'].startswith('O3M:') for p in atlas['parts'])
assert re.fullmatch('[0-9a-f]{64}', args.zip_sha256)

source_url = 'https://anatomytool.org/open3dmodel-create'
source_glb_url = 'https://caskanatomy.info/open3dmodelfiles/upper-limb/upper-limb-glb.zip'
existing = {p['sourceId'].strip().casefold() for p in atlas['parts']}
existing_ids = {p['id'] for p in atlas['parts']}
bpy.ops.import_scene.gltf(filepath=str(args.glb))
new_chunks = []
blob = bytearray()
added = []
concepts = defaultdict(list)
triangles = 0


def source_position(v):
    # Blender's glTF importer converts glTF Y-up to Blender Z-up.
    return (v.x, v.z, -v.y)


def append(fmt, values):
    while len(blob) % 4:
        blob.append(0)
    offset = len(blob)
    blob.extend(struct.pack('<' + str(len(values)) + fmt, *values))
    return offset


def flush():
    global blob
    if not blob:
        return
    name = f'o3m-upper-limb-{len(new_chunks)}.bin'
    raw = bytes(blob)
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    (args.atlas.parent / name).write_bytes(raw)
    (args.atlas.parent / (name + '.gz')).write_bytes(packed)
    new_chunks.append({'url': '/models/' + name, 'bytes': len(raw), 'gzip': '/models/' + name + '.gz', 'gzipBytes': len(packed)})
    blob = bytearray()


def system_for(group):
    suffix = group.rsplit(' - ', 1)[-1]
    return {
        'bones': 'skeletal', 'muscles': 'muscular', 'arteries': 'arterial',
        'veins': 'venous', 'nerves': 'nervous', 'cartilages': 'connective',
        'capsules, ligaments, fasciae': 'connective', 'synovia, bursae': 'connective',
    }.get(suffix)


palette = {
    'skeletal': [217, 207, 184], 'muscular': [163, 82, 78],
    'arterial': [188, 76, 76], 'venous': [104, 114, 164],
    'nervous': [225, 188, 104], 'connective': [181, 189, 178],
}


def included(group, name):
    if name.casefold() in existing:
        return False
    if group.startswith(('Back - ', 'Head and neck - ', 'Thorax - bones', 'Thorax - cartilages')):
        return False
    if not group.startswith(('Arm - ', 'Forearm - ', 'Hand and wrist - ', 'Pectoral girdle - ', 'Thorax - ')):
        return False
    if system_for(group) is None:
        return False
    # The publisher's right-sided hand has many unsuffixed ligaments, sheaths,
    # and neurovascular structures. Its joint cartilage and one mislabeled
    # ligament in the muscle group are left out to avoid duplicate surfaces.
    if group == 'Hand and wrist - cartilages' or name == 'Superficial transverse metacarpal lig.':
        return False
    return name.endswith('.r') or group.startswith('Hand and wrist - ')


for obj in sorted((o for o in bpy.data.objects if o.type == 'MESH' and o.parent), key=lambda o: o.name):
    group = obj.parent.name
    name = obj.name
    if not included(group, name):
        continue
    system = system_for(group)
    mesh = obj.data
    mesh.calc_loop_triangles()
    by_material = defaultdict(list)
    for face in mesh.loop_triangles:
        by_material[mesh.polygons[face.polygon_index].material_index].append(face)
    for mirrored in (False, True):
        side_name = (name[:-2] + ('.l' if mirrored else '.r')) if not name.endswith('.r') else (name[:-2] + '.l' if mirrored else name)
        display_name = re.sub(r'\.[rl]$', lambda m: ' (left)' if m.group() == '.l' else ' (right)', side_name)
        for slot, faces in sorted(by_material.items()):
            if not faces:
                continue
            mat = obj.material_slots[slot].material if slot < len(obj.material_slots) else None
            mat_name = mat.name if mat else 'Unassigned'
            part_id = 'O3M:' + side_name + (':' + mat_name if len(by_material) > 1 else '')
            assert part_id not in existing_ids, part_id
            existing_ids.add(part_id)
            vertices, normals, indices = [], [], []
            lookup = {}
            minimum, maximum = [math.inf] * 3, [-math.inf] * 3
            transform = obj.matrix_world
            normal_transform = transform.to_3x3().inverted_safe().transposed()
            reflected = transform.determinant() < 0
            for face in faces:
                corners = list(face.vertices)
                if reflected != mirrored:
                    corners[1], corners[2] = corners[2], corners[1]
                for src in corners:
                    target = lookup.get(src)
                    if target is None:
                        p = list(source_position(transform @ mesh.vertices[src].co))
                        n = list(source_position((normal_transform @ mesh.vertices[src].normal).normalized()))
                        if mirrored:
                            p[0], n[0] = -p[0], -n[0]
                        target = len(lookup)
                        lookup[src] = target
                        vertices.extend(p)
                        normals.extend(max(-32767, min(32767, round(v * 32767))) for v in n)
                        for axis in range(3):
                            minimum[axis] = min(minimum[axis], p[axis])
                            maximum[axis] = max(maximum[axis], p[axis])
                    indices.append(target)
            if len(blob) > 6_000_000:
                flush()
            material_id = 'O3M:' + system + ':' + mat_name
            base_color = mat.diffuse_color if mat else (0.7, 0.7, 0.7, 1)
            source_color = [round(255 * max(0, min(1, c))) for c in base_color[:3]]
            atlas['materials'][material_id] = {'color': palette[system], 'sourceColor': source_color}
            part = {
                'id': part_id, 'conceptId': 'O3M:' + side_name,
                'name': display_name + (f' · {mat_name}' if len(by_material) > 1 else ''),
                'system': system, 'chunk': len(atlas['chunks']) + len(new_chunks),
                'positions': append('f', vertices), 'normals': append('h', normals),
                'indices': append('I', indices), 'vertexCount': len(lookup),
                'indexCount': len(indices), 'bounds': [minimum, maximum],
                'sourceId': name, 'material': material_id, 'groups': [group],
                'license': 'CC BY-SA 4.0',
                'provenance': {
                    'label': 'Open 3D Model original upper-limb GLB', 'url': source_url,
                    'detail': 'Publisher mesh and material; left side mirrored from the right source' if mirrored else 'Publisher right-side mesh and material',
                    'mirrored': mirrored,
                },
            }
            atlas['parts'].append(part)
            added.append(part)
            concepts['O3M:' + side_name].append(part_id)
            triangles += len(indices) // 3

flush()
atlas['chunks'].extend(new_chunks)
for concept_id, members in concepts.items():
    name = concept_id[4:]
    display_name = re.sub(r'\.[rl]$', lambda m: ' (left)' if m.group() == '.l' else ' (right)', name)
    atlas['concepts'].append({'id': concept_id, 'name': display_name, 'elements': members})
atlas['triangles'] += triangles
atlas['version'] = 'Z-Anatomy + Open 3D Model original publisher sources'
atlas['source'] = 'Z-Anatomy + Open 3D Model'
atlas['scope'] = 'Direct Z-Anatomy body plus original Open 3D Model upper-limb additions; left-only additions are mirrored and marked.'
atlas['provenance']['additionalSources'] = [{
    'url': source_url, 'archiveUrl': source_glb_url,
    'archiveSha256': args.zip_sha256, 'sourceFile': 'upper-limb.glb',
    'sourceFileSha256': __import__('hashlib').sha256(args.glb.read_bytes()).hexdigest(),
    'sourceMeshNodes': 532, 'embeddedImages': 84,
    'notes': 'The original GLB embeds image textures. The browser atlas currently displays its source material colors through the common anatomy palette; source images are not rendered.',
}]
args.atlas.write_text(json.dumps(atlas, separators=(',', ':')))
print('HUMAN_ATLAS_O3M_BUILD=' + json.dumps({'addedSurfaces': len(added), 'addedConcepts': len(concepts), 'addedTriangles': triangles, 'newChunks': len(new_chunks), 'totalSurfaces': len(atlas['parts'])}))
