"""Build browser atlas assets directly from the official Z-Anatomy Startup.blend.

Run with Blender 4.3+ (the source file is never executed):
  blender -b Startup.blend --factory-startup --disable-autoexec --python scripts/build-z-anatomy.py -- --out /tmp/z-anatomy-atlas --source-sha256 SHA256

Use scripts/fetch-primary-models.mjs to obtain and verify the original archive.
This importer reads Blender mesh and curve objects, collections, and materials. It does
not use any third-party GLB export, catalog, alias table, or rendered asset.
"""
import argparse
import gzip
import hashlib
import json
import math
import re
import struct
import sys
from collections import defaultdict
from pathlib import Path

import bpy
from mathutils import Vector


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--source-sha256', required=True)
    parser.add_argument('--name', default='atlas-z-anatomy')
    return parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])


ARGS = arguments()
assert re.fullmatch(r'[0-9a-f]{64}', ARGS.source_sha256), 'Expected original archive SHA-256'
ARGS.out.mkdir(parents=True, exist_ok=True)

# The source is in meters, Blender Z-up, with the anatomical front toward -Y.
# Three.js uses Y-up and an anatomical front toward +Z.
def position(v):
    return (v.x, v.z, -v.y)


def source_color(material):
    if material is None:
        return [190, 190, 190]
    nodes = material.node_tree.nodes if material.use_nodes and material.node_tree else []
    gamma = next((n for n in nodes if n.type == 'GAMMA' and 'Color' in n.inputs), None)
    if gamma:
        color = gamma.inputs['Color'].default_value
    else:
        principled = next((n for n in nodes if n.type == 'BSDF_PRINCIPLED' and not n.inputs['Base Color'].is_linked), None)
        if principled:
            color = principled.inputs['Base Color'].default_value
        else:
            group = next((n for n in nodes if n.type == 'GROUP' and 'Color' in n.inputs and not n.inputs['Color'].is_linked), None)
            color = group.inputs['Color'].default_value if group else material.diffuse_color
    # Blender's material values are linear RGB; the viewer expects sRGB bytes.
    def srgb(value):
        value = max(0.0, min(1.0, float(value)))
        return round(255 * (12.92 * value if value <= 0.0031308 else 1.055 * value ** (1 / 2.4) - 0.055))
    return [srgb(v) for v in color[:3]]


def display_color(system, material_name, original):
    mat = material_name.lower()
    if system == 'integumentary':
        return original if any(term in mat for term in ('hair', 'black', 'cornea', 'iris', 'pupil')) else [191, 151, 126]
    if 'tendon' in mat:
        return [218, 195, 160]
    if 'cartilage' in mat:
        return [177, 191, 193]
    if 'ligament' in mat:
        return [211, 189, 164]
    palette = {
        'skeletal': [217, 207, 184], 'muscular': [163, 82, 78],
        'cardiac': [169, 79, 78], 'arterial': [188, 76, 76],
        'venous': [104, 114, 164], 'nervous': [225, 188, 104],
        'sensory': [191, 169, 145], 'respiratory': [166, 137, 147],
        'digestive': [166, 105, 95], 'urinary': [176, 122, 102],
        'lymphatic': [139, 166, 130], 'endocrine': [181, 129, 142],
        'reproductive': [174, 115, 114], 'connective': [181, 189, 178],
        'fascia': [202, 183, 169], 'attachments': [198, 149, 124],
    }
    return palette[system]


def system_for(obj, material_name):
    groups = {c.name for c in obj.users_collection}
    group_text = ' '.join(g for g in groups if not re.match(r'^\d+:', g)).lower()
    name = obj.name.lower()
    mat = material_name.lower()
    if '2: Muscular insertions' in groups:
        return 'attachments'
    if '9: Regions of human body' in groups or mat.startswith('skin'):
        return 'integumentary'
    if '1: Skeletal system' in groups:
        return 'skeletal'
    if '3: Joints' in groups:
        return 'connective'
    if '4: Muscular system' in groups:
        return 'fascia' if 'fascia' in mat or 'fascia' in name else 'muscular'
    if '5: Cardiovascular system' in groups:
        return 'venous' if 'vein' in mat or 'vein' in name else 'arterial' if 'arter' in mat or 'arter' in name else 'cardiac'
    if '6: Lymphoid organs' in groups:
        return 'lymphatic'
    if '7: Nervous system & Sense organs' in groups:
        return 'sensory' if any(word in group_text for word in ('eye', 'ear', 'olfact', 'vestibul', 'retina', 'lens')) else 'nervous'
    if '8: Visceral systems' in groups:
        for term, system in [('digest', 'digestive'), ('respirat', 'respiratory'), ('urin', 'urinary'), ('genital', 'reproductive'), ('reproduct', 'reproductive'), ('endocrin', 'endocrine')]:
            if term in group_text:
                return system
        if any(word in name for word in ('kidney', 'ureter', 'bladder')):
            return 'urinary'
        if any(word in name for word in ('lung', 'trachea', 'bronch')):
            return 'respiratory'
        return 'digestive'
    return 'lymphatic' if obj.name == 'Lymph node' else None


def license_for(name):
    if re.search(r'kidney(?:\.\d+)?\.[lr]$', name, re.I):
        return 'CC BY-NC 4.0'
    if re.search(r'(cochlea|vestibule)(?:\.\d+)?\.[lr]$', name, re.I):
        return 'CC BY-NC-SA 4.0'
    return 'CC BY-SA 4.0'


def display_name(name):
    name = name.strip()
    return re.sub(r'\.([lr])$', lambda m: ' (left)' if m.group(1) == 'l' else ' (right)', name)


parts = []
concepts = []
materials = {}
collections = defaultdict(list)
chunks = []
blob = bytearray()
triangles = 0
skipped = []


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
    filename = f'{ARGS.name}-{len(chunks)}.bin'
    raw = bytes(blob)
    (ARGS.out / filename).write_bytes(raw)
    compressed = gzip.compress(raw, compresslevel=9, mtime=0)
    (ARGS.out / (filename + '.gz')).write_bytes(compressed)
    chunks.append({'url': '/models/' + filename, 'bytes': len(raw), 'gzip': '/models/' + filename + '.gz', 'gzipBytes': len(compressed)})
    blob = bytearray()


# Vessels, peripheral nerves, ducts, and ureters in the original file are
# mostly curves. Convert each evaluated curve to its rendered tube mesh while
# retaining its original object name, material, collections, and transform.
depsgraph = bpy.context.evaluated_depsgraph_get()
for obj in sorted((o for o in bpy.data.objects if o.type in ('MESH', 'CURVE')), key=lambda o: o.name):
    if obj.name.endswith('.g') or any(slot.material and slot.material.name in ('Text', 'Lines', 'Planes', 'Directions', 'Movement') for slot in obj.material_slots):
        skipped.append({'name': obj.name, 'reason': 'Reference lettering or diagram'})
        continue
    evaluated = obj.evaluated_get(depsgraph) if obj.type == 'CURVE' else None
    source_mesh = evaluated.to_mesh() if evaluated else obj.data
    if hasattr(source_mesh, 'calc_loop_triangles'):
        source_mesh.calc_loop_triangles()
    if not source_mesh.loop_triangles:
        if evaluated:
            evaluated.to_mesh_clear()
        skipped.append({'name': obj.name, 'reason': 'No triangles'})
        continue
    triangles_by_material = defaultdict(list)
    for triangle in source_mesh.loop_triangles:
        triangles_by_material[source_mesh.polygons[triangle.polygon_index].material_index].append(triangle)
    node = evaluated.matrix_world if evaluated else obj.matrix_world
    normal_matrix = node.to_3x3().inverted_safe().transposed()
    reflected = node.determinant() < 0
    object_parts = []
    for slot_index, group in sorted(triangles_by_material.items()):
        material = obj.material_slots[slot_index].material if slot_index < len(obj.material_slots) else None
        material_name = material.name if material else 'Unassigned'
        system = system_for(obj, material_name)
        if system is None:
            skipped.append({'name': obj.name, 'reason': 'Non-anatomy reference object'})
            break
        lookup = {}
        positions = []
        normals = []
        indices = []
        minimum = [math.inf] * 3
        maximum = [-math.inf] * 3
        for triangle in group:
            face = list(triangle.vertices)
            if reflected:
                face[1], face[2] = face[2], face[1]
            for source_index in face:
                target = lookup.get(source_index)
                if target is None:
                    vertex = source_mesh.vertices[source_index]
                    p = position(node @ vertex.co)
                    n = position((normal_matrix @ vertex.normal).normalized())
                    target = len(lookup)
                    lookup[source_index] = target
                    positions.extend(p)
                    normals.extend(max(-32767, min(32767, round(v * 32767))) for v in n)
                    for axis in range(3):
                        minimum[axis] = min(minimum[axis], p[axis])
                        maximum[axis] = max(maximum[axis], p[axis])
                indices.append(target)
        if len(blob) > 6_000_000:
            flush()
        id = 'ZA:' + obj.name + (':' + material_name if len(triangles_by_material) > 1 else '')
        assert not any(p['id'] == id for p in parts), id
        material_id = 'ZA:' + system + ':' + material_name
        original_color = source_color(material)
        materials[material_id] = {'color': display_color(system, material_name, original_color), 'sourceColor': original_color}
        license = license_for(obj.name)
        part = {
            'id': id, 'conceptId': 'ZA:' + obj.name, 'name': display_name(obj.name) + (f' · {material_name}' if len(triangles_by_material) > 1 else ''),
            'system': system, 'chunk': len(chunks), 'positions': append('f', positions), 'normals': append('h', normals),
            'indices': append('I', indices), 'vertexCount': len(lookup), 'indexCount': len(indices),
            'bounds': [minimum, maximum], 'sourceId': obj.name, 'material': material_id,
            'groups': sorted(c.name for c in obj.users_collection), 'license': license,
            'provenance': {'label': 'Z-Anatomy original Blender model', 'url': 'https://github.com/Z-Anatomy/Models-of-human-anatomy', 'detail': license + '; original object and material, world transform baked.'},
        }
        parts.append(part)
        object_parts.append(id)
        triangles += len(indices) // 3
    if object_parts:
        concepts.append({'id': 'ZA:' + obj.name, 'name': display_name(obj.name), 'elements': object_parts})
        for collection in obj.users_collection:
            collections[collection.name].extend(object_parts)
    if evaluated:
        evaluated.to_mesh_clear()

flush()
for name, members in sorted(collections.items()):
    if len(members) > 1:
        concepts.append({'id': 'ZA:collection:' + name, 'name': name, 'elements': members})
atlas = {
    'version': 'Z-Anatomy original Blender source', 'quality': 'full', 'sex': 'male',
    'source': 'Z-Anatomy', 'scope': 'Direct Z-Anatomy source geometry; named material surfaces and native display colors.',
    'parts': parts, 'concepts': concepts, 'chunks': chunks, 'materials': materials, 'triangles': triangles,
    'provenance': {
        'sourceUrl': 'https://github.com/Z-Anatomy/Models-of-human-anatomy',
        'archiveSha256': ARGS.source_sha256,
        'sourceFile': 'Z-Anatomy/Startup.blend', 'sourceFileSha256': hashlib.sha256(Path(bpy.data.filepath).read_bytes()).hexdigest(),
        'excludedReferenceObjects': skipped,
        'notes': 'The original archive has material colors but no anatomical image texture files.'
    },
}
(ARGS.out / (ARGS.name + '.json')).write_text(json.dumps(atlas, separators=(',', ':')))
print('HUMAN_ATLAS_BUILD=' + json.dumps({'parts': len(parts), 'concepts': len(concepts), 'triangles': triangles, 'chunks': len(chunks), 'materials': len(materials), 'skipped': len(skipped)}))
