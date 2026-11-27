"""Compact runtime replay tests; no atlas outputs or source catalogues change."""
import importlib.util
import pathlib
import struct
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('runtime_attributes', ROOT / 'scripts/ta98_runtime_attributes.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class RuntimeAttributes(unittest.TestCase):
    def test_compact_shared_tissue_repack_retains_chart_and_normal_words(self):
        # Deliberately give the source fields nonzero offsets and sentinel bytes.
        fields = {'positions': struct.pack('<12f', *range(12)),
                  'normals': struct.pack('<12h', *range(-6, 6)),
                  'tissueChart': struct.pack('<12h', 0, 32767, 0, 500, 0, 32767,
                                            1000, -32767, 0, 1500, 0, -32767),
                  'normalSigns': bytes([1, 255, 1, 255]),
                  'indices': struct.pack('<6I', 0, 1, 2, 0, 2, 3)}
        source = bytearray(b'prefix-marker')
        part = {'id': 'thyroid-cell', 'vertexCount': 4, 'indexCount': 6}
        for key, values in fields.items():
            part[key] = len(source)
            source.extend(values)
            source.extend(b'padding')
        original = bytes(source)
        saved = module.packed_attributes(part, original)
        packed = bytearray()
        target = dict(part)
        for key, values in saved.items():
            while len(packed) % 4:
                packed.append(0)
            target[key] = len(packed)
            packed.extend(values)
        self.assertEqual(module.packed_attributes(target, bytes(packed)), fields)
        self.assertEqual(original, bytes(source))
        self.assertNotEqual(part['tissueChart'], target['tissueChart'])

    def test_assembly_marker_is_preserved_without_aliasing(self):
        marker = {'schema': 'ta98-thyroid-wall-v1', 'capsuleParts': ['capsule'],
                  'parenchymaParts': ['parenchyma'], 'sourceReview': {'status': 'provisional'}}
        atlas = {'version': 'fixture', 'thyroidWallAssembly': marker, 'parts': []}
        metadata = module.runtime_metadata(atlas)
        self.assertEqual(metadata['thyroidWallAssembly'], marker)
        metadata['thyroidWallAssembly']['capsuleParts'].append('different')
        self.assertEqual(marker['capsuleParts'], ['capsule'])
        self.assertNotIn('parts', metadata)

    def test_truncated_or_invalid_ranges_fail_before_output(self):
        for offset in [-1, 2.5, 1000]:
            with self.assertRaises(ValueError):
                module.packed_attributes({'id': 'bad', 'vertexCount': 3, 'indexCount': 0,
                                          'tissueChart': offset}, bytes(18))


if __name__ == '__main__':
    unittest.main()
