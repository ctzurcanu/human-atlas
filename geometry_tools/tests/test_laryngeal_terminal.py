"""Independent checks of the actual packed bilateral laryngeal repair."""
import json
from pathlib import Path
import unittest

import numpy as np

from geometry_tools.model import AtlasReader
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.oriented_rays import boundary_loops

BASE = Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')
PACKED = BASE / 'blood-installed-v1889/circulation-atlas.json'


@unittest.skipUnless(PACKED.is_file(), 'Local unredistributed source fixture')
class LaryngealTerminalRepair(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.new = json.loads(PACKED.read_text())
        cls.old_path = BASE / 'blood-installed-v1876/circulation-atlas.json'
        cls.old = json.loads(cls.old_path.read_text())
        cls.rows = json.loads((BASE / 'laryngeal-terminal-rotation-v1888/report.json').read_text())['rows']
        cls.reader = AtlasReader()
        cls.meshes = {}
        for row in cls.rows:
            pid = row['partId']
            cls.meshes[pid] = [cls.reader.read_part(a, next(p for p in a['parts'] if p['id'] == pid),
                                                   frame='model:local-ta98:assembled-metres', catalogue_path=path)
                               for a, path in [(cls.old, cls.old_path), (cls.new, PACKED)]]

    def test_skeleton_and_all_other_owners_and_terms_unchanged(self):
        changed = {r['partId'] for r in self.rows}
        allowed = {'positions', 'normals', 'normalSigns', 'uvs', 'colors', 'indices', 'chunk', 'bounds'}
        self.assertEqual([p['id'] for p in self.old['parts']], [p['id'] for p in self.new['parts']])
        for p, q in zip(self.old['parts'], self.new['parts']):
            if p['id'] not in changed:
                self.assertEqual(p, q)
            else:
                self.assertEqual({k:v for k,v in p.items() if k not in allowed},
                                 {k:v for k,v in q.items() if k not in allowed})
        self.assertEqual(self.old['concepts'], self.new['concepts'])
        self.assertEqual(self.old['chunks'], self.new['chunks'][:-1])

    def test_only_six_tip_controls_change_without_calibre_loss(self):
        for row in self.rows:
            old, new = self.meshes[row['partId']]; tip = np.asarray(row['tipRing'])
            fixed = np.ones(len(old.vertices), bool); fixed[tip] = False
            np.testing.assert_array_equal(old.vertices[fixed], new.vertices[fixed])
            np.testing.assert_array_equal(old.faces, new.faces)
            self.assertLess(np.linalg.norm(new.vertices - old.vertices, axis=1).max(), .000061)
            np.testing.assert_allclose(old.vertices[tip].mean(0), new.vertices[tip].mean(0), atol=2e-9, rtol=0)
            np.testing.assert_allclose(np.linalg.norm(old.vertices[tip,None]-old.vertices[None,tip], axis=2),
                                       np.linalg.norm(new.vertices[tip,None]-new.vertices[None,tip], axis=2),
                                       atol=4e-9, rtol=0)

    def test_actual_packed_surface_no_crossings_caps_or_extra_tunnels(self):
        for old, new in self.meshes.values():
            self.assertEqual(diagnostics(old, include_pairs=True)['self_intersection_pairs'], 29)
            diag = diagnostics(new, include_pairs=True)
            self.assertEqual(diag['self_intersection_pairs'], 0)
            self.assertEqual(diag['nonmanifold_edges'], 0)
            edges, loops, _, boundary = boundary_loops(new)
            self.assertEqual(len(loops), 2); self.assertEqual(len(edges), 12)
            self.assertEqual(boundary['inconsistent_shared_edge_winding'], 0)
            all_edges = np.unique(np.sort(np.concatenate([new.faces[:,[0,1]], new.faces[:,[1,2]],
                                                          new.faces[:,[2,0]]]),axis=1),axis=0)
            self.assertEqual(len(new.vertices)-len(all_edges)+len(new.faces), 0)

    def test_shading_repaired_locally_and_source_appearance_data_retained(self):
        for row in self.rows:
            old, new = self.meshes[row['partId']]; touched = row['normalTransport']['changedStars']
            fixed = np.ones(len(old.vertices), bool); fixed[touched] = False
            for k in old.attributes:
                if k != 'normals':
                    np.testing.assert_array_equal(old.attributes[k], new.attributes[k])
            np.testing.assert_array_equal(old.attributes['normals'][fixed], new.attributes['normals'][fixed])
            t = new.triangles; facet = np.cross(t[:,1]-t[:,0], t[:,2]-t[:,0])
            effective = new.attributes['normals'].astype(float)*new.attributes['normalSigns'][:,None]
            self.assertTrue(np.all(np.einsum('ij,ij->i', facet, effective[new.faces].mean(1)) > 0))
            for vertex in touched:
                self.assertTrue(np.all(facet[np.any(new.faces==vertex,axis=1)]@effective[vertex] > 0))

    def test_all_previously_protected_source_mouths_stay_verified(self):
        p = BASE / 'blood-installed-v1889/source-mouth-replay.json'
        self.assertTrue(p.is_file(), 'Packed finalization must replay source interfaces')
        d = json.loads(p.read_text())
        self.assertEqual(d['summary']['verified'], 260)
        self.assertEqual(d['summary']['pending'], 0)
        for interface in d['interfaces']:
            self.assertTrue(interface['sourcePortProof']['passed'])


if __name__ == '__main__':
    unittest.main()
