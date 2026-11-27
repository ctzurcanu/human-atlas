import unittest
import numpy as np
from geometry_tools.algorithms.shared_edge_collapse import collapse_short_edges
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard


def octahedron():
    v = np.array([[0., 0., 1.], [0., 0., -1.], [1e-9, 0., 1.],
                  [0., 1., 0.], [-1., 0., 0.], [0., -1., 0.]])
    f = np.array([[0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 2],
                  [1, 3, 2], [1, 4, 3], [1, 5, 4], [1, 2, 5]])
    p = v[f]
    n = np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None]
    w = np.repeat(n[:, None], 3, axis=1)
    w[[0, 3]] *= -1
    return v, f, w


class SharedEdgeCollapseTests(unittest.TestCase):
    def test_material_link_above_uint64_capacity_is_still_checked(self):
        v,f,w=octahedron();vertices=[];surfaces={};cells={};words={}
        for i in range(65):
            sid=str(i);vertices.append(v+[i*4.,0.,0.]);surfaces[sid]=f+i*len(v)
            cells[sid]=[{'surface':sid,'side':1}];words[sid]=w.copy()
            if i<64:words[sid][[0,3]]*=-1
        result=collapse_short_edges(np.concatenate(vertices),surfaces,cells,words,
            maximum_steps=4,wall_seconds=5.)
        self.assertTrue(result['evidence']['passed'])
        self.assertEqual(len(result['evidence']['contractions']),1)
        self.assertEqual(len(result['surfaces']['64']),6)
        np.testing.assert_array_equal(result['surfaces']['0'],surfaces['0'])

    def test_global_guard_preserves_the_closed_surface_during_contraction(self):
        v,f,w=octahedron();guard=LocalCrossingGuard(v,f,[],maximum_displacement_m=5e-8)
        result=collapse_short_edges(v,{'s':f},{'t':[{'surface':'s','side':1}]},{'s':w},global_guard=guard)
        self.assertTrue(result['evidence']['passed']);self.assertEqual(len(guard.remaining),0)
        self.assertEqual(guard.evidence()['removed_faces'],2)
        np.testing.assert_array_equal(guard.faces[guard.active],result['surfaces']['s'])

    def test_locked_source_controls_cannot_be_removed(self):
        v,f,w=octahedron()
        result=collapse_short_edges(v,{'s':f},{'t':[{'surface':'s','side':1}]},{'s':w},locked=[0,2])
        self.assertFalse(result['evidence']['passed']);self.assertEqual(result['evidence']['removed_physical_faces'],0)
        np.testing.assert_array_equal(result['surfaces']['s'],f)

    def test_subnanometre_edge_removal_preserves_closed_boundary(self):
        v, f, w = octahedron()
        source = f.copy()
        result = collapse_short_edges(v, {'skin': f}, {'tissue': [{'surface': 'skin', 'side': 1}]}, {'skin': w})
        self.assertTrue(result['evidence']['passed'])
        self.assertEqual(result['evidence']['removed_physical_faces'], 2)
        out = result['surfaces']['skin']
        edges = np.sort(np.concatenate([out[:, [0, 1]], out[:, [1, 2]], out[:, [2, 0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        self.assertTrue(np.all(counts == 2))
        self.assertEqual(len(np.unique(out))-len(counts)+len(out), 2)
        np.testing.assert_array_equal(f, source)
        self.assertLessEqual(result['evidence']['maximum_accepted_cluster_displacement_m'], 5e-8)

    def test_reciprocal_interface_uses_one_survivor_map(self):
        v, f, w = octahedron()
        cells = {'a': [{'surface': 'shared', 'side': 1}], 'b': [{'surface': 'shared', 'side': -1}]}
        result = collapse_short_edges(v, {'shared': f}, cells, {'shared': w})
        self.assertTrue(result['evidence']['passed'])
        self.assertEqual(len(result['source_faces']['shared']), 6)
        self.assertEqual(len(result['evidence']['contractions']), 1)

    def test_no_global_distance_weld_across_an_open_boundary(self):
        v, f, w = octahedron()
        result = collapse_short_edges(v, {'skin': f[:4]}, {'tissue': [{'surface': 'skin', 'side': 1}]}, {'skin': w[:4]})
        self.assertFalse(result['evidence']['passed'])
        self.assertEqual(result['evidence']['removed_physical_faces'], 0)
        np.testing.assert_array_equal(result['surfaces']['skin'], f[:4])

    def test_displacement_bound_prevents_larger_contraction(self):
        v, f, w = octahedron()
        result = collapse_short_edges(v, {'skin': f}, {'tissue': [{'surface': 'skin', 'side': 1}]}, {'skin': w}, maximum_displacement_m=1e-12)
        self.assertFalse(result['evidence']['passed'])
        self.assertEqual(result['evidence']['removed_physical_faces'], 0)

    def test_one_ring_search_preserves_the_source_bound(self):
        v, f, w = octahedron()
        result = collapse_short_edges(v, {'skin': f}, {'tissue': [{'surface': 'skin', 'side': 1}]}, {'skin': w}, candidate_face_rings=1)
        self.assertTrue(result['evidence']['passed'])
        self.assertLessEqual(result['evidence']['maximum_accepted_cluster_displacement_m'], 5e-8)
        self.assertEqual(result['evidence']['candidate_face_rings'], 1)

    def test_tetrahedral_extinction_rejected(self):
        v = np.array([[0., 0., 0.], [1e-9, 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        f = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
        p = v[f]
        w = np.repeat(np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])[:, None], 3, axis=1)
        w[[0, 1]] *= -1
        r = collapse_short_edges(v, {'s': f}, {'t': [{'surface': 's', 'side': 1}]}, {'s': w})
        self.assertEqual(r['evidence']['removed_physical_faces'], 0)
        self.assertFalse(r['evidence']['passed'])

    def test_declared_pair_is_replayed_instead_of_assumed_to_cross(self):
        v, f, w = octahedron()
        # Both shading failures remain actionable. The common-edge pair is
        # touching, so it must not be recorded as a declared crossing.
        r = collapse_short_edges(v, {'s': f}, {'t': [{'surface': 's', 'side': 1}]}, {'s': w}, defect_face_pairs=[[1, 2]])
        self.assertTrue(r['evidence']['passed'])
        self.assertEqual(r['evidence']['initial_declared_defect_pairs'], 0)

    def test_existing_source_cluster_survives_the_next_stage(self):
        v, f, w = octahedron()
        v=np.r_[v, [[0., 0., 1.+4e-8]]]
        r=collapse_short_edges(v, {'s': f}, {'t': [{'surface': 's', 'side': 1}]}, {'s': w}, initial_clusters={0: [0, 6]})
        self.assertTrue(r['evidence']['passed'])
        cluster=next(q for q in r['evidence']['source_clusters'] if 6 in q['members'])
        self.assertIn(0, cluster['members'])
        self.assertIn(2, cluster['members'])
        self.assertLessEqual(np.linalg.norm(v[cluster['members']]-v[cluster['kept']],axis=1).max(),5e-8)

    def test_explicit_membership_relaxation_still_preserves_both_closed_cells(self):
        v, f, w=octahedron()
        v=np.r_[v, [[2.,0.,1.],[0.,2.,1.],[0.,0.,3.]]]
        other=np.array([[0,7,6],[0,6,8],[0,8,7],[6,7,8]])
        p=v[other];ow=np.repeat(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])[:,None],3,axis=1)
        surfaces={'skin':f,'other':other};cells={'a':[{'surface':'skin','side':1}],'b':[{'surface':'other','side':1}]};words={'skin':w,'other':ow}
        strict=collapse_short_edges(v,surfaces,cells,words)
        self.assertFalse(strict['evidence']['passed'])
        explicit=collapse_short_edges(v,surfaces,cells,words,preserve_material_membership=False)
        self.assertTrue(explicit['evidence']['passed'])
        np.testing.assert_array_equal(explicit['surfaces']['other'],other)
        for faces in explicit['surfaces'].values():
            edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
            self.assertTrue(np.all(np.unique(edges,axis=0,return_counts=True)[1]==2))


if __name__ == '__main__':
    unittest.main()
