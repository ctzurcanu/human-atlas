import unittest
import tempfile
from pathlib import Path
import numpy as np

from geometry_tools.bladder_rugae import fold_height, fixed_source_fields, ambient_relief_height


class BladderRugaeFieldTests(unittest.TestCase):
    def test_fixed_trigone_and_neck_are_literal_while_body_has_bounded_relief(self):
        center = np.array([0., .919, .002])
        points = np.c_[np.linspace(-.03, .03, 1001), np.full(1001, .925), np.full(1001, .03)]
        height = fold_height(points, np.full(len(points), .01), center=center)
        self.assertGreater(height.max(), .0003)
        self.assertLessEqual(height.max(), .00035)
        np.testing.assert_array_equal(fold_height(points, np.zeros(len(points)), center=center), np.zeros(len(points)))
        points[:, 1] = .898
        np.testing.assert_array_equal(fold_height(points, np.full(len(points), .01), center=center), np.zeros(len(points)))

    def test_protected_edge_transition_has_no_jump_or_linear_corner(self):
        candidates=np.c_[np.linspace(-.03,.03,1001),np.full(1001,.925),np.full(1001,.03)]
        best=fold_height(candidates,np.ones(len(candidates)),center=[0.,.919,.002]).argmax()
        point = candidates[best:best+1]
        distances = np.array([0., 1e-6, 2e-6, 1e-5])
        p = np.repeat(point, len(distances), axis=0)
        height = fold_height(p, distances, center=[0., .919, .002])
        self.assertEqual(height[0], 0.)
        self.assertGreater(height[1],0.)
        self.assertLess(height[1], 1e-12)
        self.assertLess(height[2], 10*height[1])

    def test_unbounded_or_unphysical_inputs_are_rejected(self):
        p = np.array([[0., .92, .03]])
        for options in [{'height_m': .01}, {'spacing_m': .0001}, {'transition_m': 0.}]:
            with self.assertRaises(ValueError):
                fold_height(p, np.ones(1), center=[0., .919, .002], **options)
        with self.assertRaises(ValueError):
            fold_height(p, -np.ones(1), center=[0., .919, .002])

    def test_static_refinement_restores_native_corner_words(self):
        with tempfile.TemporaryDirectory()as folder:
            path=Path(folder)/'native.npz'
            old=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
            words=np.array([[[100,0,32767],[0,100,32767],[0,0,32767]]],dtype='<i2')
            np.savez(path,encoded_positions=old,faces=np.array([[0,1,2]]),words=words)
            p=np.r_[old,[[.5,0,0]]];surfaces={'lining':np.array([[0,3,2],[3,1,2]])}
            lineage={}
            fields,evidence=fixed_source_fields(p,p,surfaces,np.zeros(len(p),bool),
                dict(source=str(path),frame='fixture',surfaces={'lining':'faces'},corner_words={'lining':'words'}),
                parent_hints={'lining':np.array([0,0])},parent_faces_out=lineage)
            np.testing.assert_array_equal(fields['lining'][0,0],words[0,0])
            np.testing.assert_array_equal(fields['lining'][1,1],words[0,1])
            self.assertEqual(evidence[0]['static_faces'],2)
            np.testing.assert_array_equal(lineage['lining'],[0,0])

    def test_near_centroid_does_not_authorize_outside_source_corners(self):
        with tempfile.TemporaryDirectory()as folder:
            path=Path(folder)/'native.npz';old=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
            np.savez(path,encoded_positions=old,faces=np.array([[0,1,2]]),words=np.tile([0,0,32767],(1,3,1)).astype('<i2'))
            p=np.r_[old,[[-.1,0,0],[1.1,0,0],[.2,.5,0]]]
            with self.assertRaisesRegex(ValueError,'complete native source witness'):
                fixed_source_fields(p,p,{'lining':np.array([[3,4,5]])},np.zeros(len(p),bool),
                    dict(source=str(path),frame='fixture',surfaces={'lining':'faces'},corner_words={'lining':'words'}),
                    parent_hints={'lining':np.array([0])})

    def test_clinical_protection_and_outer_covering_are_independent(self):
        center=[0.,.919,.002]
        points=np.c_[np.linspace(-.03,.03,1001),np.full(1001,.925),np.full(1001,.03)]
        result=ambient_relief_height(points,np.full(len(points),.01),np.full(len(points),.003),center=center)
        self.assertGreater(result.max(),.0003)
        self.assertLessEqual(result.max(),.00035)
        for clinical,covering in [(0.,.003),(.01,0.)]:
            np.testing.assert_array_equal(ambient_relief_height(points,np.full(len(points),clinical),np.full(len(points),covering),center=center),np.zeros(len(points)))
        thin=ambient_relief_height(points,np.full(len(points),.01),np.full(len(points),.0001),center=center)
        self.assertLess(thin.max(),.000003)

    def test_nearly_tangent_average_keeps_outward_direction_after_int16(self):
        from geometry_tools.algorithms.surface_normals import grouped_normal_words
        p=np.array([[0.,0.,0.],[1e-8,0.,0.],[0.,.001,0.],[0.,0.,1.]])
        f=np.array([[0,1,2],[0,2,3]])
        fields,proof=grouped_normal_words(p,{'wall':f},{'wall':['same','same']})
        t=p[f];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        self.assertTrue(np.all(np.einsum('tki,ti->tk',fields['wall'].astype(float),normal)>0))
        self.assertEqual(proof['quantized_hard_corner_fallbacks'],2)


if __name__ == '__main__':
    unittest.main()
