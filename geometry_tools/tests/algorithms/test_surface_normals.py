import unittest
import numpy as np

from geometry_tools.algorithms.surface_normals import grouped_normal_words,source_normal_words
from geometry_tools.algorithms.surface_normals import transport_z_shear_words


class SurfaceNormalTests(unittest.TestCase):
    def test_affine_shear_uses_inverse_transpose_and_keeps_translation_words_literal(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[1,1,0]])
        f=np.array([[0,1,2],[1,3,2]]);w=np.tile([0,0,32767],(4,1)).astype('i2')
        out,proof=transport_z_shear_words(p,f,w,.2*p[:,0]-.3*p[:,1])
        expected=np.array([-.2,.3,1.]);expected=np.rint(expected/np.linalg.norm(expected)*32767).astype('i2')
        np.testing.assert_array_equal(out,np.tile(expected,(4,1)))
        self.assertFalse(proof['positions_modified'])
        raw=np.tile([103,201,32761],(4,1)).astype('i2')
        translated,_=transport_z_shear_words(p,f,raw,np.full(4,.1))
        np.testing.assert_array_equal(translated,raw)

    def test_changed_vertical_chart_requires_an_explicit_different_field(self):
        p=np.array([[0.,0,0],[0,1,0],[0,0,1]])
        f=np.array([[0,1,2]]);w=np.tile([32767,0,0],(3,1)).astype('i2')
        with self.assertRaisesRegex(ValueError,'vertical chart'):
            transport_z_shear_words(p,f,w,np.array([0.,.1,0.]))

    def test_partition_vertices_sample_one_source_field_independently_of_triangulation(self):
        source=np.array([[0.,0,0],[1,0,0],[0,1,0]])
        faces=np.array([[0,1,2]]);words=np.array([[[0,0,32767],[0,10000,31000],[10000,0,31000]]],dtype='i2')
        v=np.r_[source,[[.5,0,0],[.2,.2,0]]];before=v.copy()
        surfaces={'left':np.array([[0,3,4]]),'right':np.array([[3,1,4]])}
        out,proof=source_normal_words(v,surfaces,{'left':np.array([0]),'right':np.array([0])},source,faces,words)
        np.testing.assert_array_equal(v,before)
        np.testing.assert_array_equal(out['left'][0,1],out['right'][0,0])
        np.testing.assert_array_equal(out['left'][0,2],out['right'][0,2])
        self.assertGreater(out['left'][0,1,1],0)
        self.assertEqual(proof['hard_corner_fallbacks'],0)
        self.assertFalse(proof['geometric_positions_modified'])

    def test_a_semantic_chart_split_does_not_leave_a_shading_crease(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,0],[0,1,0],[0,0,1]],float)
        surfaces={'first':np.array([[0,1,2]]),'second':np.array([[3,4,5]])}
        words,_=grouped_normal_words(v,surfaces,{'first':['wall'],'second':['wall']})
        expected=np.array([23170,0,23170],dtype='i2')
        np.testing.assert_array_equal(words['first'][0,0],expected)
        np.testing.assert_array_equal(words['second'][0,0],expected)
        np.testing.assert_array_equal(words['first'][0,2],words['second'][0,1])
        hard,_=grouped_normal_words(v,surfaces,{'first':['outer'],'second':['rim']})
        np.testing.assert_array_equal(hard['first'][0,0],[0,0,32767])
        np.testing.assert_array_equal(hard['second'][0,0],[32767,0,0])


if __name__=='__main__':
    unittest.main()
