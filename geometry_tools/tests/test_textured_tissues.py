"""Keep texture seams when a finite tissue bank shares physical coordinates."""
import unittest
import numpy as np

from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.tissue import normal_words, decode_packed


class TexturedTissueTests(unittest.TestCase):
    def test_one_geometric_vertex_retains_its_distinct_source_uv_charts(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
        f=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
        surfaces={'shell':f};cells={'gland':[dict(surface='shell',side=1,neighbor=None)]}
        uv=np.arange(24,dtype='<f4').reshape(4,3,2)/23
        packed=pack_finite_cells(v,surfaces,normal_words(v,surfaces),cells,{'gland':'gland'},corner_uvs={'shell':uv})
        decoded=decode_packed(packed,'fixture',[0,0,0])['gland']
        np.testing.assert_array_equal(decoded.vertices[decoded.faces],v[f])
        np.testing.assert_array_equal(decoded.attributes['uvs'][decoded.faces],uv)
        same=np.flatnonzero(np.all(decoded.vertices==[0,0,0],axis=1))
        self.assertGreaterEqual(len(same),3)
        self.assertEqual(len(np.unique(decoded.attributes['uvs'][same],axis=0)),3)


if __name__=='__main__':
    unittest.main()
