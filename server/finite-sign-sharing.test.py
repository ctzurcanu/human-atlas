"""Peelable tissues share orientation buffers without changing used normals."""
import unittest
import numpy as np
from geometry_tools.tissue import layer_arrangement,normal_words
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells

class FiniteSignSharing(unittest.TestCase):
    def test_nonadjacent_layers_share_buffers_and_interfaces_stay_opposed(self):
        square=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])
        fields=np.array([square+[0,0,h]for h in [3.,2.,1.,0.]])
        positions,surfaces,cells,_=layer_arrangement(fields,np.array([[0,1,2],[0,2,3]]),['outer','middle','inner'])
        packed=pack_finite_cells(positions,surfaces,normal_words(positions,surfaces),cells,{n:n for n in cells})
        self.assertEqual(packed.cells['outer']['normalSigns'],packed.cells['inner']['normalSigns'])
        self.assertNotEqual(packed.cells['outer']['normalSigns'],packed.cells['middle']['normalSigns'])
        self.assertEqual(len({p['normalSigns']for p in packed.cells.values()}),2)
        for name,uses in cells.items():
            p=packed.cells[name];signs=np.ndarray((p['vertexCount'],),'i1',buffer=packed.data,offset=p['normalSigns'])
            for use in uses:self.assertTrue(np.all(signs[packed.surfaces[use['surface']]]==use['side']))

if __name__=='__main__':unittest.main()
