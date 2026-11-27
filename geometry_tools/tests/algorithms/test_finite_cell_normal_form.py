import unittest
import numpy as np
from geometry_tools.algorithms.finite_cell_normal_form import normalize_cell_boundaries


class FiniteCellNormalFormTests(unittest.TestCase):
    def inputs(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        surfaces={'first':np.array([[0,1,2]]),'second':np.array([[0,2,1]]),
                  'film-rest':np.array([[0,2,3]])}
        words={'first':np.tile([0,0,32767],(1,3,1)).astype(np.int16),
               'second':np.tile([0,0,-32767],(1,3,1)).astype(np.int16),
               'film-rest':np.tile([32767,0,0],(1,3,1)).astype(np.int16)}
        uvs={sid:v[f][:,:,:2].astype(np.float32)for sid,f in surfaces.items()}
        cells={'a':[dict(surface='first',side=-1,neighbor='film')],
               'b':[dict(surface='second',side=-1,neighbor='film')],
               'film':[dict(surface='first',side=1,neighbor='a'),
                       dict(surface='second',side=1,neighbor='b'),dict(surface='film-rest',side=1)]}
        return v,surfaces,cells,words,uvs

    def test_exact_collapsed_film_gives_one_reciprocal_material_interface(self):
        v,s,c,w,u=self.inputs()
        result=normalize_cell_boundaries(v,s,c,w,u)
        self.assertEqual(result['evidence']['opposed_cell_boundary_pairs_cancelled'],{'film':1})
        first=result['cells']['a'][0];second=result['cells']['b'][0]
        self.assertEqual(first['surface'],second['surface'])
        self.assertEqual(first['neighbor'],'b');self.assertEqual(second['neighbor'],'a')
        self.assertEqual(first['side'],-second['side'])
        self.assertEqual(len(result['cells']['film']),1)
        self.assertFalse(result['evidence']['coordinates_moved'])
        self.assertTrue(result['evidence']['each_material_oriented_triangle_coefficient_preserved'])
        self.assertEqual(result['evidence']['outward_normal_field_failures'],[])

    def test_repeated_positive_material_area_is_rejected(self):
        v,s,c,w,u=self.inputs()
        c['a'].append(dict(c['a'][0]))
        with self.assertRaisesRegex(ValueError,'Repeated positive-area'):
            normalize_cell_boundaries(v,s,c,w,u)

    def test_opposed_free_faces_cancel_without_removing_live_material(self):
        v,s,_,w,u=self.inputs()
        cells={'material':[dict(surface=sid,side=1)for sid in s]}
        before=v.copy()
        result=normalize_cell_boundaries(v,s,cells,w,u)
        self.assertEqual(result['evidence']['ownerless_physical_faces_cancelled'],1)
        self.assertEqual(result['evidence']['opposed_cell_boundary_pairs_cancelled'],{'material':1})
        self.assertEqual(sum(len(f)for f in result['surfaces'].values()),1)
        use=result['cells']['material'][0]
        np.testing.assert_array_equal(result['surfaces'][use['surface']],s['film-rest'])
        np.testing.assert_array_equal(v,before)
        self.assertTrue(result['evidence']['each_material_oriented_triangle_coefficient_preserved'])

    def test_distinct_parallel_coordinates_are_not_welded(self):
        v,s,c,w,u=self.inputs()
        v=np.r_[v,v[:3]+[0.,0.,1e-10]]
        s['second']=np.array([[4,6,5]])
        result=normalize_cell_boundaries(v,s,c,w,u)
        self.assertEqual(result['evidence']['opposed_cell_boundary_pairs_cancelled'],{})
        self.assertEqual(len(result['cells']['film']),3)


if __name__=='__main__':unittest.main()
