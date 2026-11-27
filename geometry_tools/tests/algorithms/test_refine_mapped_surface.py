import unittest
from collections import Counter
import numpy as np
from geometry_tools.algorithms.refine_mapped_surface import refine_mapped_surface


class MappedRefinementTests(unittest.TestCase):
    def test_shared_edge_refinement_preserves_rims_source_facets_and_original_controls(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])
        f=np.array([[0,1,2],[0,2,3]])
        def field(q):return q+np.column_stack([np.zeros(len(q)),np.zeros(len(q)),q[:,0]*q[:,1]])
        old=field(p);r=refine_mapped_surface(p,old,f,[(0,2)],field)
        np.testing.assert_array_equal(r['positions'][:4],old)
        np.testing.assert_array_equal(r['source_positions'][:4],p)
        self.assertEqual(r['shared_edge_controls'],1);self.assertEqual(r['new_controls'],3)
        def rims(f):
            c=Counter(tuple(sorted([int(a),int(b)]))for t in f for a,b in zip(t,np.roll(t,-1)));return sorted(e for e,n in c.items()if n==1)
        self.assertEqual(rims(r['faces']),rims(f))
        source_tri=r['source_positions'][r['faces']]
        cross=np.cross(source_tri[:,1]-source_tri[:,0],source_tri[:,2]-source_tri[:,0])
        self.assertTrue(np.all(cross[:,2]>0));self.assertAlmostEqual(np.linalg.norm(cross,axis=1).sum()/2,1.)
        np.testing.assert_allclose(r['positions'][4:],field(r['source_positions'][4:]),atol=0,rtol=0)
        self.assertEqual(set(r['parent_faces']),{0,1})

    def test_protected_rims_and_invalid_edges_fail_before_mapping(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        f=np.array([[0,1,2]])
        for edges in [[(0,1)],[(0,9)]]:
            with self.assertRaises(ValueError):refine_mapped_surface(p,p,f,edges,lambda q:q)


if __name__=='__main__':unittest.main()
