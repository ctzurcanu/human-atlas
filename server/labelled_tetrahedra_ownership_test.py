import unittest
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.algorithms.labelled_tetrahedra import partition_labelled_tetrahedra


class IdenticalAffineOwnershipTests(unittest.TestCase):
    def setUp(self):
        self.positions=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]])
        self.tetrahedra=np.array([[0,1,2,3],[0,2,1,4]])

    @staticmethod
    def volume(positions,faces):
        q=positions[faces]
        return np.einsum('ij,ij->i',q[:,0],np.cross(q[:,1],q[:,2])).sum()/6

    def test_positive_volume_ambiguity_requires_caller_policy(self):
        with self.assertRaisesRegex(ValueError,'duplicate positive volume'):
            partition_labelled_tetrahedra(self.positions,self.tetrahedra[:1],np.ones((5,2)))

    def test_identical_cell_has_one_owner_and_unchanged_controls(self):
        result=partition_labelled_tetrahedra(self.positions,self.tetrahedra[:1],np.ones((5,2)),
                                           identical_score_priority=[1,0])
        self.assertEqual(len(result['cells'][0]),0)
        self.assertAlmostEqual(self.volume(result['positions'],result['cells'][1]),1/6)
        np.testing.assert_array_equal(result['positions'],self.positions)
        self.assertEqual(result['evidence']['identical_affine_score_ownership'],
                         [dict(tetrahedron=0,identical_labels=[0,1],owner=1)])

    def test_tied_cell_and_neighbor_keep_one_opposed_interface(self):
        scores=np.ones((5,2));scores[4]=[2.,0.]
        result=partition_labelled_tetrahedra(self.positions,self.tetrahedra,scores,
                                           identical_score_priority=[1,0])
        first,second=result['cells'][0],result['cells'][1]
        shared=set(map(tuple,np.sort(first,axis=1)))&set(map(tuple,np.sort(second,axis=1)))
        self.assertEqual(shared,{(0,1,2)})
        a=next(f for f in first if set(f)=={0,1,2});b=next(f for f in second if set(f)=={0,1,2})
        self.assertTrue(any(np.array_equal(np.roll(a,i),b[[0,2,1]])for i in range(3)))
        for faces in [first,second]:self.assertAlmostEqual(self.volume(result['positions'],faces),1/6)

    def test_shared_face_tie_without_volume_tie_needs_no_priority(self):
        scores=np.ones((5,2));scores[3]=[0.,2.];scores[4]=[2.,0.]
        result=partition_labelled_tetrahedra(self.positions,self.tetrahedra,scores)
        self.assertEqual(result['evidence']['identical_affine_score_ownership'],[])
        self.assertAlmostEqual(sum(self.volume(result['positions'],f)for f in result['cells'].values()),1/3)

    def test_invalid_priority_rejected(self):
        for order in [[0,0],[0],[0,1,2],[False,1],[0.,1.]]:
            with self.subTest(order=order),self.assertRaisesRegex(ValueError,'every integer label'):
                partition_labelled_tetrahedra(self.positions,self.tetrahedra,np.ones((5,2)),
                                             identical_score_priority=order)

    def test_iterable_policy_is_retained_for_replay(self):
        result=partition_labelled_tetrahedra(self.positions,self.tetrahedra[:1],np.ones((5,2)),
                                           identical_score_priority=iter([1,0]))
        self.assertEqual(result['evidence']['identical_affine_score_priority'],[1,0])


if __name__=='__main__':
    unittest.main()
