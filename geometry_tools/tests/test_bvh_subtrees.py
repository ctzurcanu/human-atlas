import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH


class BVHSubtreeTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(187)
        self.mesh=Mesh(rng.random((96,3)),np.arange(96).reshape(-1,3),'test')
        self.tree=BVH(self.mesh,leaf_size=4)

    def pairs(self,root=None):
        return {tuple(p) for batch in self.tree.candidate_batches(self.tree,0.,selfcheck=True,root_pair=root)for p in batch}

    def test_three_self_subtrees_cover_every_full_candidate_exactly_once(self):
        left,right=map(int,self.tree.children[0])
        groups=[self.pairs((left,left)),self.pairs((left,right)),self.pairs((right,right))]
        self.assertEqual(set.union(*groups),self.pairs())
        for i in range(3):
            for j in range(i):self.assertFalse(groups[i]&groups[j])

    def test_explicit_root_matches_default_query(self):
        self.assertEqual(self.pairs((0,0)),self.pairs())

    def test_invalid_root_cannot_silently_skip_geometry(self):
        with self.assertRaises(ValueError):self.pairs((-1,0))
        with self.assertRaises(ValueError):self.pairs((0,len(self.tree.low)))


if __name__=='__main__':unittest.main()
