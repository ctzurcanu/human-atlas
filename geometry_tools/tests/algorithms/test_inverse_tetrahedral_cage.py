import unittest
import numpy as np
from geometry_tools.algorithms.tetrahedral_cage import RegularTetrahedralCage
from geometry_tools.algorithms.inverse_tetrahedral_cage import inverse_cage_points


class InverseCageTests(unittest.TestCase):
    def test_random_controls_round_trip_across_cell_boundaries(self):
        c=RegularTetrahedralCage([0,0,0],.01,[5,5,5]);d=np.zeros_like(c.nodes)
        d[~c.boundary]=np.random.default_rng(4).normal(0,.0008,((~c.boundary).sum(),3))
        self.assertTrue(c.check(d)['passed'])
        p=np.random.default_rng(6).uniform(.0001,.0399,(500,3));target=c.map(p,d)
        r=inverse_cage_points(c,target,d)
        np.testing.assert_allclose(r['source_points'],p,atol=3e-11,rtol=0)
        self.assertLessEqual(r['maximum_residual_m'],2e-11)
    def test_exterior_and_fixed_boundary_are_literal(self):
        c=RegularTetrahedralCage([0,0,0],.01,[4,4,4]);d=np.zeros_like(c.nodes);d[~c.boundary]=[.002,0,0]
        p=np.array([[-1,.02,.02],[.03,.015,.014],[.04,1,-2]])
        np.testing.assert_array_equal(inverse_cage_points(c,p,d)['source_points'],p)
    def test_inverted_field_and_exhausted_iteration_budget_are_rejected(self):
        c=RegularTetrahedralCage([0,0,0],.01,[4,4,4]);d=np.zeros_like(c.nodes)
        d[c.node_ids([1,1,1])]=[.05,0,0]
        with self.assertRaisesRegex(ValueError,'positive'):
            inverse_cage_points(c,[[.015,.015,.015]],d)
        with self.assertRaisesRegex(ValueError,'budgets'):
            inverse_cage_points(c,[[.015,.015,.015]],np.zeros_like(d),maximum_iterations=0)

    def test_reciprocal_obstacle_gradient_uses_the_inverse_transpose(self):
        c=RegularTetrahedralCage([0,0,0],.01,[4,4,4]);d=np.zeros_like(c.nodes)
        d[~c.boundary]=[.0006,-.0003,.0004]
        p=np.array([[.0141,.0178,.0123]]);target=c.map(p,d)
        r=inverse_cage_points(c,target,d);normal=np.array([.4,-.2,.7])
        covector=-np.linalg.solve(r['jacobians'][0].T,normal)
        node=int(r['locations']['node_ids'][0,1]);weight=float(r['locations']['weights'][0,1])
        analytic=weight*covector
        for axis in range(3):
            plus=d.copy();minus=d.copy();plus[node,axis]+=1e-7;minus[node,axis]-=1e-7
            a=inverse_cage_points(c,target,plus)['source_points'][0]@normal
            b=inverse_cage_points(c,target,minus)['source_points'][0]@normal
            self.assertAlmostEqual((a-b)/2e-7,analytic[axis],places=8)


if __name__=='__main__':unittest.main()
