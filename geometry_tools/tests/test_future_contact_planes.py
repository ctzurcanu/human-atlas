import unittest
import numpy as np
from geometry_tools.algorithms.future_contact_planes import future_contact_rows


class FutureContactPlanesTests(unittest.TestCase):
    def test_complete_plane_protects_corners_not_moving_in_old_probe(self):
        v=np.array([[0,-.001,0],[0,.001,0],[0,0,.001],
                    [.001,-.001,0],[.001,.001,0],[.001,0,.001]])
        f=np.array([[0,1,2],[3,4,5]]);owner=np.arange(6)
        probe=np.zeros((6,3));probe[3,0]=-.002
        fixed=np.array([True,True,True,False,False,False])
        feedback=dict(vertices=v,faces=f,owner=owner,pairs=np.array([[0,1]]),probe_shift=probe)
        partial,_,_,_=future_contact_rows(v,f,owner,feedback,fixed,gap_m=.00005)
        self.assertEqual(set(partial[:,1]),{3})
        groups,normals,rhs,proof=future_contact_rows(v,f,owner,feedback,fixed,gap_m=.00005,complete_pair_planes=True)
        self.assertEqual(set(groups[:,1]),{3,4,5})
        self.assertTrue(proof['completePairPlanes'])
        # A new objective may now move every corner, including those that
        # were safe in the rejected probe. All must stay on the same side.
        direction=np.zeros((6,3));direction[3:,0]=-.002
        from geometry_tools.algorithms.clearance_direction_projection import project_cell_direction
        for backend in ['osqp','clarabel']:
            with self.subTest(backend=backend):
                q,projection=project_cell_direction(v,f,owner,fixed,direction,groups,normals,rhs,solver_backend=backend)
                self.assertTrue(np.all(q[3:,0]>=-.00095-2e-9))
                np.testing.assert_array_equal(q[:3],np.zeros((3,3)))
                self.assertLessEqual(projection['maximumContactPlaneViolation_m'],2e-9)

    def test_available_gap_is_not_discarded(self):
        v=np.array([[0,-.001,0],[0,.001,0],[0,0,.001],
                    [.001,-.001,0],[.001,.001,0],[.001,0,.001]])
        f=np.array([[0,1,2],[3,4,5]]);owner=np.array([0,0,0,1,1,1]);shift=np.array([[0,0,0],[-.002,0,0]])
        feedback=dict(vertices=v,faces=f,owner=owner,pairs=np.array([[0,1]]),probe_shift=shift)
        groups,normals,rhs,p=future_contact_rows(v,f,owner,feedback,[True,False],gap_m=.00005)
        self.assertGreater(len(rhs),0);self.assertTrue(np.all(rhs<0))
        np.testing.assert_allclose(normals,[np.array([1.,0,0])]*len(rhs))
        np.testing.assert_allclose(rhs,-.00095,atol=1e-14)
        from geometry_tools.algorithms.clearance_direction_projection import project_cell_direction
        q,proof=project_cell_direction(v,f,owner,np.array([True,False]),shift,groups,normals,rhs)
        np.testing.assert_array_equal(q[0],shift[0])
        self.assertGreaterEqual(q[1,0],-.00095-2e-9)
        self.assertLessEqual(proof['maximumContactPlaneViolation_m'],2e-9)
        self.assertGreater(proof['linearConstraints']['removedExactDuplicateRows'],0)
        self.assertEqual(proof['sourceFaceOrientationRowsExcluded'],2)
        with self.assertRaises(ValueError):future_contact_rows(v+.0001,f,owner,feedback,[True,False],gap_m=.00005)


if __name__=='__main__':unittest.main()
