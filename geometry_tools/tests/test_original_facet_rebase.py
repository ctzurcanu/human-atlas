import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.model import Mesh
from geometry_tools.algorithms.original_facet_pose_limits import original_facet_limits

class OriginalFacetRebase(unittest.TestCase):
    def setup_geometry(self,points):
        faces=np.array([[0,1,2]])
        baseline=Mesh(np.array([[0.,0,0],[.01,0,0],[0,.01,0]]),faces,'test:metres',{}, {})
        current=Mesh(points,faces,'test:metres',{}, {})
        return baseline,current,np.arange(3),points.copy(),np.full(3,.01),np.zeros((6,3))

    def test_rebase_cannot_increase_original_strain_budget(self):
        baseline,current,g,c,l,x=self.setup_geometry(np.array([[0.,0,0],[.012,0,0],[0,.01,0]]))
        reset=original_facet_limits(baseline,current,g,c,l)
        protected=original_facet_limits(baseline,current,g,c,l,baseline=baseline)
        self.assertEqual(len(reset.violations(x)[0]),0)
        self.assertEqual(len(protected.violations(x)[0]),1)
        self.assertAlmostEqual(protected.budget[0],.075)

    def test_rebase_cannot_reset_original_facet_turn(self):
        points=np.array([[0.,0,0],[.01,0,0],[0,.01,0]])@Rotation.from_euler('x',15,degrees=True).as_matrix().T
        baseline,current,g,c,l,x=self.setup_geometry(points)
        reset=original_facet_limits(baseline,current,g,c,l)
        protected=original_facet_limits(baseline,current,g,c,l,baseline=baseline)
        self.assertEqual(len(reset.orientation_violations(x)[0]),0)
        self.assertEqual(len(protected.orientation_violations(x)[0]),1)

if __name__=='__main__':unittest.main()
