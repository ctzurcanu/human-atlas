import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.model import Mesh
from geometry_tools.algorithms.packed_pose_limits import PackedPoseLimits
from geometry_tools.algorithms.original_facet_pose_limits import OriginalFacetPoseLimits
from geometry_tools.algorithms.original_facet_pose_limits import original_facet_limits

class PackedPoseTests(unittest.TestCase):
    def fixture(self):
        mesh=Mesh(np.array([[0.,0,0],[.01,0,0],[0,.01,0]]),np.array([[0,1,2]]),'test:metres',{}, {})
        groups=np.arange(3);centers=mesh.vertices.copy();levers=np.ones(3)*.01
        model=OriginalFacetPoseLimits(mesh,mesh,groups,centers,levers)
        source={'activeProjections':[{'partId':'test','carrierControlIds':[0,1,2]}]}
        helper=PackedPoseLimits([mesh],[source],{'test':{}},[model],[groups],centers,levers)
        return helper,model,mesh

    def test_fixed_coordinates_survive_literal_and_encoding_is_replayed(self):
        helper,model,mesh=self.fixture();x=np.zeros((6,3))
        exact,packed=helper.vertices(x)
        np.testing.assert_array_equal(packed[0],mesh.vertices)
        x[1,0]=1.23456789e-6
        exact,packed=helper.vertices(x)
        np.testing.assert_array_equal(packed[0][[0,2]],mesh.vertices[[0,2]])
        self.assertEqual(packed[0][1,0],float(np.float32(exact[0][1,0])))
        self.assertNotEqual(packed[0][1,0],exact[0][1,0])

    def test_correction_rows_do_not_mutate_physical_shape_budgets(self):
        helper,model,mesh=self.fixture();x=np.zeros((6,3));x[1,0]=.002
        budget=model.budget.copy();A,b,p=helper.rows(x)
        self.assertGreater(p['violatingExactConstraints'],0)
        self.assertGreater(len(b),0)
        np.testing.assert_array_equal(model.budget,budget)
        self.assertGreater(np.max(b-A@x.ravel()),0)

    def test_total_registration_turn_and_step_turn_are_separate_limits(self):
        _,_,baseline=self.fixture();points=baseline.vertices@Rotation.from_euler('x',15,degrees=True).as_matrix().T
        current=Mesh(points,baseline.faces,baseline.frame,{},{});g=np.arange(3);levers=np.full(3,.01)
        global_model=original_facet_limits(baseline,current,g,points,levers,baseline=baseline)
        step=OriginalFacetPoseLimits(baseline,current,g,points,levers)
        source={'activeProjections':[{'partId':'test','carrierControlIds':[0,1,2]}]}
        helper=PackedPoseLimits([current],[source],{'test':{}},[global_model],[g],points,levers,
            global_turn_degrees=20,step_models=[step],step_turn_degrees=10)
        x=np.zeros((6,3));self.assertEqual(helper.rows(x)[2]['violatingExactConstraints'],0)
        target=points@Rotation.from_euler('x',12,degrees=True).as_matrix().T;x[:3]=target-points
        self.assertGreater(helper.rows(x)[2]['violatingExactConstraints'],0)
        np.testing.assert_array_equal(global_model.budget,np.full(1,.075))

if __name__=='__main__':unittest.main()
