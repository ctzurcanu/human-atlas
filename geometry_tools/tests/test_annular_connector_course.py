import unittest
import numpy as np
from geometry_tools.commands.pose_annular_connector import guided_course


class ConnectorCourseTests(unittest.TestCase):
    def test_fixed_mouths_and_tangents_are_literal_with_round_interior(self):
        guide=np.array([[0.,0.,0.],[.003,.01,0.],[.004,.02,.002],[.005,.03,.003]])
        start=np.array([0.,-.005,0.]);end=np.array([.005,.03,.004])
        first=np.array([0.,1.,0.]);last=np.array([0.,0.,1.])
        path,tangent=guided_course(guide,start,end,first,last,np.linspace(0.,1.,20))
        np.testing.assert_array_equal(path[[0,-1]],[start,end])
        np.testing.assert_array_equal(tangent[[0,-1]],[first,last])
        self.assertTrue(np.isfinite(path).all())
        self.assertTrue((np.linalg.norm(np.diff(path,axis=0),axis=1)>0).all())
        np.testing.assert_allclose(np.linalg.norm(tangent,axis=1),1.,atol=1e-14)


if __name__=='__main__':unittest.main()
