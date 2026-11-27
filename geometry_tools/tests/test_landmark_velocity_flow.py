import unittest
import numpy as np
from geometry_tools.algorithms.landmark_velocity_flow import prepare_landmark_velocity,advance_landmark_flow


class LandmarkFlowTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0,0,0],[.01,0,0],[0,.01,0],[0,0,.01]],float)

    def test_identity_is_literal_and_fixed_roots_stay_literal(self):
        bank,_=prepare_landmark_velocity(self.p,self.p)
        points=np.r_[self.p,[[0,1.2,0]]]
        q,proof=advance_landmark_flow(points,bank)
        self.assertTrue(proof['complete'],proof);np.testing.assert_array_equal(q,points)

    def test_landmarks_follow_the_same_field_and_a_saved_stage_resumes(self):
        target=self.p+[.003,.002,0]
        bank,_=prepare_landmark_velocity(self.p,target,kernel_width_m=[.018,.006],kernel_weights=[1.,.1],regularization=1e-9)
        points=np.r_[self.p,[[0,1.2,0]]]
        full,proof=advance_landmark_flow(points,bank)
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_allclose(full[:4],target,atol=1e-7);np.testing.assert_array_equal(full[-1],points[-1])
        part,first=advance_landmark_flow(points,bank,maximum_steps=1)
        self.assertFalse(first['complete']);self.assertTrue(first['maximumStepsExhausted'])
        resumed,last=advance_landmark_flow(part,bank,start_time=first['endTime'])
        self.assertTrue(last['complete'],last);np.testing.assert_allclose(resumed,full,atol=1e-7)

    def test_cubic_transport_matches_landmarks_and_separate_point_batches(self):
        target=self.p*np.array([1.12,.94,1.04])+[.002,.001,-.001]
        bank,_=prepare_landmark_velocity(self.p,target,regularization=1e-9)
        points=np.r_[self.p,[[.003,.003,.003],[0,1.2,0]]]
        full,proof=advance_landmark_flow(points,bank,temporal_interpolation='cubic')
        self.assertTrue(proof['complete'],proof)
        self.assertEqual(proof['temporalCoefficientInterpolation'],'cubic')
        np.testing.assert_allclose(full[:4],target,atol=1e-7)
        pieces=[]
        for p in [points[:3],points[3:]]:
            q,check=advance_landmark_flow(p,bank,temporal_interpolation='cubic')
            self.assertTrue(check['complete'],check);pieces.append(q)
        np.testing.assert_allclose(np.concatenate(pieces),full,atol=1e-7)
        np.testing.assert_array_equal(full[-1],points[-1])
        with self.assertRaises(ValueError):
            advance_landmark_flow(points,bank,temporal_interpolation='unspecified')

    def test_declared_stage_finishes_without_claiming_whole_flow(self):
        target=self.p*np.array([1.1,.96,1.02])+[.002,.001,0]
        bank,_=prepare_landmark_velocity(self.p,target,regularization=1e-9)
        points=np.r_[self.p,[[0,1.2,0]]]
        full,check=advance_landmark_flow(points,bank,temporal_interpolation='cubic')
        self.assertTrue(check['wholeFlowComplete'])
        stage,first=advance_landmark_flow(points,bank,end_time=.375,temporal_interpolation='cubic')
        self.assertTrue(first['complete']);self.assertFalse(first['timedOut'])
        self.assertFalse(first['wholeFlowComplete']);self.assertEqual(first['endTime'],.375)
        resumed,last=advance_landmark_flow(stage,bank,start_time=.375,temporal_interpolation='cubic')
        self.assertTrue(last['wholeFlowComplete']);np.testing.assert_allclose(resumed,full,atol=1e-7)
        np.testing.assert_array_equal(resumed[-1],points[-1])
        with self.assertRaises(ValueError):advance_landmark_flow(points,bank,start_time=.5,end_time=.25)


if __name__=='__main__':unittest.main()
