import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.monotone_uncross import monotone_uncross


class MonotoneSelfTests(unittest.TestCase):
    def test_core_uncrossing_preserves_other_source_owner(self):
        v=np.array([[-.001,-.001,0.],[.001,-.001,0.],[0.,.001,0.],
                    [-.0005,0.,-.0000002],[.0005,0.,.000002],[0.,-.0005,.000002]])
        v+=np.array([.3,.8,.07]);ref=v.copy();ref[3,2]=.070002
        f=np.array([[0,1,2],[3,4,5]])
        m=Mesh(v,f,'test:metres');r=Mesh(ref,f,'test:metres')
        q,proof=monotone_uncross(m,r,np.arange(6),[0,1,2,4,5],face_owners=[0,1],maximum_accepted_steps=2)
        self.assertTrue(proof['complete'],proof)
        self.assertEqual(proof['initialSelfPairs'],1)
        self.assertEqual(proof['finalSelfPairs'],0)
        np.testing.assert_array_equal(q[[0,1,2,4,5]],v[[0,1,2,4,5]])
        self.assertLessEqual(proof['maximumSourceMotion_m'],.00001)

    def test_deferred_pair_updates_equal_every_step_full_replay(self):
        v=np.array([[-.001,-.001,0.],[.001,-.001,0.],[0.,.001,0.],
                    [-.0005,0.,-.0000002],[.0005,0.,.000002],[0.,-.0005,.000002]])
        v=np.r_[v,v+[.01,0.,0.]]+np.array([.3,.8,.07]);ref=v.copy();ref[[3,9],2]=.070002
        f=np.array([[0,1,2],[3,4,5],[6,7,8],[9,10,11]])
        m=Mesh(v,f,'test:metres');r=Mesh(ref,f,'test:metres')
        protected=[i for i in range(12)if i not in[3,9]]
        kwargs=dict(face_owners=[0,1,2,3],maximum_accepted_steps=2)
        a,pa=monotone_uncross(m,r,np.arange(12),protected,**kwargs)
        b,pb=monotone_uncross(m,r,np.arange(12),protected,defer_full_replay=True,**kwargs)
        np.testing.assert_array_equal(a,b)
        self.assertTrue(pa['complete'] and pb['complete'])
        self.assertEqual(pb['initialSelfPairs'],2)
        self.assertEqual(pb['finalSelfPairs'],0)

    def test_opposed_source_sections_keep_both_shapes(self):
        v=np.array([[-.001,-.001,0.],[.001,-.001,0.],[0.,.001,0.],
                    [-.0005,0.,-.0000002],[.0005,0.,.000002],[0.,-.0005,.000002]])
        ref=v.copy();ref[3,2]=.000002
        f=np.array([[0,1,2],[3,4,5]])
        q,proof=monotone_uncross(Mesh(v,f,'test:metres'),Mesh(ref,f,'test:metres'),
              np.arange(6),[],rigid_groups=[[0,1,2],[3,4,5]],coupled_groups=True,
              face_owners=[0,1],defer_full_replay=True,maximum_accepted_steps=2)
        self.assertTrue(proof['complete'],proof)
        self.assertTrue(any(s['coupledSourceSections']for s in proof['acceptedSteps']),proof)
        for group in [[0,1,2],[3,4,5]]:
            delta=q[group]-v[group]
            np.testing.assert_allclose(delta,np.broadcast_to(delta[0],delta.shape),atol=1e-15)
        np.testing.assert_allclose((q-v).mean(0),0.,atol=1e-15)

    def test_fixed_interface_locks_its_entire_rigid_section(self):
        v=np.array([[-.001,-.001,0.],[.001,-.001,0.],[0.,.001,0.],
                    [-.0005,0.,-.0000002],[.0005,0.,.000002],[0.,-.0005,.000002]])
        ref=v.copy();ref[3,2]=.000002
        f=np.array([[0,1,2],[3,4,5]])
        q,proof=monotone_uncross(Mesh(v,f,'test:metres'),Mesh(ref,f,'test:metres'),
              np.arange(6),[0],rigid_groups=[[0,1,2],[3,4,5]],
              face_owners=[0,1],maximum_accepted_steps=2)
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_array_equal(q[:3],v[:3])
        delta=q[3:]-v[3:]
        np.testing.assert_allclose(delta,np.broadcast_to(delta[0],delta.shape),atol=1e-15)


if __name__=='__main__':unittest.main()
