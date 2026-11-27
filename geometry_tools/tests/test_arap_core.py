import unittest
import numpy as np
from geometry_tools.algorithms.arap_core import arap_core


class ArapCoreTests(unittest.TestCase):
    def test_known_rigid_motion_has_zero_strain_with_fixed_attachments(self):
        p=np.array([[x,y,0.]for y in range(3)for x in range(3)])
        edges=np.array([(i,j)for i in range(9)for j in range(i+1,9)
                        if np.linalg.norm(p[i]-p[j])<=1.5])
        a=.4
        R=np.array([[np.cos(a),0,np.sin(a)],[0,1,0],[-np.sin(a),0,np.cos(a)]])
        wanted=p@R.T+[.2,-.3,.4]
        fixed=np.array([0,1,2,3,5,6,7,8])
        initial=wanted.copy(); initial[4]+=[.1,.05,.1]
        q,r=arap_core(p,edges,fixed,wanted[fixed],initial,maximum_steps=80)
        np.testing.assert_array_equal(q[fixed],wanted[fixed])
        np.testing.assert_allclose(q,wanted,atol=1e-7)
        self.assertLess(r['energy'][-1],r['energy'][0])

    def test_all_fixed_source_is_literal(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0]])
        q,r=arap_core(p,[[0,1],[1,2],[0,2]],[0,1,2],p,p)
        np.testing.assert_array_equal(q,p)
        self.assertEqual(r['iterations'],0)

    def test_native_deadline_rejects_incomplete_candidate(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0]])
        with self.assertRaises(TimeoutError):
            arap_core(p,[[0,1],[1,2],[0,2]],[0,1],p[:2],p,wall_seconds=1e-12)

    def test_cached_fixed_rotations_preserve_the_complete_solve(self):
        p=np.array([[x*.01,y*.01,.002*np.sin(x*.3+y*.4)]for y in range(9)for x in range(11)])
        edges=np.array([(i,j)for i in range(len(p))for j in range(i+1,len(p))
                        if np.linalg.norm(p[i,:2]-p[j,:2])<.0142])
        free=np.array([4*11+5,4*11+6,5*11+5]);fixed=np.setdiff1d(np.arange(len(p)),free)
        target=p.copy();target[:,2]+=.003*np.exp(-np.sum((p[:,:2]-[.05,.045])**2,axis=1)/.0004)
        initial=target.copy();initial[free,2]+=.001
        cached,proof=arap_core(p,edges,fixed,target[fixed],initial,maximum_steps=80)
        ordinary,reference=arap_core(p,edges,fixed,target[fixed],initial,maximum_steps=80,cache_fixed_rotations=False)
        np.testing.assert_allclose(cached,ordinary,atol=1e-12,rtol=0)
        np.testing.assert_allclose(proof['energy'],reference['energy'],atol=1e-14,rtol=0)
        self.assertGreater(proof['cached_rotation_controls'],len(p)//2)
        self.assertEqual(reference['cached_rotation_controls'],0)


if __name__=='__main__':unittest.main()
