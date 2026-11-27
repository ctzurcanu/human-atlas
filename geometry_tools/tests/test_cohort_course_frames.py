import unittest
import numpy as np
from geometry_tools.algorithms.cohort_course_frames import align_cohort_course_frames


class CourseFrameTests(unittest.TestCase):
    def test_middle_section_faces_its_actual_course_and_root_stays_literal(self):
        theta=np.arange(6)*2*np.pi/6
        ring=np.c_[.001*np.cos(theta),.001*np.sin(theta),np.zeros(6)]
        source=np.concatenate([ring+[0,0,z]for z in [0.,.01,.02]])
        faces=[]
        for k in range(2):
            for i in range(6):
                a=6*k+i;b=6*k+(i+1)%6;c=b+6;d=a+6
                faces.extend([[a,b,c],[a,c,d]])
        sections=[np.arange(6*k,6*k+6)for k in range(3)]
        current=source+np.repeat([[0,0,0],[.005,0,0],[.015,0,0]],6,axis=0)
        q,proof=align_cohort_course_frames(source,current,np.array(faces),sections,source)
        np.testing.assert_array_equal(q[:6],current[:6])
        np.testing.assert_allclose(q[sections[1]].mean(0),current[sections[1]].mean(0),atol=1e-15)
        _,_,vh=np.linalg.svd(q[sections[1]]-q[sections[1]].mean(0))
        tangent=np.array([.015,0,.02]);tangent/=np.linalg.norm(tangent)
        self.assertGreater(abs(vh[-1]@tangent),1-1e-12)
        for s in sections:
            np.testing.assert_allclose(np.linalg.norm(q[s][:,None]-q[s][None],axis=2),np.linalg.norm(current[s][:,None]-current[s][None],axis=2),atol=1e-14)
        self.assertEqual(proof['sourceFacesRemoved'],0)


if __name__=='__main__':unittest.main()
