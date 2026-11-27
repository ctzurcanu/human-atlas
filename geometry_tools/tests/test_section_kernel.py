import unittest
import numpy as np
from geometry_tools.algorithms.section_kernel import repair_section_kernel,repair_section_patch

class SectionKernel(unittest.TestCase):
    def test_coupled_patch_handles_swapped_interior_controls(self):
        p=np.array([[x,y]for x in range(4)for y in range(4)],float);f=[]
        for x in range(3):
            for y in range(3):
                a=x*4+y;b=a+4;c=a+1;d=b+1;f.extend([[a,b,d],[a,d,c]])
        f=np.array(f);fixed=np.flatnonzero(np.any((p==0)|(p==3),axis=1));p[[5,6]]=p[[6,5]]
        result,proof=repair_section_patch(p,f,fixed,minimum_double_area=.05,maximum_displacement_m=.9)
        self.assertTrue(proof['passed']);np.testing.assert_array_equal(result[fixed],p[fixed])
        self.assertLessEqual(np.linalg.norm(result-p,axis=1).max(),.9)

    def test_fixed_boundary_stays_literal_while_free_interior_is_repaired(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.],[1.1,.5]])
        f=np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
        result,proof=repair_section_kernel(p,f,range(4),minimum_double_area=.05,maximum_displacement_m=.3)
        self.assertTrue(proof['passed']);np.testing.assert_array_equal(result[:4],p[:4])
        self.assertLessEqual(np.linalg.norm(result-p,axis=1).max(),.3)
        self.assertGreaterEqual(proof['minimum_double_area_m2'],.05)

    def test_insufficient_distance_is_a_failure_not_a_relaxed_constraint(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.],[1.1,.5]])
        f=np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
        result,proof=repair_section_kernel(p,f,range(4),minimum_double_area=.05,maximum_displacement_m=.01)
        self.assertFalse(proof['passed']);np.testing.assert_array_equal(result,p)

if __name__=='__main__':unittest.main()
