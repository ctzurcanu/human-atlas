import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.incremental_source_metric_rows import incremental_source_metric_rows

class IncrementalSourceMetricTests(unittest.TestCase):
    def mesh(self,v):return Mesh(np.asarray(v,float),np.array([[0,1,2]]),'fixture:metres')
    def test_a_rigid_bend_does_not_consume_stretch_budget(self):
        a=self.mesh([[0,0,0],[1,0,0],[0,1,0]])
        angle=np.pi/4;rotation=np.array([[1,0,0],[0,np.cos(angle),-np.sin(angle)],[0,np.sin(angle),np.cos(angle)]])
        b=self.mesh(a.vertices@rotation.T)
        _,offset,p=incremental_source_metric_rows(a,b,np.arange(3),3)
        self.assertLess(np.linalg.norm(offset),1e-12)
        self.assertEqual(p['sourceFacetBudgets'],[.075])
    def test_inherited_strain_budget_is_not_reset_or_added(self):
        a=self.mesh([[0,0,0],[1,0,0],[0,1,0]])
        b=self.mesh([[0,0,0],[1.2,0,0],[0,1,0]])
        C,offset,p=incremental_source_metric_rows(a,b,np.arange(3),3)
        self.assertAlmostEqual(p['sourceFacetBudgets'][0],.2)
        self.assertAlmostEqual(np.linalg.norm(offset),.075)
        d=np.zeros((3,3));d[1,0]=.01
        self.assertGreater(np.linalg.norm((C@d.ravel()).reshape(-1,6)+offset),.075)
        d[1,0]=-.01
        self.assertLess(np.linalg.norm((C@d.ravel()).reshape(-1,6)+offset),.075)
if __name__=='__main__':unittest.main()
