import unittest

import numpy as np

from geometry_tools.algorithms.declared_mouth_area_budgets import reconcile_declared_mouths
from geometry_tools.model import Mesh
from geometry_tools.tests.test_declared_mouth_batch import DeclaredMouthBatchTests


class SeparateMouthAreaBudgetTests(unittest.TestCase):
    def source(self, area_ratio):
        meshes,pair=DeclaredMouthBatchTests().fixtures()
        m=meshes['second'];p=m.vertices.copy();p[:,:2]*=np.sqrt(area_ratio)
        meshes['second']=Mesh(p,m.faces,m.frame,m.attributes,m.source)
        return meshes,[pair]

    def solve(self, area_ratio):
        meshes,pairs=self.source(area_ratio)
        return reconcile_declared_mouths(meshes,pairs,maximum_gap_m=3e-5,
            maximum_motion_m=2e-5,maximum_relative_area_change=.025,
            maximum_source_area_discrepancy=.1,maximum_facet_rotation_degrees=10)

    def test_source_discrepancy_does_not_equal_physical_area_change(self):
        _,proof=self.solve(1.036)
        self.assertLess(max(proof['pairs'][0]['relativeSectionAreaChanges']),.025)
        self.assertTrue(proof['separateInputAndOutputAreaBudgets'])

    def test_larger_input_budget_never_authorizes_larger_physical_change(self):
        with self.assertRaisesRegex(ValueError,'physical budget'):
            self.solve(1.08)


if __name__=='__main__':unittest.main()
