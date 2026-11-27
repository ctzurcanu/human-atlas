import importlib.util
from pathlib import Path
import unittest

import numpy as np

from geometry_tools.model import Mesh


path=Path(__file__).resolve().parents[2]/'scripts/install-hand-source-baseline.py'
spec=importlib.util.spec_from_file_location('hand_source_baseline_installer',path)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ComponentMembershipTests(unittest.TestCase):
    def setUp(self):
        self.v=np.array([[0,0,0],[1,0,0],[0,1,0],[3,0,0],[4,0,0],[3,1,0]],float)
        self.f=np.array([[0,1,2],[3,4,5]])

    def test_spatial_reordering_does_not_change_source_component_membership(self):
        a=Mesh(self.v,self.f,'test:metres')
        moved=self.v.copy();moved[:3,0]+=6;moved[3:,0]-=6
        b=Mesh(moved,self.f,'test:metres')
        np.testing.assert_array_equal(module.component_membership(a),module.component_membership(b))

    def test_new_physical_weld_is_detected_even_with_literal_source_faces(self):
        a=Mesh(self.v,self.f,'test:metres')
        moved=self.v.copy();moved[3]=moved[0]
        b=Mesh(moved,self.f,'test:metres')
        self.assertFalse(np.array_equal(module.component_membership(a),module.component_membership(b)))


if __name__=='__main__':unittest.main()
