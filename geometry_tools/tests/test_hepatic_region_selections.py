import importlib.util
from pathlib import Path
import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

spec=importlib.util.spec_from_file_location('hepatic_region_labels',Path(__file__).resolve().parents[2]/'server/hepatic-region-labels.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class RegionSelections(unittest.TestCase):
    def test_point_pinches_are_removed_without_reconstructing_source(self):
        v=np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
        f=np.array([[4,0,1],[4,1,2],[4,2,3],[4,3,0],[5,1,0],[5,2,1],[5,3,2],[5,0,3]])
        labels=np.array([0,1,0,1,0,1,0,1]);before_v=v.copy();before_f=f.copy()
        result,evidence=module.regularize_labels(v,f,labels)
        self.assertGreater(evidence['initial_pinched_fans'],0)
        self.assertEqual(evidence['remaining_pinched_fans'],0)
        self.assertEqual(set(result),{0,1})
        np.testing.assert_array_equal(v,before_v);np.testing.assert_array_equal(f,before_f)
        for color in [0,1]:self.assertEqual(topology(Mesh(v,f[result==color],'fixture:metres'))['invalid_vertex_links'],[])


if __name__=='__main__':unittest.main()
