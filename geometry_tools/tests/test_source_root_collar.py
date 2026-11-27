import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_root_collar import union_source_root


def fixture(offset=1e-8):
    t=np.arange(6)*np.pi/3;root=np.c_[np.cos(t),np.zeros(6),np.sin(t)]*.0001
    vertices=np.r_[root,root+[0,-.0006,0],root+[0,.0006,0],root+[.0006,.0006,0],root]
    vertices[26:30,1]+=offset
    roots=[list(range(6)),list(range(6)),list(range(24,30))];faces=[];cells=[]
    for source,target in zip(roots,[list(range(6,12)),list(range(12,18)),list(range(18,24))]):
        start=len(faces)
        for i in range(6):
            j=(i+1)%6;faces.extend([[source[i],source[j],target[j]],[source[i],target[j],target[i]]])
        cells.append(list(range(start,len(faces))))
    faces=np.array(faces)
    # Positive-Y collars require opposite face winding to the negative-Y one.
    for ids in cells[1:]:faces[ids]=faces[ids][:,[0,2,1]]
    return Mesh(vertices,faces,'fixture-metres'),cells,roots


class SourceRootCollarTest(unittest.TestCase):
    def test_bounded_partial_source_chart_yields_three_exact_external_ports(self):
        m,cells,roots=fixture();v=m.vertices.copy()
        out,p=union_source_root(m,cells,roots,roots[0],maximum_correction_m=2e-7,maximum_root_area_relative_change=.005)
        np.testing.assert_array_equal(m.vertices,v);np.testing.assert_array_equal(out.vertices[:len(v)],v)
        self.assertEqual(len(p['rootSeamEdits']),4)
        self.assertEqual(p['carrierAudit']['mouths'],3)
        self.assertTrue(p['carrierAudit']['topological_continuity_pass'])
        self.assertEqual(p['rootClosuresRetained'],0)
        self.assertEqual(p['emittedLumenDisks'],0)

    def test_larger_source_discrepancy_is_rejected_without_mutation(self):
        m,cells,roots=fixture(1e-6);v=m.vertices.copy()
        with self.assertRaisesRegex(ValueError,'construction bound'):
            union_source_root(m,cells,roots,roots[0],maximum_correction_m=2e-7,maximum_root_area_relative_change=.005)
        np.testing.assert_array_equal(m.vertices,v)

    def test_no_majority_chart_cannot_be_inferred_from_proximity(self):
        m,cells,roots=fixture();roots=[roots[0],roots[2],roots[2]]
        with self.assertRaisesRegex(ValueError,'Two literal complete roots'):
            union_source_root(m,cells,roots,list(range(6)),maximum_correction_m=2e-7,maximum_root_area_relative_change=.005)


if __name__=='__main__':unittest.main()
