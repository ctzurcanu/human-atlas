import unittest
import numpy as np
from geometry_tools.algorithms.bounded_float32 import untangle
from geometry_tools.algorithms.tube_sections import TubeGeometryError
from geometry_tools.distance import triangle_pairs,exact_transverse_crossing


class BoundedFloat32Tests(unittest.TestCase):
    def fixture(self):
        p=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[.1,.900000001,0]],float)
        f=np.array([[0,1,2],[1,3,4]]);q=p.astype(np.float32).astype(float)
        w=np.broadcast_to(np.array([0,0,32767],np.int16),(2,3,3)).copy()
        return p,q,f,w

    def test_contact_is_repaired_with_fixed_clinical_controls(self):
        p,q,f,w=self.fixture()
        result,evidence=untangle(p,q,f,f[None],w,mutable_controls=[4],new_steps=8,maximum_displacement_m=1e-6,time_limit=3)
        np.testing.assert_array_equal(result[:4],q[:4])
        self.assertNotEqual(triangle_pairs(result[f[:1]],result[f[1:]])['relation'][0],'coplanar-overlap')
        self.assertLessEqual(evidence['maximum_coordinate_displacement_m'],1e-6)

    def test_no_selected_control_cannot_hide_a_failure(self):
        p,q,f,w=self.fixture()
        with self.assertRaises(TubeGeometryError):untangle(p,q,f,f[None],w,mutable_controls=[])

    def test_single_axis_budget_still_requires_actual_separation(self):
        p,q,f,w=self.fixture()
        result,evidence=untangle(p,q,f,f[None],w,mutable_controls=[4],new_steps=8,
            maximum_displacement_m=1e-6,time_limit=3,separation_axis_limit=1)
        np.testing.assert_array_equal(result[:4],q[:4])
        self.assertNotIn(triangle_pairs(result[f[:1]],result[f[1:]])['relation'][0],
                         ['crossing','coplanar-overlap'])
        self.assertEqual(evidence['remaining_diagnosed_crossings'],0)

    def test_feasibility_mode_retains_independent_geometry_guards(self):
        p,q,f,w=self.fixture();source=p.copy();stored=q.copy()
        result,evidence=untangle(p,q,f,f[None],w,mutable_controls=[4],new_steps=8,
            maximum_displacement_m=1e-6,time_limit=3,feasibility_only=True)
        np.testing.assert_array_equal(p,source)
        np.testing.assert_array_equal(q,stored)
        np.testing.assert_array_equal(result[:4],stored[:4])
        np.testing.assert_array_equal(result.astype(np.float32).astype(float),result)
        self.assertLessEqual(np.abs(result-stored).max(),1e-6)
        tri=result[f];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        self.assertTrue(np.all(np.einsum('tki,ti->tk',w,normal)>0))
        self.assertFalse(exact_transverse_crossing(tri[0],tri[1]))
        self.assertNotEqual(triangle_pairs(tri[:1],tri[1:])['relation'][0],'coplanar-overlap')
        self.assertEqual(evidence['objective_mode'],'feasibility')
        self.assertFalse(evidence['full_mesh_acceptance'])

    def test_feasibility_mode_cannot_change_immutable_controls(self):
        p,q,f,w=self.fixture()
        with self.assertRaises(TubeGeometryError):
            untangle(p,q,f,f[None],w,mutable_controls=[],feasibility_only=True)

    def test_feasibility_mode_rejects_inward_source_fields(self):
        p,q,f,w=self.fixture()
        with self.assertRaises(TubeGeometryError):
            untangle(p,q,f,f[None],-w,mutable_controls=[4],feasibility_only=True)


if __name__=='__main__':unittest.main()
