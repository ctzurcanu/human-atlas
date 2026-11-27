import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from geometry_tools.algorithms.bounded_float32 import untangle
from geometry_tools.algorithms.tube_sections import TubeGeometryError


class BoundedOrientationTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0.,0.,0.],[.001,.001,0.],[.002,.002+1e-10,0.],[.003,.004,0.]])
        self.q=self.p.astype(np.float32).astype(float)
        self.f=np.array([[0,1,2],[1,2,3]])
        self.words=np.broadcast_to(np.array([0,0,32767],np.int16),(2,3,3)).copy()

    def test_repairs_the_collapse_and_preserves_every_fixed_control(self):
        result,evidence=untangle(self.p,self.q,self.f,None,self.words,mutable_controls=[2],
                                original_controls=[0,1,3],additional_controls=[2],time_limit=2.)
        np.testing.assert_array_equal(result[[0,1,3]],self.q[[0,1,3]])
        t=result[self.f]
        self.assertTrue(np.all(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2]>0))
        np.testing.assert_array_equal(result.astype(np.float32).astype(float),result)
        self.assertTrue(evidence['orientation_only'])
        self.assertEqual(evidence['diagnosed_pairs_checked'],0)
        self.assertFalse(evidence['full_mesh_acceptance'])

    def test_immutable_control_cannot_be_repaired(self):
        with self.assertRaises(TubeGeometryError):
            untangle(self.p,self.q,self.f,None,self.words,mutable_controls=[],additional_controls=[2])

    def test_requires_actual_diagnosed_controls(self):
        with self.assertRaises(TubeGeometryError):
            untangle(self.p,self.q,self.f,None,self.words,mutable_controls=[2])

    def test_rejects_inward_source_fields(self):
        with self.assertRaises(TubeGeometryError):
            untangle(self.p,self.q,self.f,None,-self.words,mutable_controls=[2],additional_controls=[2])

    def test_rejected_collapsed_solver_proposal_is_excluded_then_rechecked(self):
        from scipy.optimize import milp
        calls=[]
        def invalid_first(*args,**kwargs):
            calls.append(True)
            if len(calls)>1:return milp(*args,**kwargs)
            # The nearest-grid coordinates collapse the first triangle. A
            # linear solver's success flag cannot establish exact geometry.
            bounds=kwargs['bounds'];values=np.zeros_like(bounds.lb)
            integer=np.flatnonzero(kwargs['integrality'])
            values[integer]=(bounds.lb[integer]+bounds.ub[integer])/2
            return SimpleNamespace(success=True,status=0,x=values,message='Collapsed proposal')
        with patch('scipy.optimize.milp',side_effect=invalid_first):
            result,evidence=untangle(self.p,self.q,self.f,None,self.words,
                mutable_controls=[2],original_controls=[0,1,3],additional_controls=[2],
                time_limit=3.,maximum_rounds=8,feasibility_only=True)
        np.testing.assert_array_equal(result[[0,1,3]],self.q[[0,1,3]])
        t=result[self.f]
        self.assertTrue(np.all(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2]>0))
        self.assertGreater(evidence['solver_rounds'],1)
        self.assertGreater(evidence['rejected_proposals'][0]['inward_or_collapsed_faces'],0)
        self.assertFalse(evidence['full_mesh_acceptance'])

    def test_round_budget_is_finite_and_explicit(self):
        for value in (0,33,True,1.5):
            with self.subTest(value=value),self.assertRaises(TubeGeometryError):
                untangle(self.p,self.q,self.f,None,self.words,mutable_controls=[2],
                    additional_controls=[2],maximum_rounds=value)

    def test_declared_reference_halfspace_is_checked_after_rounding(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.5,.500000001,0.]])
        q=p.astype(np.float32).astype(float);f=np.array([[0,1,3],[0,3,2]])
        words=np.broadcast_to(np.array([0,0,32767],np.int16),(2,3,3)).copy()
        result,evidence=untangle(p,q,f,None,words,mutable_controls=[3],additional_controls=[3],
            maximum_displacement_m=1e-6,time_limit=2.,
            point_halfspaces=[dict(control=3,normal=[0.,1.,0.],point=[0.,.5,0.],minimum=1e-8)])
        np.testing.assert_array_equal(result[:3],q[:3])
        self.assertGreaterEqual(result[3,1]-.5,1e-8)
        self.assertTrue(evidence['point_halfspace_checks'][0]['passed'])
        self.assertFalse(evidence['full_mesh_acceptance'])

    def test_halfspace_cannot_move_an_immutable_reference(self):
        with self.assertRaises(TubeGeometryError):
            untangle(self.p,self.q,self.f,None,self.words,mutable_controls=[2],additional_controls=[2],
                point_halfspaces=[dict(control=0,normal=[0.,1.,0.],point=[0.,0.,0.],minimum=1e-8)])

    def test_rejects_invalid_halfspace_inputs(self):
        for item in [dict(control=99,normal=[1,0,0],point=[0,0,0]),
                     dict(control=2,normal=[0,0,0],point=[0,0,0]),
                     dict(control=2,normal=[1,0,0],point=[np.nan,0,0])]:
            with self.subTest(item=item),self.assertRaises(TubeGeometryError):
                untangle(self.p,self.q,self.f,None,self.words,mutable_controls=[2],additional_controls=[2],point_halfspaces=[item])


if __name__=='__main__':unittest.main()
