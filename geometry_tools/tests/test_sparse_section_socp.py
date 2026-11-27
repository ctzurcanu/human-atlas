import json
import unittest
import numpy as np
from scipy import sparse
from geometry_tools.model import Mesh
from geometry_tools.algorithms.section_pair_separation import triangle_axis_sparse_rows,triangle_axis_rows
from geometry_tools.algorithms.sparse_section_socp import solve_sparse_section_socp,solve_independent_sparse_sections
from geometry_tools.algorithms.section_course_socp import solve_weighted_section_socp


class SparseSectionTests(unittest.TestCase):
    def test_independent_components_match_the_unsplit_objective(self):
        S=sparse.csr_matrix((0,2));A=sparse.csr_matrix([[1.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,1.,0.]])
        b=np.array([.001,.002]);full,_=solve_sparse_section_socp(S,A,b,maximum_motion_m=.003)
        separate,p=solve_independent_sparse_sections(S,A,b,maximum_motion_m=.003)
        self.assertTrue(p['feasible']);self.assertEqual(p['componentCount'],2)
        np.testing.assert_allclose(separate,full,atol=1e-9,rtol=0)
        _,p=solve_independent_sparse_sections(sparse.csr_matrix([[1.,-1.]]),A,b,maximum_motion_m=.003)
        self.assertEqual(p['componentCount'],1,'Objective coupling cannot be dropped')

    def test_failed_component_never_emits_a_whole_candidate(self):
        S=sparse.csr_matrix((0,2));A=sparse.csr_matrix([[1.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,1.,0.]])
        d,p=solve_independent_sparse_sections(S,A,[.001,.004],maximum_motion_m=.003)
        self.assertIsNone(d);self.assertFalse(p['feasible']);self.assertEqual(p['status'],'PrimalInfeasible')
        json.dumps(p,allow_nan=False)

    def test_same_coupled_calibrated_solution_as_dense_solver(self):
        S=sparse.csr_matrix([[1.,-1.]])
        H=np.eye(2)+(S.T@S).toarray();K=np.linalg.inv(H)
        w=np.array([[.25,.75]]);n=np.array([[1.,0.,0.]]);b=np.array([.0025]);A=sparse.csr_matrix((w[:,:,None]*n[:,None]).reshape(1,6))
        expected,p=solve_weighted_section_socp(K,w,n,b,maximum_motion_m=.003)
        actual,q=solve_sparse_section_socp(S,A,b,maximum_motion_m=.003)
        self.assertTrue(p['feasible']);self.assertTrue(q['feasible'])
        np.testing.assert_allclose(actual,expected,atol=1e-9,rtol=0)
        self.assertFalse(q['denseFaceByStationMatrixAllocated'])

    def test_infeasible_sparse_problem_emits_no_geometry_and_valid_json(self):
        d,p=solve_sparse_section_socp(sparse.csr_matrix((0,1)),sparse.csr_matrix([[1.,0.,0.]]),[.004],maximum_motion_m=.003)
        self.assertIsNone(d);self.assertFalse(p['feasible']);self.assertEqual(p['status'],'PrimalInfeasible')
        json.dumps(p,allow_nan=False)

    def test_sparse_station_reduction_preserves_the_full_vertex_feasible_set(self):
        a=Mesh(np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]),np.array([[0,1,2]]),'fixture')
        b=Mesh(a.vertices+[0.,0.,2.],a.faces,a.frame)
        sa=np.array([0,0,-1]);sb=np.array([1,1,1]);wa=np.array([[1.,0.],[1.,0.],[0.,0.]]);wb=np.tile([0.,1.],(3,1))
        W,N,rhs=triangle_axis_rows(a,b,wa,wb,[[0,0]],[0,0,1],gap_m=.1)
        A,bounds,p=triangle_axis_sparse_rows(a,b,sa,sb,[[0,0]],[0,0,1],station_count=2,gap_m=.1)
        self.assertLess(len(bounds),len(rhs));self.assertTrue(p['denseMatrixNotAllocated'])
        for value in [-3.,-1.,0.,1.,3.]:
            d=np.array([[value,0.,value],[0.,0.,-value]])
            self.assertEqual(bool(np.all(np.einsum('ni,ni->n',W@d,N)>=rhs)),bool(np.all(A@d.ravel()>=bounds)))


if __name__=='__main__':unittest.main()
