import copy
import json
import math
import unittest

from geometry_tools.solver_diagnostics import json_safe_solver_diagnostics


class SolverDiagnosticSerializationTests(unittest.TestCase):
    def test_infeasible_status_and_nonfinite_objectives_survive_strict_json(self):
        source={'status':'PrimalInfeasible','feasible':False,
                'solverPrimalObjective':float('nan'),
                'solverDualObjective':float('inf'),
                'leadingHalfspaceDualRows':[{'row':7,'dual':float('-inf')}],
                'elapsedSeconds':.04}
        result=json_safe_solver_diagnostics(source)
        encoded=json.dumps(result,allow_nan=False)
        reloaded=json.loads(encoded)
        self.assertEqual(reloaded['status'],'PrimalInfeasible')
        self.assertFalse(reloaded['feasible'])
        self.assertIsNone(reloaded['solverPrimalObjective'])
        self.assertIsNone(reloaded['solverDualObjective'])
        self.assertIsNone(reloaded['leadingHalfspaceDualRows'][0]['dual'])
        self.assertEqual(reloaded['nonfiniteDiagnosticFields'],[
            {'path':'/solverPrimalObjective','classification':'nan'},
            {'path':'/solverDualObjective','classification':'+infinity'},
            {'path':'/leadingHalfspaceDualRows/0/dual','classification':'-infinity'}])
        self.assertTrue(reloaded['nonfiniteDiagnosticValues'])
        self.assertTrue(math.isnan(source['solverPrimalObjective']))

    def test_finite_solved_receipt_keeps_status_values_and_does_not_mutate(self):
        source={'status':'Solved','feasible':True,'iterations':3,
                'maximumResidualMeters':1e-10,'perRow':[0.,1.,None]}
        before=copy.deepcopy(source)
        result=json_safe_solver_diagnostics(source)
        self.assertEqual(source,before)
        for key,value in before.items():self.assertEqual(result[key],value)
        self.assertFalse(result['nonfiniteDiagnosticValues'])
        self.assertEqual(result['nonfiniteDiagnosticFields'],[])
        json.dumps(result,allow_nan=False)

    def test_rejects_candidate_array_or_missing_status_instead_of_repairing_it(self):
        with self.assertRaises(ValueError):
            json_safe_solver_diagnostics({'vertices':[[0.,float('nan'),0.]]})
        with self.assertRaises(ValueError):
            json_safe_solver_diagnostics({'status':'Solved','nonfiniteDiagnosticValues':True})


if __name__=='__main__':unittest.main()
