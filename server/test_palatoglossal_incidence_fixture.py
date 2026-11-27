"""The retained PG witness must reject false planes and retain true ancestry."""
import json,sys,unittest
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(Path(__file__).parent))
from palatoglossal_incidence_fixture import common_parent_arrangement
from geometry_tools.algorithms.exact_source_intersections import contains,replay_boolean_vertices
W=ROOT/'.local-models/geometry-workbench/tongue-completion'

class ActualFailedIncidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=json.loads((W/'palatoglossal-mucosal-exact-incidence-v3/evidence.json').read_text())
        cls.fixture=cls.report['reunited_original_parent_fixture']
        cls.tris=[tuple(tuple(Q(x)for x in p)for p in t)for t in cls.fixture['corrected_fragment_triangles']]
        cls.source=tuple(tuple(Q(x)for x in p)for p in cls.fixture['exact_original_parent_vertices'])
    def test_actual_reported_binary_planes_still_rejected(self):
        z=np.load(W/'palatoglossal-mucosal-entry-v2-owner-checkpoint/failed-incidence-witness.npz');t=z['triangles'];point=z['point']
        # Every synthetic output face is incident to the exact retained failed
        # control; other controls are literal sources and cannot mask it.
        vertices=np.vstack([point,t[0,:2],t[1,:2],t[2,:2]])
        with self.assertRaisesRegex(ValueError,'No exact incident source intersection'):
            replay_boolean_vertices(vertices,np.array([[0,1,2],[0,3,4],[0,5,6]]),np.array([1,2,2]),np.array([0,0,1]),1,t[0],np.array([[0,1,2]]),2,t[1:].reshape(-1,3),np.array([[0,1,2],[3,4,5]]))
    def test_source_bound_reunion_has_exact_common_segment(self):
        axes,poly=common_parent_arrangement(self.tris,[('tongue-mucosa',405)]*3,self.source)
        self.assertEqual(len(poly),2)
        q=tuple(Q(x)for x in self.fixture['selected_common_control'])
        self.assertTrue(all(contains(q,t)for t in self.tris))
        self.assertTrue(contains(q,self.source))
        self.assertEqual(self.fixture['common_polygon_signed_chart_area'],'0')
    def test_distinct_ancestor_planes_cannot_be_reunited(self):
        with self.assertRaisesRegex(ValueError,'Distinct original source planes'):
            common_parent_arrangement(self.tris,[('tongue-mucosa',405),('tongue-mucosa',406),('tongue-mucosa',405)],self.source)
    def test_unreplayed_fragment_cannot_claim_original_plane(self):
        z=np.load(W/'palatoglossal-mucosal-entry-v2-owner-checkpoint/failed-incidence-witness.npz')
        binary=[tuple(tuple(Q.from_float(float(x))for x in p)for p in t)for t in z['triangles']]
        with self.assertRaisesRegex(ValueError,'has not been reconstructed'):
            common_parent_arrangement(binary,[('tongue-mucosa',405)]*3,self.source)
    def test_one_literal_encoding_has_analytical_half_ulp_bound(self):
        q=[Q(x)for x in self.fixture['selected_common_control']];p=self.fixture['selected_common_control_float']
        self.assertLess(self.fixture['selected_control_distance_from_boolean_point_m'],2e-10)
        for exact,value,bound in zip(q,p,self.fixture['analytical_encoding_bound_per_coordinate_m']):
            self.assertLessEqual(abs(Q.from_float(value)-exact),Q.from_float(bound))
        self.assertFalse(self.report['qualifying_candidate']);self.assertFalse(self.report['installed'])

if __name__=='__main__':unittest.main()
