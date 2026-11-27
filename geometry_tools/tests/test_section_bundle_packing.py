import unittest
import numpy as np
import shapely
from geometry_tools.algorithms.section_bundle_packing import pack_section_polygons


class SectionBundlePackingTest(unittest.TestCase):
    def test_primary_spatial_selection_precedes_irrecoverable_beam_truncation(self):
        square=np.array([[-.1,-.1],[.1,-.1],[.1,.1],[-.1,.1]])
        kwargs=dict(reference_pose=[0,0,0],maximum_displacement_m=20.,beam_width=2,allowed_masks=[[True,True,True],[True,False,False]])
        poses=np.array([[0.,0.,0.],[.1,0.,0.],[10.,0.,0.]])
        self.assertFalse(pack_section_polygons([square,square],poses,[],**kwargs)['feasible'])
        r=pack_section_polygons([square,square],poses,[],spread_primary_candidates=True,**kwargs)
        self.assertTrue(r['feasible']);self.assertEqual(r['selectedPoseIndices'],[2,0])

    def test_independent_profiles_fit_without_shrinking(self):
        square=np.array([[-.0004,-.0004],[.0004,-.0004],[.0004,.0004],[-.0004,.0004]])
        # Two separated source profiles cannot retain their 4 mm spacing in
        # this 3 mm corridor, but can move individually without radius loss.
        polygons=[square+[-.002,0],square+[.002,0]]
        obstacles=[np.array([[-.01,-.01],[-.0015,-.01],[-.0015,.01],[-.01,.01]]),np.array([[.0015,-.01],[.01,-.01],[.01,.01],[.0015,.01]])]
        poses=np.c_[np.linspace(-.003,.003,61),np.zeros(61),np.zeros(61)]
        result=pack_section_polygons(polygons,poses,obstacles,reference_pose=[0,0,0],beam_width=64)
        self.assertTrue(result['feasible'])
        self.assertGreaterEqual(result['minimumTissueClearance_m'],.00025)
        self.assertGreaterEqual(result['minimumMutualClearance_m'],.0001)
        self.assertLess(result['maximumAreaRelativeError'],1e-12)
        self.assertNotEqual(result['selectedPoseIndices'][0],result['selectedPoseIndices'][1])

    def test_no_invalid_closure_or_scaling(self):
        with self.assertRaises(ValueError):
            pack_section_polygons([[[0,0],[1,1],[0,1],[1,0]]],[[0,0,0]],[],reference_pose=[0,0,0])
        result=pack_section_polygons([[[0,0],[1,0],[1,1],[0,1]]],[[0,0,0]],[[[-1,-1],[2,-1],[2,2],[-1,2]]],reference_pose=[0,0,0])
        self.assertFalse(result['feasible'])
    def test_clear_endpoints_cannot_swap_through_each_other(self):
        square=np.array([[-.0004,-.0004],[.0004,-.0004],[.0004,.0004],[-.0004,.0004]])
        kwargs=dict(reference_pose=[0,0,0],allowed_masks=[[False,True],[True,False]])
        polygons=[square+[-.002,0],square+[.002,0]];poses=[[-.004,0,0],[.004,0,0]]
        self.assertTrue(pack_section_polygons(polygons,poses,[],**kwargs)['feasible'])
        self.assertFalse(pack_section_polygons(polygons,poses,[],preserve_translation_homotopy=True,**kwargs)['feasible'])
    def test_optional_cost_and_companion_distance_change_the_selected_layout(self):
        square=np.array([[-.0002,-.0002],[.0002,-.0002],[.0002,.0002],[-.0002,.0002]])
        poses=np.array([[0,0,0],[.001,0,0],[.005,0,0]])
        costs=np.array([[0,1,1],[1,1,0]],float)
        kwargs=dict(reference_pose=[0,0,0],extra_pose_costs=costs)
        distant=pack_section_polygons([square,square],poses,[],**kwargs)
        self.assertEqual(distant['selectedPoseIndices'],[0,2])
        nearby=pack_section_polygons([square,square],poses,[],pair_centre_limits=[(0,1,.002)],**kwargs)
        self.assertEqual(nearby['selectedPoseIndices'],[0,1])
    def test_invalid_optional_cost_is_rejected(self):
        with self.assertRaises(ValueError):
            pack_section_polygons([[[0,0],[1,0],[0,1]]],[[0,0,0]],[],reference_pose=[0,0,0],extra_pose_costs=[[float('nan')]])
    def test_explicit_compartment_connector_rejects_separated_pockets(self):
        square=np.array([[-.0002,-.0002],[.0002,-.0002],[.0002,.0002],[-.0002,.0002]])
        poses=np.array([[-.002,0,0],[.002,0,0],[-.002,.002,0]])
        masks=[[True,False,False],[False,True,False]]
        barrier=shapely.box(-.0005,-.001,.0005,.001)
        self.assertTrue(pack_section_polygons([square,square],poses,[],reference_pose=[0,0,0],allowed_masks=masks)['feasible'])
        self.assertFalse(pack_section_polygons([square,square],poses,[],reference_pose=[0,0,0],allowed_masks=masks,centre_connector_barriers=[(0,1,barrier)])['feasible'])
        masks[1]=[False,False,True]
        self.assertTrue(pack_section_polygons([square,square],poses,[],reference_pose=[0,0,0],allowed_masks=masks,centre_connector_barriers=[(0,1,barrier)])['feasible'])
    def test_secondary_variations_do_not_exhaust_primary_route_states(self):
        square=np.array([[-.0001,-.0001],[.0001,-.0001],[.0001,.0001],[-.0001,.0001]])
        poses=np.array([[0,0,0],[.005,0,0],[.001,0,0],[.00125,0,0],[.0015,0,0],[.00175,0,0],[.006,0,0]])
        masks=[[True,True,False,False,False,False,False],[False,False,True,True,True,True,True]]
        kwargs=dict(reference_pose=[0,0,0],allowed_masks=masks,beam_width=4,return_configurations=True)
        ordinary=pack_section_polygons([square,square],poses,[],**kwargs);self.assertEqual({c['poseIndices'][0]for c in ordinary['configurations']},{0})
        diverse=pack_section_polygons([square,square],poses,[],primary_diversity_m=.001,maximum_per_primary_cell=1,**kwargs);self.assertEqual({c['poseIndices'][0]for c in diverse['configurations']},{0,1})


if __name__=='__main__':
    unittest.main()
