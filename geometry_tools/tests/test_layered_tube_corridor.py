import unittest
import numpy as np
from geometry_tools.algorithms.layered_tube_corridor import layered_tube_corridor,route_station_correspondence


class LayeredTubeCorridorTests(unittest.TestCase):
    def test_containing_ball_fast_path_only_skips_already_clear_outlines(self):
        from geometry_tools.algorithms.section_configuration_field import SectionConfigurationField
        ring=np.array([[.001,0,0],[0,0,.001],[-.001,0,0],[0,0,-.001]])
        calls=[]
        def gap(p):
            calls.append(len(p));return .003-abs(p[:,0])
        field=SectionConfigurationField(ring[None],gap,ball_clearance=gap)
        value=field(np.array([[0,0,0],[.0025,0,0]]),np.zeros(2))
        self.assertEqual(field.boundingBallAccepted,1)
        self.assertEqual(calls,[2,9])
        self.assertGreater(value[0],0);self.assertLess(value[1],0)

    def test_clear_outline_cannot_encircle_a_blocked_course_centre(self):
        from geometry_tools.algorithms.section_configuration_field import SectionConfigurationField
        theta=np.arange(6)*np.pi/3
        ring=np.c_[.001*np.cos(theta),np.zeros(6),.001*np.sin(theta)]
        def gap(p):return np.linalg.norm(p[:,[0,2]],axis=1)-.0002
        self.assertTrue(np.all(gap(ring)>0))
        field=SectionConfigurationField(ring[None],gap,ball_clearance=gap)
        self.assertLess(field(np.zeros((1,3)),np.zeros(1))[0],0)
    def test_literal_flat_outline_fits_slot_without_reducing_its_wide_dimension(self):
        from geometry_tools.algorithms.section_configuration_field import SectionConfigurationField
        y=np.linspace(0,.01,11);guide=np.c_[np.zeros(len(y)),y,np.zeros(len(y))]
        theta=np.arange(6)*np.pi/3
        outline=np.c_[.001*np.cos(theta),np.zeros(6),.00015*np.sin(theta)]
        def gap(p):return .0003-abs(p[:,2])
        radii=np.full(len(y),.001)
        route,proof=layered_tube_corridor(guide,radii,gap,half_width_m=.001)
        self.assertIsNone(route)
        field=SectionConfigurationField(np.broadcast_to(outline,(len(y),6,3)),gap)
        route,proof=layered_tube_corridor(guide,radii,gap,half_width_m=.001,
                                          configuration_clearance=field)
        self.assertTrue(proof['complete']);self.assertTrue(proof['actualSectionOutlineProposal'])
        np.testing.assert_array_equal(route[[0,-1]],guide[[0,-1]])
        self.assertFalse(proof['sourceRadiiReduced'])
        np.testing.assert_array_equal(field.offsets[0],outline)

    def test_section_configuration_midpoints_reject_sheet_tunnelling(self):
        from geometry_tools.algorithms.section_configuration_field import SectionConfigurationField
        guide=np.c_[np.zeros(11),np.arange(11)*.001,np.zeros(11)]
        theta=np.arange(6)*np.pi/3
        ring=np.c_[.0001*np.cos(theta),np.zeros(6),.0001*np.sin(theta)]
        def gap(p):return np.abs(p[:,1]-.0055)-.00001
        field=SectionConfigurationField(np.broadcast_to(ring,(11,6,3)),gap)
        route,proof=layered_tube_corridor(guide,np.full(11,.0001),gap,half_width_m=.001,
                                          configuration_clearance=field)
        self.assertIsNone(route);self.assertEqual(proof['status'],'NoContinuousSampledCorridor')
    def test_saved_detour_preserves_source_stations_instead_of_total_length_fraction(self):
        guide=np.c_[np.zeros(5),np.arange(5)*.001,np.zeros(5)]
        # A lateral step at the first station lengthens only the proximal part.
        route=np.array([[0,0,0],[.0005,0,0],[.0005,.001,0],[0,.002,0],[0,.003,0],[0,.004,0]])
        mapped,proof=route_station_correspondence(guide,route,grid_m=.0005,half_width_m=.001)
        self.assertTrue(proof['sourceStationLastVisitsOrdered']);self.assertEqual(mapped[0],0.)
        self.assertGreater(mapped[1],.0014)
        self.assertNotAlmostEqual(mapped[1]/mapped[-1],.25,places=3)
        shifted=route.copy();shifted[2,0]+=.00001
        with self.assertRaisesRegex(ValueError,'literal adjacent-station'):
            route_station_correspondence(guide,shifted,grid_m=.0005,half_width_m=.001)

    def test_unfixed_donor_mouth_can_move_but_is_explicitly_reported_for_recoupling(self):
        y=np.linspace(0,.01,21);guide=np.c_[np.zeros(len(y)),y,np.zeros(len(y))]
        def gap(p):return abs(p[:,0])-.0005
        route,proof=layered_tube_corridor(guide,np.full(len(y),.0002),gap,
            half_width_m=.002,fixed_prefix=0,fixed_distal=False)
        self.assertTrue(proof['complete']);self.assertGreater(proof['proposedEndpointDisplacements_m'][0],.0005)
        self.assertEqual(proof['fixedPrefixStations'],0);self.assertFalse(proof['fixedDistal'])
        self.assertFalse(proof['anatomicalAcceptance'])

    def test_route_goes_around_fixed_box_with_original_endpoints_and_radius(self):
        y=np.linspace(0,.02,41);guide=np.c_[np.zeros(len(y)),y,np.zeros(len(y))]
        def clearance(p):
            # Exact signed distance to a fixed axis-aligned box.
            d=abs(p-[0,.01,0])-[.001,.002,.001]
            outside=np.linalg.norm(np.maximum(d,0),axis=1)
            return outside+np.minimum(d.max(1),0)
        radius=np.full(len(y),.0002)
        route,proof=layered_tube_corridor(guide,radius,clearance,half_width_m=.003)
        self.assertTrue(proof['complete']);np.testing.assert_array_equal(route[:2],guide[:2])
        np.testing.assert_array_equal(route[-1],guide[-1])
        self.assertGreater(np.linalg.norm(route[len(y)//2,[0,2]]),.001)
        self.assertTrue(np.all(clearance(route)>=radius));self.assertFalse(proof['sourceRadiiReduced'])

    def test_midpoint_screen_rejects_tunneling_through_an_unbroken_sheet(self):
        y=np.linspace(0,.02,21);guide=np.c_[np.zeros(len(y)),y,np.zeros(len(y))]
        def gap(p):return np.abs(p[:,1]-.0105)
        route,proof=layered_tube_corridor(guide,np.full(len(y),.0001),gap,half_width_m=.001)
        self.assertIsNone(route);self.assertEqual(proof['status'],'NoContinuousSampledCorridor')

    def test_reversible_graph_preserves_mouths_and_never_crosses_fixed_sheet(self):
        y=np.linspace(0,.02,21);guide=np.c_[np.zeros(len(y)),y,np.zeros(len(y))]
        def gap(p):return np.abs(p[:,1]-.0105)
        route,proof=layered_tube_corridor(guide,np.full(len(y),.0001),gap,
            half_width_m=.001,allow_lateral_detours=True)
        self.assertIsNone(route);self.assertEqual(proof['status'],'DisconnectedSampledClearanceGraph')
        def clear(p):return np.full(len(p),.002)
        route,proof=layered_tube_corridor(guide,np.full(len(y),.0001),clear,
            half_width_m=.001,allow_lateral_detours=True)
        self.assertTrue(proof['complete']);self.assertTrue(proof['lateralDetoursAllowed'])
        np.testing.assert_array_equal(route[:2],guide[:2]);np.testing.assert_array_equal(route[-1],guide[-1])
        self.assertFalse(proof['sourceRadiiReduced'])

    def test_graph_searches_all_explicitly_free_mouth_states(self):
        y=np.linspace(0,.02,41)
        guide=np.c_[np.clip((y-.006)/.008,0,1)*.002,y,np.zeros(len(y))]
        def gap(p):
            first=np.minimum(.0004-np.abs(p[:,0]),.009-p[:,1])
            second=.0004-np.abs(p[:,0]-.002)
            return np.maximum(first,second)
        route,proof=layered_tube_corridor(guide,np.full(len(y),.0001),gap,
            half_width_m=.003,allow_lateral_detours=True,fixed_prefix=0)
        self.assertTrue(proof['complete']);self.assertGreater(proof['eligibleStartStates'],1)
        self.assertGreater(route[0,0],.0015)
        self.assertGreater(proof['proposedEndpointDisplacements_m'][0],.0015)
        np.testing.assert_array_equal(route[-1],guide[-1])


if __name__=='__main__':unittest.main()
