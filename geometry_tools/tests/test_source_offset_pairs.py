"""Finite nonjoining source pairs share motion without becoming verified joins."""
import unittest
import numpy as np
from geometry_tools.algorithms.coupled_source_displacement import solve_coupled_source_displacement
from geometry_tools.algorithms.project_coupled_source_sections import propose_section_targets
from geometry_tools.source_port_dependencies import check_source_mouth
from geometry_tools.tests.test_coupled_source_displacement import pair,declaration
from geometry_tools.tests.test_project_coupled_source_sections import pose,targets


def fixture():
    meshes,rows=pair(gap=8e-5,alias=True)
    return meshes,[dict(rows[0],max_source_gap_m=8.1e-5,max_source_relative_area_error=.001)]


class SourceOffsetPairTests(unittest.TestCase):
    def test_two_sided_source_energy_allocates_free_boundary_motion_without_join_claim(self):
        meshes,pairs=fixture();shift=np.array([.0002,0.,0.]);original={p:m.vertices.copy()for p,m in meshes.items()}
        q,proof=solve_coupled_source_displacement(meshes,[],offset_pairs=pairs,
            prescribed={'b':declaration(meshes['b'],np.arange(18,24),shift)},
            fixed={'a':list(range(6))},declared_target_tolerance_m=1e-12,energy='source-triangle-fem')
        da=q['a'][18:24]-meshes['a'].vertices[18:24];db=q['b'][:6]-meshes['b'].vertices[:6]
        np.testing.assert_allclose(da,np.tile(shift/2,(6,1)),atol=5e-18,rtol=0)
        np.testing.assert_allclose(da,db,atol=5e-18,rtol=0)
        np.testing.assert_allclose(q['b'][:6]-q['a'][18:24],meshes['b'].vertices[:6]-meshes['a'].vertices[18:24],atol=3e-18,rtol=0)
        for pid,m in meshes.items():np.testing.assert_array_equal(m.vertices,original[pid])
        self.assertEqual(proof['interfaces'],[]);self.assertTrue(proof['offsetPairProofs'][0]['unjoined']);self.assertFalse(proof['offsetPairProofs'][0]['sourceJoinPassed'])
        self.assertFalse(check_source_mouth(meshes['a'],meshes['b'],pairs[0]['firstRing'],pairs[0]['secondRing'])['passed'])

    def test_common_proper_rigid_projection_rotates_offset_and_keeps_length_not_vector(self):
        meshes,pairs=fixture();q={pid:pose(m.vertices,scale=1.02)[0]for pid,m in meshes.items()}
        prescribed,_,proof=propose_section_targets(meshes,[],{},q,{},offset_pairs=pairs,maximum_uniform_scale=1.05)
        a,b=targets(meshes['a'],prescribed['a']),targets(meshes['b'],prescribed['b']);_,rotation=pose(meshes['a'].vertices)
        old=meshes['b'].vertices[:6]-meshes['a'].vertices[18:24];new=b[:6]-a[18:24]
        np.testing.assert_allclose(new,old@rotation.T,atol=4e-18,rtol=0)
        self.assertGreater(np.linalg.norm(new-old),1e-5);self.assertEqual(proof['rigidCells'],1)
        self.assertEqual(proof['cells'][0]['scale'],1.);self.assertEqual(proof['interfaces'],[])
        self.assertFalse(proof['offsetPairProofs'][0]['originalOffsetVectorPreserved']);self.assertFalse(proof['offsetPairProofs'][0]['sourceJoinPassed'])

    def test_offset_target_tolerance_cannot_weaken_verified_interface_tolerance(self):
        meshes,offset=fixture();verified,interfaces=pair(gap=2e-7);meshes.update(c=verified['a'],d=verified['b']);interfaces=[dict(interfaces[0],first='c',second='d')]
        shift=np.array([.0001,0.,0.]);prescribed={
            'a':declaration(meshes['a'],np.arange(18,24),shift),
            'b':declaration(meshes['b'],np.arange(6),shift+[2e-5,0,0])}
        q,proof=solve_coupled_source_displacement(meshes,interfaces,offset_pairs=offset,prescribed=prescribed,fixed={},declared_target_tolerance_m=1e-12,declared_offset_target_tolerance_m=3e-5)
        self.assertGreater(proof['maximumPrescribedReplayResidual_m'],1e-5)
        prescribed.update(c=declaration(meshes['c'],np.arange(18,24),shift),d=declaration(meshes['d'],np.arange(6),shift+[2e-9,0,0]))
        with self.assertRaisesRegex(ValueError,'agreement tolerance'):
            solve_coupled_source_displacement(meshes,interfaces,offset_pairs=offset,prescribed=prescribed,fixed={},declared_target_tolerance_m=1e-12,declared_offset_target_tolerance_m=3e-5)

    def test_incomplete_or_out_of_bounds_or_verified_overlap_pairs_reject_in_both_helpers(self):
        meshes,pairs=fixture();proposed={pid:m.vertices.copy()for pid,m in meshes.items()}
        wrong=[dict(pairs[0],max_source_gap_m=1e-6),dict(pairs[0],firstRing=list(range(18,23)),secondRing=list(range(5))),dict(pairs[0],max_source_relative_area_error=float('nan'))]
        for row in wrong:
            for mode in ['field','section']:
                with self.subTest(row=row,mode=mode),self.assertRaises(ValueError):
                    if mode=='field':solve_coupled_source_displacement(meshes,[],offset_pairs=[row],prescribed={},fixed={},declared_target_tolerance_m=1e-12)
                    else:propose_section_targets(meshes,[],{},proposed,{},offset_pairs=[row],maximum_uniform_scale=1.)
        tiny,verified=pair(gap=2e-7);overlap=[dict(verified[0],max_source_gap_m=1e-6,max_source_relative_area_error=.01)]
        with self.assertRaisesRegex(ValueError,'overlap'):
            solve_coupled_source_displacement(tiny,verified,offset_pairs=overlap,prescribed={},fixed={},declared_target_tolerance_m=1e-12)
        with self.assertRaisesRegex(ValueError,'overlap'):
            propose_section_targets(tiny,verified,{}, {p:m.vertices for p,m in tiny.items()}, {},offset_pairs=overlap,maximum_uniform_scale=1.)


if __name__=='__main__':unittest.main()
