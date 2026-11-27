import unittest
import numpy as np
from geometry_tools.algorithms.bone_radial_field import measured_bone_radial_field,cached_bone_radial_field
from geometry_tools.tests.test_sparse_group_clearance import cube


class MeasuredRadialFieldTests(unittest.TestCase):
    def test_actual_inner_and_outer_surfaces_determine_mid_corridor(self):
        axis=np.array([[0.,0.,0.],[0.,0.,.03],[0.,0.,.06],[0.,0.,.08]])
        box=lambda r:cube([-r,-r,-.03],[r,r,.12])
        field,proof=measured_bone_radial_field(box(.004),box(.008),axis,axis,[0,1,0],[0,1,0],
            source_skin=box(.012),target_skin=box(.016))
        scale=field(np.array([1.5]),np.array([[.008,0,0]]))
        self.assertAlmostEqual(float(scale[0])*.008,.012,delta=.00005)
        self.assertGreater(proof['pairedBoneSkinIntervals']['measuredPairedCorridors'],20000)
        self.assertFalse(proof['anatomicalAcceptance'])

    def test_one_sided_envelope_is_rejected(self):
        axis=np.array([[0.,0.,0.],[0.,0.,.03],[0.,0.,.06]])
        box=cube([-.004,-.004,-.02],[.004,.004,.08])
        with self.assertRaises(ValueError):
            measured_bone_radial_field(box,box,axis,axis,[0,1,0],[0,1,0],source_skin=box)

    def test_bone_interior_does_not_inherit_a_reversed_annulus(self):
        axis=np.array([[0.,0.,0.],[0.,0.,.03],[0.,0.,.06],[0.,0.,.08]])
        box=lambda r:cube([-r,-r,-.03],[r,r,.12])
        field,proof=measured_bone_radial_field(box(.004),box(.002),axis,axis,[0,1,0],[0,1,0],
            source_skin=box(.008),target_skin=box(.020))
        radii=np.array([.001,.002,.004,.006,.008])
        mapped=radii*field(np.full(5,1.5),np.c_[radii,np.zeros((5,2))])
        self.assertTrue(np.all(mapped>0))
        self.assertTrue(np.all(np.diff(mapped)>0))
        self.assertAlmostEqual(mapped[1],.001,delta=.00002)
        self.assertGreater(proof['pairedBoneSkinIntervals']['interiorBoneQuerySamples'],0)
        self.assertTrue(proof['pairedBoneSkinIntervals']['sourceBoneEntryAnatomyNotInferred'])
        cached=cached_bone_radial_field(field.sample_bank)
        np.testing.assert_allclose(cached(np.full(5,1.5),np.c_[radii,np.zeros((5,2))]),
            field(np.full(5,1.5),np.c_[radii,np.zeros((5,2))]),atol=1e-14)


if __name__=='__main__':unittest.main()
