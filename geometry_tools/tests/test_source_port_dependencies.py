import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.source_port_dependencies import check_source_mouth


class SourcePortDependencyTests(unittest.TestCase):
    def pair(self):
        p=np.array([[0.,0,0],[.001,0,0],[.001,.001,0],[0,.001,0]])
        f=np.array([[0,1,2],[0,2,3]])
        return Mesh(p,f,'fixture'),Mesh(p+[0,0,2e-7],f[:,[0,2,1]],'fixture')

    def test_full_opposed_boundaries_match_within_the_declared_tolerance(self):
        a,b=self.pair();r=check_source_mouth(a,b,range(4),range(4))
        self.assertTrue(r['passed']);self.assertFalse(r['wholeLumenPatencyEstablished'])

    def test_partial_boundary_is_not_promoted_to_an_interface(self):
        a,b=self.pair()
        with self.assertRaises(ValueError):check_source_mouth(a,b,[0,1,2],[0,1,2])

    def test_detached_or_wrongly_wound_mouth_is_not_verified(self):
        a,b=self.pair();self.assertFalse(check_source_mouth(a,Mesh(b.vertices+.001,b.faces,b.frame),range(4),range(4))['passed'])
        self.assertFalse(check_source_mouth(a,Mesh(b.vertices,a.faces,b.frame),range(4),range(4))['passed'])

    def test_mismatched_coordinate_frames_are_rejected(self):
        a,b=self.pair()
        with self.assertRaises(ValueError):check_source_mouth(a,Mesh(b.vertices,b.faces,'other'),range(4),range(4))


if __name__=='__main__':unittest.main()
