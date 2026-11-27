import unittest
import numpy as np
import manifold3d as md
from geometry_tools.model import Mesh
from geometry_tools.algorithms.oriented_shell_regions import occupied_shell_regions

def shells(solid):
    result=[]
    for part in solid.decompose():
        z=part.to_mesh64();result.append(Mesh(z.vert_properties[:,:3],z.tri_verts,'fixture:metres'))
    return result

class OrientedRegions(unittest.TestCase):
    def test_hollow_box_is_one_region_not_two_shell_volumes(self):
        a=md.Manifold.cube([4,4,4]);hole=md.Manifold.cube([2,2,2]).translate([1,1,1]);r=occupied_shell_regions(shells(a-hole),np.array([[.5,.5,.5],[2,2,2],[5,5,5.]]))
        self.assertEqual(len(r),1);self.assertEqual(r[0]['query_indices'],[0]);self.assertAlmostEqual(r[0]['volume_m3'],56)
    def test_positive_island_inside_negative_cavity_is_separate(self):
        a=md.Manifold.cube([6,6,6]);hole=md.Manifold.cube([4,4,4]).translate([1,1,1]);island=md.Manifold.cube([1,1,1]).translate([2,2,2]);r=occupied_shell_regions(shells((a-hole)+island),np.array([[.5,.5,.5],[1.5,1.5,1.5],[2.5,2.5,2.5]]))
        self.assertEqual(len(r),2);self.assertEqual(sorted(row['query_indices']for row in r),[[0],[2]]);self.assertAlmostEqual(sum(row['volume_m3']for row in r),153)
    def test_disjoint_solids_remain_separate(self):
        a=md.Manifold.cube([1,1,1])+md.Manifold.cube([1,1,1]).translate([2,0,0]);r=occupied_shell_regions(shells(a),np.array([[.5,.5,.5],[2.5,.5,.5],[1.5,.5,.5]]));self.assertEqual(sorted(row['query_indices']for row in r),[[0],[1]])

if __name__=='__main__':unittest.main()
