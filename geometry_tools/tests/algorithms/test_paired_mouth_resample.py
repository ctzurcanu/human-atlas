import unittest
import numpy as np
import manifold3d as md
from geometry_tools.algorithms.paired_mouth_resample import conform_paired_mouth
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

class PairedMouthTests(unittest.TestCase):
    def test_unequal_cycles_keep_native_geometry_and_volume(self):
        solid=md.Manifold.cylinder(1.,1.,1.,8)-md.Manifold.cylinder(1.,.5,.5,5)
        mesh=solid.to_mesh64();v=np.asarray(mesh.vert_properties)[:,:3];f=np.asarray(mesh.tri_verts,dtype=np.int64);top=f[np.all(abs(v[f,2]-1.)<1e-12,axis=1)];loops=boundary_loops(top);self.assertEqual(sorted(map(len,loops)),[5,8])
        points,faces,cap,rings,center,proof=conform_paired_mouth(v,f,top,loops,[0.,0.,1.])
        self.assertTrue(np.array_equal(points[:len(v)],v));self.assertEqual(len(rings[0]),len(rings[1]));self.assertTrue(proof['nativeControlsLiteral'])
        for old,new in zip(loops,rings):self.assertTrue(set(old)<=set(new))
        after=md.Manifold(md.Mesh64(np.ascontiguousarray(points),faces.astype(np.uint64)));self.assertEqual(after.status(),md.Error.NoError);self.assertAlmostEqual(after.volume(),solid.volume(),places=12)
        self.assertEqual({frozenset(r)for r in boundary_loops(cap)},{frozenset(r)for r in rings})

if __name__=='__main__':unittest.main()
