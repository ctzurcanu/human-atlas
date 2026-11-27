"""Independent checks of the explicit PL patch deviation certificate."""
import importlib.util
import unittest
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('diagonal_check',Path(__file__).resolve().parents[2]/'scripts/install-current-head-diagonal-partials.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class SurfaceBound(unittest.TestCase):
    def setUp(self):
        self.old=np.array([[0,1,2],[0,2,3]])
        self.new=np.array([[1,2,3],[1,3,0]])

    def test_planar_surface_is_identical(self):
        v=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]])
        p=module.patch_check(v,self.old,self.new,1e-8)
        self.assertEqual(p['symmetricSurfaceDistanceUpperBound_m'],0)

    def test_warped_quad_bound_covers_both_projected_surfaces(self):
        v=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,.02]])
        p=module.patch_check(v,self.old,self.new,.03)
        t=v[self.old];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]).sum(0);normal/=np.linalg.norm(normal)
        distances=[]
        # Independently intersect rays through barycentric samples with the
        # opposite two triangles, in both directions of correspondence.
        for first,second in [(self.old,self.new),(self.new,self.old)]:
            for triangle in v[first]:
                for i in range(21):
                    for j in range(21-i):
                        point=triangle[0]+i/20*(triangle[1]-triangle[0])+j/20*(triangle[2]-triangle[0])
                        hits=[]
                        for other in v[second]:
                            u,w,h=np.linalg.solve(np.column_stack([other[1]-other[0],other[2]-other[0],-normal]),point-other[0])
                            if u>=-1e-10 and w>=-1e-10 and u+w<=1+1e-10:hits.append(abs(h))
                        self.assertTrue(hits);distances.append(min(hits))
        self.assertLessEqual(max(distances),p['symmetricSurfaceDistanceUpperBound_m']+1e-12)
        self.assertGreater(max(distances),.009)

    def test_reversed_or_changed_boundary_is_rejected(self):
        v=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]])
        for f in [self.new[:,::-1],np.array([[0,1,2],[1,2,3]])]:
            with self.assertRaises(AssertionError):module.patch_check(v,self.old,f,.001)

if __name__=='__main__':unittest.main()
