import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.directed_local_exit import directed_local_exit_planes

class DirectedLocalExitTest(unittest.TestCase):
    def index(self):
        theta=np.arange(12)*2*np.pi/12
        v=np.array([[.001*np.cos(t),y,.001*np.sin(t)]for y in [0.,.02]for t in theta]);f=[]
        for k in range(12):
            a=k;b=(k+1)%12;c=a+12;d=b+12;f.extend([[a,c,b],[b,c,d]])
        return SurfaceIndex(Mesh(v,np.array(f),'test:metres'))

    def test_one_consistent_direction_exits_both_sides_of_an_open_tube(self):
        index=self.index();before=index.mesh.geometry_sha;points=np.array([[.0002,.01,.0001],[-.0012,.01,.0001]])
        exits,proof=directed_local_exit_planes(index,points,np.array([0,1]),np.array([1.,0.,0.]),maximum_ray_m=.003)
        self.assertTrue(np.all(exits[:,0]>.0008))
        self.assertLess(float(index.closest_points(exits)[0]['distance_m'].max()),3e-9)
        self.assertEqual(proof['verifiedLocalProbes'],2)
        self.assertFalse(proof['closedSolidContainmentEstablished'])
        self.assertEqual(index.mesh.geometry_sha,before)

    def test_ray_bound_rejects_an_incomplete_exit(self):
        with self.assertRaisesRegex(ValueError,'finite ray bound'):
            directed_local_exit_planes(self.index(),np.array([[-.0012,.01,.0001]]),np.array([0]),np.array([1.,0.,0.]),maximum_ray_m=.00001)

if __name__=='__main__':unittest.main()
