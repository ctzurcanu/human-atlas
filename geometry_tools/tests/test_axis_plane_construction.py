import unittest
import numpy as np
from geometry_tools.implicit_cells import Arrangement
from geometry_tools.model import Mesh
from geometry_tools.audit import topology


class AxisPlaneConstructionTests(unittest.TestCase):
    def test_nearly_flat_construction_keeps_the_exact_plane_and_closed_cells(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.3,.3,1e-17]])
        fields=np.c_[np.full(4,-1.),p[:,0]-.25]
        a=Arrangement(p,fields,lambda x:'right' if x[1]>0 else 'left',axis_planes={1:(0,.25,1)})
        v,s,c=a.build(np.array([[0,1,2,3]]),skip_fields=[0])
        for name,uses in c.items():
            f=np.concatenate([s[u['surface']]if u['side']==1 else s[u['surface']][:,[0,2,1]]for u in uses])
            self.assertTrue(topology(Mesh(v,f,'test:metres'))['closed_oriented_manifold'],name)
        for point,value in zip(a.points,a.values):self.assertEqual(value[1],point[0]-.25)
        shared=[sid for sid in s if sum(any(u['surface']==sid for u in uses)for uses in c.values())==2]
        self.assertTrue(shared)
        for sid in shared:np.testing.assert_array_equal(v[np.unique(s[sid]),0],.25)


if __name__=='__main__':unittest.main()
