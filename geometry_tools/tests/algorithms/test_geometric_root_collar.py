import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.algorithms.geometric_root_collar import subdivide_geometric_collar,planar_open_wall_volume
from geometry_tools.algorithms.root_collar_transport import transport_root_collar
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface
from geometry_tools.tests.algorithms.test_root_collar_transport import finite_pipe


class GeometricCollarTests(unittest.TestCase):
    def prepared(self):
        p,f,root=finite_pipe(n=9)
        z=clip_scalar_surface(p,f,p[:,2],.75)
        ff=np.r_[z['negative_faces'],z['positive_faces']]
        region=np.r_[np.ones(len(z['negative_faces']),bool),np.zeros(len(z['positive_faces']),bool)]
        return p,f,root,subdivide_geometric_collar(z['positions'],ff,region,
            normal=[0,0,1],origin=[0,0,0],start=.2,end=.75)

    def test_actual_paired_sections_are_coplanar_and_transport_is_isometric(self):
        p,f,root,z=self.prepared();q=z['vertices'];faces=z['faces']
        np.testing.assert_array_equal(q[:len(p)],p)
        r=Rotation.from_euler('y',6,degrees=True).as_matrix()
        moved=transport_root_collar(q,z['station'],r,[.01,0,0],[0,0,0])
        np.testing.assert_array_equal(moved[root],p[root])
        for profile in z['profiles']:
            ids=np.concatenate(profile['rings']);self.assertLess(np.ptp(q[ids,2]),1e-12)
            np.testing.assert_allclose(np.linalg.norm(moved[ids,None]-moved[ids][None],axis=2),
                np.linalg.norm(q[ids,None]-q[ids][None],axis=2),rtol=0,atol=2e-14)
        area=lambda t:np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1).sum()/2
        self.assertAlmostEqual(area(q[faces]),area(p[f]),places=12)

    def test_analytic_annulus_closure_matches_exact_polygonal_wall_volume(self):
        p,f,root,z=self.prepared();q=z['vertices'];wall=z['faces'][z['root_region']]
        rings=z['profiles'][-1]['rings'];volume=planar_open_wall_volume(q,wall,rings,origin=[0,0,0])
        exact=9*np.sin(2*np.pi/9)/2*(1-.8**2)*.75
        self.assertAlmostEqual(volume,exact,places=13)
        r=Rotation.from_euler('xyz',[35,-15,20],degrees=True).as_matrix();m=q@r.T+[10,5,-4]
        self.assertAlmostEqual(planar_open_wall_volume(m,wall,rings,origin=[10,5,-4]),exact,places=12)

    def test_offset_inner_outer_sections_cannot_masquerade_as_one_plane(self):
        _,_,_,z=self.prepared();q=z['vertices'].copy();rings=z['profiles'][-1]['rings']
        q[rings[0],2]+=.02
        with self.assertRaisesRegex(ValueError,'coplanar'):
            planar_open_wall_volume(q,z['faces'][z['root_region']],rings,origin=[0,0,0])


if __name__=='__main__':unittest.main()
