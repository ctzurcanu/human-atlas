"""Analytical distances and an actual near-obstacle tube detour."""
from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.algorithms.section_course_clearance import closest_triangles,smooth_course_metric,project_course_halfspaces,propose_section_detour

def triangle_meets_cube(triangle,radius=.005):
    """Independent convex clipping, including edge/vertex contacts."""
    poly=list(triangle)
    for axis in range(3):
        for sign in [-1,1]:
            if not poly:return False
            output=[]
            for a,b in zip(poly,poly[1:]+poly[:1]):
                da,db=sign*a[axis]-radius,sign*b[axis]-radius
                if da<=0:output.append(a)
                if (da<0<db)or(db<0<da):output.append(a+(b-a)*da/(da-db))
            poly=output
    return bool(poly)

class SectionCourseTests(unittest.TestCase):
    def test_plane_edge_and_corner_distances(self):
        triangles=np.array([[[0.,0.,0.],[2.,0.,0.],[0.,2.,0.]]])
        points=np.array([[.5,.5,3.],[2.,2.,1.],[-1.,0.,0.],[0.,-1.,2.]])
        expected=np.array([[.5,.5,0.],[1.,1.,0.],[0.,0.,0.],[0.,0.,0.]])
        actual,face=closest_triangles(points,triangles,batch_size=2)
        np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-14)
        np.testing.assert_array_equal(face,np.zeros(4,int))

    def test_fixed_port_cannot_be_released(self):
        K=smooth_course_metric(9,fixed=[0,8])
        d,_=project_course_halfspaces(K,[4],[[0.,0.,1.]],[.003])
        np.testing.assert_array_equal(d[[0,8]],np.zeros((2,3)))
        self.assertAlmostEqual(d[4,2],.003,places=14)
        with self.assertRaises(ValueError):project_course_halfspaces(K,[0],[[1.,0.,0.]],[.001])

    def test_nonzero_area_sliver_is_not_misclassified_as_degenerate(self):
        triangle=np.array([[[0.,0.,0.],[1.,0.,0.],[1.,1e-9,0.]]])
        point=np.array([[.75,.25e-9,.003]])
        actual,_=closest_triangles(point,triangle)
        np.testing.assert_allclose(actual,[[.75,.25e-9,0.]],rtol=0.,atol=1e-14)

    def test_actual_tube_triangles_clear_cube_without_scaling_sections(self):
        cube=.005*np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
        cf=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        theta=np.arange(12)*2*np.pi/12;centers=np.c_[np.full(41,.0055),np.linspace(-.025,.025,41),np.zeros(41)]
        rings=np.arange(41*12).reshape(41,12);v=(centers[:,None]+.001*np.c_[np.cos(theta),np.zeros(12),np.sin(theta)][None]).reshape(-1,3)
        f=[]
        for j in range(40):
            for k in range(12):
                a=j*12+k;b=j*12+(k+1)%12;c=b+12;d=a+12;f.extend([[a,c,b],[a,d,c]])
        f=np.array(f);source=v.copy();obstacle=cube[cf]
        contains=lambda p:(np.all(abs(p)<.005,axis=1),np.zeros(len(p),bool))
        # This fixture approaches the cube's +X face. Convex clipping tests the
        # complete actual tube triangles, including the aligned-edge case that
        # a proper-face-crossing count alone can miss.
        def hit(a,b):
            faces=[j for j,t in enumerate(b)if triangle_meets_cube(t)]
            return dict(properOrCoplanarIntersections=len(faces),examples=[(6,j)for j in faces])
        # Some crossings run exactly along source triangle edges, so a proper
        # face-pair count can be zero despite this known occupied overlap.
        self.assertTrue(np.any(np.all(abs(v)<.005,axis=1)))
        moved,_,proof=propose_section_detour(v,f,rings,obstacle,contains=contains,crossings=hit,origin=[0.,0.,0.],fixed_stations=[0,1,39,40],maximum_motion_m=.002)
        self.assertIsNone(proof['rejectedReason']);self.assertEqual(hit(obstacle,moved[f])['properOrCoplanarIntersections'],0)
        self.assertFalse(np.any(np.all(abs(moved)<.005,axis=1)))
        np.testing.assert_array_equal(v,source);np.testing.assert_array_equal(moved[rings[[0,1,39,40]]],v[rings[[0,1,39,40]]])
        self.assertLess(proof['maximumMotionMeters'],.002)

if __name__=='__main__':unittest.main()
