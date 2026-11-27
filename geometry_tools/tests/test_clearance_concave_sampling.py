import unittest
import numpy as np
from shapely import Polygon,constrained_delaunay_triangles
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.rigid_section_clearance import fit_section_clearance

class ConcaveSamplingTests(unittest.TestCase):
    def test_edge_sample_detects_concave_bridge_missed_by_controls_and_centroid(self):
        xy=np.array([[0,0],[3,0],[3,1],[1,1],[1,3],[0,3]],float)*.01
        v=np.r_[np.c_[xy,np.full(6,-.01)],np.c_[xy,np.full(6,.01)]];faces=[]
        for t in constrained_delaunay_triangles(Polygon(xy)).geoms:
            ids=[np.linalg.norm(xy-p,axis=1).argmin()for p in np.array(t.exterior.coords)[:3]]
            u=xy[ids[1]]-xy[ids[0]];w=xy[ids[2]]-xy[ids[0]]
            if u[0]*w[1]-u[1]*w[0]<0:ids=ids[::-1]
            faces.extend([ids[::-1],[i+6 for i in ids]])
        for i in range(6):
            j=(i+1)%6;faces.extend([[i,j,j+6],[i,j+6,i+6]])
        index=SurfaceIndex(Mesh(v,np.array(faces),'test:metres'))
        points=np.array([[.002,.025,0],[.002,.002,0],[.025,.002,0]]);f=np.array([[0,1,2]])
        fixed,old=fit_section_clearance(points,f,[],[(index,-1,'concave skin',.0001)],maximum_steps=1,wall_seconds=4.)
        np.testing.assert_array_equal(fixed,points);self.assertEqual(old['initial'][0]['violatingSamples'],0)
        moved,proof=fit_section_clearance(points,f,[],[(index,-1,'concave skin',.0001)],surface_barycentrics=[[1/3]*3,[.5,0,.5]],maximum_steps=1,wall_seconds=4.)
        self.assertGreater(proof['initial'][0]['violatingSamples'],0)
        self.assertGreater(np.linalg.norm(moved-points),0)
        self.assertLess(proof['steps'][0]['afterEnergy_m2'],proof['steps'][0]['energy_m2'])

    def test_invalid_surface_samples_are_rejected(self):
        with self.assertRaises(ValueError):fit_section_clearance(np.eye(3),np.array([[0,1,2]]),[],[],surface_barycentrics=[[.5,.5,.5]])
