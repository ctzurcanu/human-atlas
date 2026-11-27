"""Shared-edge refinement keeps closed material topology and occupied volume."""
import unittest
import numpy as np
from geometry_tools.algorithms.refine_mapped_surface import refine_mapped_surface
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

class SharedRefinement(unittest.TestCase):
    def test_every_edge_selection_of_closed_tetrahedron(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        edges=np.unique(np.sort(np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]],axis=1),axis=0)
        before=topology(Mesh(p,f,'fixture'))
        for mask in range(1,64):
            selection=edges[[bool(mask&(1<<i))for i in range(6)]]
            q=refine_mapped_surface(p,p,f,selection,lambda x:x,centroid_fan=False)
            after=topology(Mesh(q['positions'],q['faces'],'fixture'))
            self.assertTrue(after['closed_oriented_manifold'],mask)
            self.assertAlmostEqual(after['oriented_triangle_integral_m3'],before['oriented_triangle_integral_m3'],delta=1e-14,msg=mask)
            self.assertTrue(np.array_equal(q['source_positions'][:4],p))
            src=p[f[q['parent_faces']]];child=q['source_positions'][q['faces']]
            self.assertTrue(np.all(np.einsum('fi,fi->f',np.cross(src[:,1]-src[:,0],src[:,2]-src[:,0]),np.cross(child[:,1]-child[:,0],child[:,2]-child[:,0]))>0),mask)

if __name__=='__main__':unittest.main()
