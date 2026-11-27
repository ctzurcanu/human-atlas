"""A declared construction resolution must not weaken finite cell checks."""
import unittest
import numpy as np
import manifold3d as md
from geometry_tools.algorithms.finite_tissue_boolean import subtract_occupied


class ConstructionResolutionTests(unittest.TestCase):
    @staticmethod
    def mesh(s):
        m=s.to_mesh64();return np.asarray(m.vert_properties)[:,:3],np.asarray(m.tri_verts,int)

    def cut(self,resolution):
        v,f=self.mesh(md.Manifold.cube((1.,1.,1.)));rv,rf=self.mesh(md.Manifold.cube((1.,1.,1.)).translate((.5,0.,0.)))
        return subtract_occupied(v,{'outer':f},{'outer':dict(id='outer',kind='outer')},{'host':f},
            {'host':[('outer',i,1)for i in range(len(f))]},{'host':[dict(surface='outer',side=1)]},
            rv,rf,construction_resolution_meters=resolution)

    def test_declared_precision_keeps_a_half_cube_and_its_actual_contact(self):
        for resolution in [None,1e-8]:
            r=self.cut(resolution);p=r['positions'][r['cellFaces']['host']]
            actual=float(np.einsum('ij,ij->i',p[:,0],np.cross(p[:,1],p[:,2])).sum()/6)
            self.assertAlmostEqual(actual,.5,places=13)
            self.assertEqual(len(r['vascularContacts']),1)
            self.assertEqual(r['evidence']['explicitConstructionResolutionMeters'],resolution)

    def test_nonfinite_or_nonpositive_construction_resolution_is_rejected(self):
        for resolution in [0.,-1.,float('nan'),float('inf')]:
            with self.assertRaisesRegex(ValueError,'positive measured construction resolution'):
                self.cut(resolution)


if __name__=='__main__':unittest.main()
