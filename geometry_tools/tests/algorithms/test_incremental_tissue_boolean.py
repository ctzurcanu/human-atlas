"""Later tissue additions must preserve the ancestry of earlier additions."""
import unittest
import numpy as np
import manifold3d as md
from geometry_tools.algorithms.finite_tissue_boolean import subtract_occupied,barycentric
from geometry_tools.algorithms.source_face_complement import conform_surface_boundary_edges
from geometry_tools.audit import topology
from geometry_tools.model import Mesh


def arrays(solid):
    m=solid.to_mesh64()
    return np.asarray(m.vert_properties)[:,:3],np.asarray(m.tri_verts,int)


class IncrementalTissueTests(unittest.TestCase):
    def test_old_contact_names_do_not_rebind_to_the_new_cutter(self):
        v,f=arrays(md.Manifold.cube([1.,1.,1.]))
        cv,cf=arrays(md.Manifold.hull_points([[.5,-1.,-1.],[.5,2.,-1.],[.5,-1.,2.],[2.,-1.,-1.]]))
        self.assertGreater(len(f),len(cf))
        for sid in ['vascular-contact:older-host:older-ligament','vascular-free:older-ligament']:
            result=subtract_occupied(v,{sid:f},{sid:dict(id=sid,kind='outer')},{'host':f},
                {'host':[(sid,i,1)for i in range(len(f))]},
                {'host':[dict(surface=sid,side=1)]},cv,cf,partition_free_outer=True)
            chart=result['charts'][sid]
            actual=result['positions'][result['surfaces'][sid]]
            replay=np.einsum('tck,tkj->tcj',chart['barycentrics'],v[f[chart['parentFaces']]])
            np.testing.assert_allclose(actual,replay,atol=1e-14,rtol=0)
            self.assertGreaterEqual(chart['parentFaces'].max(),len(cf))
            self.assertTrue(result['evidence']['vascularOuterReplay']['closedClassifierReplay'])

    def test_hanging_edge_conformity_retains_source_face_ancestry(self):
        v,f=arrays(md.Manifold.cube([1.,1.,1.]))
        a,b,c=f[0]
        mid=len(v)
        v=np.r_[v,[(v[a]+v[b])*.5]]
        source=np.r_[[[a,mid,c],[mid,b,c]],f[1:]]
        conformed,proof=conform_surface_boundary_edges(v,{'shell':source},
            {'whole':[dict(surface='shell',side=1)]},record_face_parents=True)
        output=conformed['shell']
        self.assertTrue(topology(Mesh(v,output,'fixture'))['closed_oriented_manifold'])
        self.assertFalse(proof['coordinatesMoved'])
        self.assertEqual(proof['positiveAreaCapsAdded'],0)
        parents=proof['sourceFaceParents']['shell']
        self.assertEqual(len(output),len(parents))
        areas=np.zeros(len(source))
        for triangle,parent in zip(output,parents):
            p=v[triangle]
            w=barycentric(p,v[source[parent]])
            np.testing.assert_allclose(w@v[source[parent]],p,atol=1e-14,rtol=0)
            self.assertGreaterEqual(w.min(),-1e-14)
            areas[parent]+=np.linalg.norm(np.cross(p[1]-p[0],p[2]-p[0]))/2
        p=v[source]
        expected=np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)/2
        np.testing.assert_allclose(areas,expected,atol=1e-14,rtol=0)


if __name__=='__main__':unittest.main()
