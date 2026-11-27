import unittest,tempfile,json
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.distance import triangle_pairs,point_triangle
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts,patch_area_bounds
from geometry_tools.registry import shape_descriptor
from geometry_tools.audit import topology,winding_contains,audit


def tetra(shift=(0,0,0)):
    return Mesh(np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])+shift,np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]]),'fixture')


class CoreTests(unittest.TestCase):
    def test_indexed_reader_world_offsets_and_attributes(self):
        with tempfile.TemporaryDirectory()as td:
            root=Path(td);(root/'public').mkdir();data=bytearray();record={'id':'fixture','chunk':0,'vertexCount':4,'indexCount':3,'geometryOrigin':[1,2,3],'sourceOffset':[.1,.2,.3]}
            for name,array in [('positions',np.array([[0,0,0],[1,0,0],[0,1,0],[999,999,999]],'<f4')),('indices',np.array([0,1,2],'<u4')),('normals',np.array([[0,0,32767]]*4,'<i2')),('normalSigns',np.array([1,-1,1,1],'i1')),('uvs',np.array([[0,0],[1,0],[0,1],[0,0]],'<f4'))]:
                while len(data)%array.dtype.itemsize:data.append(0)
                record[name]=len(data);data.extend(array.tobytes())
            (root/'public/mesh.bin').write_bytes(data);reader=AtlasReader(root);mesh=reader.read_part({'chunks':[{'url':'/mesh.bin'}]},record,frame='fixture')
            np.testing.assert_allclose(mesh.bounds,[[1.1,2.2,3.3],[2.1,3.2,3.3]])
            np.testing.assert_array_equal(mesh.attributes['normalSigns'],[1,-1,1]);self.assertEqual(mesh.attributes['uvs'].dtype,np.dtype('<f4'))
            restored=reader.reload(mesh.source,'fixture');self.assertEqual(restored.geometry_sha,mesh.geometry_sha)
            (root/'public/mesh.bin').write_bytes(data+b'bad')
            with self.assertRaises(ValueError):reader.reload(mesh.source,'fixture')

    def test_triangle_interior_distance_and_crossings(self):
        a=np.array([[[-1.,0,0],[1,0,0],[0,-1,0]]]);b=np.array([[[0,-.5,.01],[0,.5,.01],[0,0,1.]]])
        result=triangle_pairs(a,b);self.assertAlmostEqual(result['distance_m'][0],.01)
        # The closest pair is edge/face interior; no pair of vertices has this distance.
        self.assertGreater(np.linalg.norm(a[0,:,None]-b[0,None,:],axis=2).min(),.4)
        b[:,:,2]-=.02;r=triangle_pairs(a,b);self.assertEqual(r['distance_m'][0],0);self.assertEqual(r['relation'][0],'crossing')
        np.testing.assert_allclose(np.einsum('ni,nij->nj',r['first_barycentric'],a),r['first'],atol=1e-14)
        np.testing.assert_allclose(np.einsum('ni,nij->nj',r['second_barycentric'],b),r['second'],atol=1e-14)

    def test_coplanar_overlap_not_shared_edge(self):
        a=np.array([[[0.,0,0],[1,0,0],[0,1,0]]]);r=triangle_pairs(a,a)
        self.assertEqual(r['relation'][0],'coplanar-overlap');self.assertAlmostEqual(r['coplanar_overlap_area_m2'][0],.5)
        b=np.array([[[0.,0,0],[1,0,0],[0,-1,0]]]);self.assertEqual(triangle_pairs(a,b)['relation'][0],'touching')
        deg=np.array([[[2.,0,0],[2,0,0],[2,0,0]]]);r=triangle_pairs(a,deg);self.assertTrue(r['degenerate'][0]);self.assertAlmostEqual(r['distance_m'][0],1)
        segment=np.array([[[.2,.2,-1],[.2,.2,1],[.2,.2,-1]]]);r=triangle_pairs(segment,a);self.assertAlmostEqual(r['distance_m'][0],0);self.assertTrue(r['degenerate'][0])

    def test_bvh_has_no_missing_pairs_and_nearest_matches_brute_force(self):
        rng=np.random.default_rng(514);v=rng.normal(size=(90,3));f=np.arange(90).reshape(-1,3);a=Mesh(v,f,'fixture');b=Mesh(v+.21,f,'fixture');ba,bb=BVH(a,3),BVH(b,5)
        i,j=np.meshgrid(np.arange(30),np.arange(30),indexing='ij');pairs=np.c_[i.ravel(),j.ravel()];r=triangle_pairs(a.triangles[pairs[:,0]],b.triangles[pairs[:,1]])
        expected=set(map(tuple,pairs[r['distance_m']<=.08]));actual,_=pair_contacts(a,b,.08,first_bvh=ba,second_bvh=bb)
        self.assertEqual(expected,set(map(tuple,actual['pairs'])));self.assertAlmostEqual(ba.closest_mesh(bb)['distance_m'],r['distance_m'].min())
        expected={tuple(p)for p,d in zip(pairs,triangle_pairs(a.triangles[pairs[:,0]],a.triangles[pairs[:,1]])['distance_m'])if p[0]<p[1]and d<=.02}
        got,_=pair_contacts(a,a,.02,first_bvh=ba,second_bvh=ba,selfcheck=True);self.assertEqual(expected,set(map(tuple,got['pairs'])))

    def test_patch_area_is_bounded_and_budget_is_explicit(self):
        a=Mesh(np.array([[0.,0,0],[1,0,0],[0,1,0]]),np.array([[0,1,2]]),'fixture')
        b=Mesh(np.array([[0.,0,.01],[1,0,.01],[0,1,.01]]),a.faces,'fixture')
        area=patch_area_bounds(a,BVH(b),[0],.02);self.assertAlmostEqual(area['lower_m2'],.5);self.assertEqual(area['uncertainty_m2'],0)
        c=Mesh(np.array([[0.,0,0],[.5,0,0],[0,.5,0]]),a.faces,'fixture');r=patch_area_bounds(a,BVH(c),[0],.0001,absolute_area_error_m2=.005,max_depth=7)
        self.assertLessEqual(r['lower_m2'],.5*(.5+.0001*np.sqrt(2))**2+1e-14);self.assertGreaterEqual(r['upper_m2'],.125-1e-14)
        _,summary=pair_contacts(tetra(),tetra(),.01,max_pairs=1);self.assertTrue(summary['truncated']);self.assertFalse(summary['complete_within_pair_budget'])

    def test_descriptor_is_similarity_invariant(self):
        a=tetra();q=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
        b=Mesh(2.1*a.vertices@q.T+[10,-2,3],a.faces,a.frame)
        np.testing.assert_allclose(shape_descriptor(a),shape_descriptor(b),atol=1e-12)

    def test_topology_detects_holes_components_and_winding(self):
        m=tetra();t=topology(m);self.assertTrue(t['closed_oriented_manifold']);self.assertEqual(t['components'][0]['closed_orientable_genus'],0)
        r=winding_contains(m,[[.1,.1,.1],[2,2,2],[0,0,0]]);self.assertEqual(r['inside'],[True,False,False]);self.assertTrue(r['boundary_ambiguous'][2])
        openm=Mesh(m.vertices,m.faces[:3],m.frame);self.assertEqual(topology(openm)['open_edges'],3)
        with self.assertRaises(ValueError):winding_contains(openm,[[.1,.1,.1]])
        self.assertTrue(audit(m,expected_components=1,expected_genus=[0],require_closed=True)['passed'])


if __name__=='__main__':unittest.main()
