"""A lumen may expose a protected component but must not erase it."""
import unittest
import numpy as np
import manifold3d as md
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.protected_assembly_carve import carve_with_protected_components


class ProtectedCarveTests(unittest.TestCase):
    def fixture(self):
        # Two adjacent cubes with a deliberately identical diagonal at x=1.
        p = np.array([[x,y,z] for x in (0.,1.,2.) for y in (0.,1.) for z in (0.,1.)])
        face = np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],
                         [0,4,5],[0,5,1],[2,3,7],[2,7,6],
                         [0,2,6],[0,6,4],[1,5,7],[1,7,3]])
        shared = face[2:4]
        a = np.delete(face, [2,3], axis=0)
        b = np.delete(face+4, [0,1], axis=0)
        surfaces = dict(a=a,b=b,shared=shared)
        cells = dict(a=[dict(surface='a',side=1),dict(surface='shared',side=1,neighbor='b')],
                     b=[dict(surface='b',side=1),dict(surface='shared',side=-1,neighbor='a')])
        normals = {}
        for name, f in surfaces.items():
            t=p[f]; n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]); n/=np.linalg.norm(n,axis=1)[:,None]
            normals[name]=np.repeat(np.rint(n*32767).astype('<i2')[:,None],3,axis=1)
        uvs={name:np.zeros((len(f),3,2),'<f4') for name,f in surfaces.items()}
        cutter=md.Manifold.cube((1.,.5,.5)).translate((.5,.25,.25)).to_mesh64()
        return p, surfaces, cells, normals, uvs, cutter

    def test_cut_host_exposes_shared_face_without_cutting_neighbour(self):
        p, surfaces, cells, normals, uvs, cutter = self.fixture()
        result=carve_with_protected_components(p,surfaces,cells,normals,uvs,['a'],
            cutter.vert_properties[:,:3],np.asarray(cutter.tri_verts))
        q, ss=result['positions'],result['surfaces']
        for name, expected in [('a',.875),('b',1.)]:
            f=np.concatenate([ss[u['surface']] if u['side']==1 else ss[u['surface']][:,[0,2,1]] for u in result['cells'][name]])
            solid,_=make_solid(q,f)
            self.assertAlmostEqual(solid.volume(),expected,places=12)
            edges=np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]
            _,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
            self.assertTrue(np.all(counts==2))
        exposed=[u for u in result['cells']['b'] if '/exposed:' in u['surface']]
        self.assertTrue(exposed)
        def area(f):
            t=q[f];return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1).sum()*.5
        self.assertAlmostEqual(sum(area(ss[u['surface']]) for u in exposed),.25,places=12)
        self.assertAlmostEqual(area(ss['shared']),.75,places=12)
        np.testing.assert_array_equal(q[:len(p)],p)
        self.assertEqual(result['evidence']['computational_fluid_cells_emitted'],0)

    def test_ancestry_mode_rejects_any_surviving_cap_fragment(self):
        p, surfaces, cells, normals, uvs, cutter = self.fixture()
        v, f = cutter.vert_properties[:, :3], np.asarray(cutter.tri_verts)
        outside_caps = np.flatnonzero(np.all(v[f, 0] == 1.5, axis=1))
        self.assertEqual(len(outside_caps), 2)
        result = carve_with_protected_components(p, surfaces, cells, normals, uvs,
            ['a'], v, f, classifier_cap_faces=outside_caps,
            prove_cap_exclusion_by_output_ancestry=True)
        self.assertEqual(result['evidence']['classifier_cap_contributions'], 0)
        with self.assertRaisesRegex(ValueError, 'Virtual cap contributed'):
            carve_with_protected_components(p, surfaces, cells, normals, uvs,
                ['a'], v, f, classifier_cap_faces=np.arange(len(f)),
                prove_cap_exclusion_by_output_ancestry=True)


if __name__=='__main__':unittest.main()
