"""Independent packed-source and failure-boundary checks for tip repairs."""
import json
from pathlib import Path
import unittest
import numpy as np

from geometry_tools.algorithms.terminal_profile_repair import terminal_profile_pose
from geometry_tools.model import AtlasReader
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.oriented_rays import boundary_loops

BASE=Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')
CURRENT=BASE/'blood-installed-v1905/circulation-atlas.json'


@unittest.skipUnless(CURRENT.is_file(),'Local unredistributed source fixture')
class TerminalProfileRepair(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.before_path=BASE/'blood-installed-v1889/circulation-atlas.json'
        cls.before=json.loads(cls.before_path.read_text());cls.after=json.loads(CURRENT.read_text())
        cls.spec=json.loads((BASE/'head-terminal-batch-spec-v1902.json').read_text())
        reader=AtlasReader();cls.meshes={}
        for p in cls.spec['profiles']:
            cls.meshes[p['partId']]=[reader.read_part(a,next(x for x in a['parts']if x['id']==p['partId']),
                frame='model:local-ta98:assembled-metres',catalogue_path=path)
                for a,path in [(cls.before,cls.before_path),(cls.after,CURRENT)]]

    def test_exact_other_anatomy_and_identity_preservation(self):
        ids=set(self.meshes);allowed={'positions','normals','normalSigns','uvs','colors','indices','chunk','bounds'}
        self.assertEqual([p['id']for p in self.before['parts']],[p['id']for p in self.after['parts']])
        for a,b in zip(self.before['parts'],self.after['parts']):
            if a['id']not in ids:self.assertEqual(a,b)
            else:self.assertEqual({k:v for k,v in a.items()if k not in allowed},{k:v for k,v in b.items()if k not in allowed})
        self.assertEqual(self.before['concepts'],self.after['concepts'])
        self.assertEqual(self.before['chunks'],self.after['chunks'][:-1])

    def test_shape_centroid_and_all_non_tip_controls(self):
        for profile in self.spec['profiles']:
            old,new=self.meshes[profile['partId']];tip=np.asarray(profile['tipIds']);fixed=np.ones(len(old.vertices),bool);fixed[tip]=False
            np.testing.assert_array_equal(old.vertices[fixed],new.vertices[fixed])
            np.testing.assert_array_equal(old.faces,new.faces)
            np.testing.assert_array_equal(old.source['source_vertex_ids'],new.source['source_vertex_ids'])
            np.testing.assert_allclose(old.vertices[tip].mean(0),new.vertices[tip].mean(0),atol=2e-9,rtol=0)
            np.testing.assert_allclose(np.linalg.norm(old.vertices[tip,None]-old.vertices[None,tip],axis=2),
                np.linalg.norm(new.vertices[tip,None]-new.vertices[None,tip],axis=2),atol=4e-9,rtol=0)
            self.assertLess(np.linalg.norm(new.vertices-old.vertices,axis=1).max(),.000045)

    def test_both_whole_meshes_have_no_crossings_and_retain_all_original_openings(self):
        for old,new in self.meshes.values():
            self.assertEqual(diagnostics(old,include_pairs=True)['self_intersection_pairs'],30)
            after=diagnostics(new,include_pairs=True);self.assertEqual(after['self_intersection_pairs'],0)
            self.assertEqual(after['nonmanifold_edges'],0)
            _,a,_,aa=boundary_loops(old);_,b,_,bb=boundary_loops(new)
            self.assertEqual({frozenset(x)for x in a},{frozenset(x)for x in b})
            for k in ['boundary_edges','closed_boundary_loops','inconsistent_shared_edge_winding']:self.assertEqual(aa[k],bb[k])

    def test_source_shading_assets_and_all_protected_joins_survive(self):
        for profile in self.spec['profiles']:
            old,new=self.meshes[profile['partId']];tip=profile['tipIds'];active=np.unique(old.faces[np.any(np.isin(old.faces,tip),axis=1)])
            fixed=np.ones(len(old.vertices),bool);fixed[active]=False
            for k in old.attributes:
                if k!='normals':np.testing.assert_array_equal(old.attributes[k],new.attributes[k])
            np.testing.assert_array_equal(old.attributes['normals'][fixed],new.attributes['normals'][fixed])
            t=new.triangles;n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);effective=new.attributes['normals'].astype(float)*new.attributes['normalSigns'][:,None]
            for vertex in active:self.assertTrue(np.all(n[np.any(new.faces==vertex,axis=1)]@effective[vertex]>0))
        proof=json.loads((CURRENT.parent/'source-mouth-replay.json').read_text())
        self.assertEqual(proof['summary']['verified'],260);self.assertEqual(proof['summary']['pending'],0)
        for edge in proof['interfaces']:self.assertTrue(edge['sourcePortProof']['passed'])

    def test_an_arbitrary_internal_profile_is_rejected(self):
        old=next(iter(self.meshes.values()))[0];_,loops,_,_=boundary_loops(old)
        boundary=set(int(i)for loop in loops for i in loop)
        internal=[i for i in range(len(old.vertices))if i not in boundary][:6]
        with self.assertRaises((AssertionError,ValueError,StopIteration)):
            terminal_profile_pose(old,internal)


if __name__=='__main__':unittest.main()
