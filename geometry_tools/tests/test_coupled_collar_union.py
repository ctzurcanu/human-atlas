import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.coupled_collar_union import union_coupled_collars
from geometry_tools.open3d_backend import diagnostics

class CoupledCollarUnionTest(unittest.TestCase):
    def source(self,attributes=False):
        theta=np.arange(6)*2*np.pi/6
        centers=[[-1,0,0],[1,0,0],[-3,1.5,0],[-3,-1.5,0],[3,1.5,0],[3,-1.5,0]]
        v=np.concatenate([np.c_[np.full(6,c[0]),c[1]+.25*np.cos(theta),c[2]+.25*np.sin(theta)] for c in centers])*.001
        faces=[];cells=[]
        for a,b in [(0,1),(0,2),(0,3),(1,4),(1,5)]:
            first=len(faces)
            for i in range(6):
                j=(i+1)%6;faces.extend([[a*6+i,a*6+j,b*6+j],[a*6+i,b*6+j,b*6+i]])
            if centers[b][0]<centers[a][0]:faces[first:]=[f[::-1]for f in faces[first:]]
            cells.append(list(range(first,len(faces))))
        fields={'uvs':np.c_[v[:,1],v[:,2]].astype('<f4')}if attributes else {}
        return Mesh(v,np.asarray(faces),'fixture:metres',fields),cells,[list(range(6)),list(range(6,12))]

    def test_two_connected_junctions_keep_four_openings(self):
        m,cells,roots=self.source();q,p=union_coupled_collars(m,cells,roots)
        self.assertTrue(p['carrierAudit']['topological_continuity_pass']);self.assertEqual(p['carrierAudit']['component_euler'],[-2]);self.assertEqual(len(p['sourcePortRawIds']),4)
        self.assertEqual(p['rootClosuresRetained'],0);self.assertEqual(p['emittedLumenDisks'],0);self.assertTrue(set(p['newFaceSourceLineage'])<=set(range(len(m.faces))))
        self.assertEqual(diagnostics(q)['self_intersection_pairs'],0);np.testing.assert_array_equal(q.vertices[:len(m.vertices)],m.vertices)

    def test_source_fields_and_controls_are_literal_prefix(self):
        m,cells,roots=self.source(True);q,p=union_coupled_collars(m,cells,roots)
        np.testing.assert_array_equal(q.attributes['uvs'][:len(m.vertices)],m.attributes['uvs']);np.testing.assert_array_equal(q.vertices[:len(m.vertices)],m.vertices)
        np.testing.assert_allclose(q.attributes['uvs'][len(m.vertices):],q.vertices[len(m.vertices):,1:3],atol=1e-9)

    def test_missing_declared_root_is_rejected(self):
        m,cells,roots=self.source()
        with self.assertRaises((AssertionError,ValueError)):union_coupled_collars(m,cells,roots[:1])

if __name__=='__main__':unittest.main()
