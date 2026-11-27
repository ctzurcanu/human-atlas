import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_ring_cells import explicit_ring_cells
from geometry_tools.algorithms.coupled_collar_union import union_coupled_collars
from geometry_tools.open3d_backend import diagnostics

class SourceRingCellsTest(unittest.TestCase):
    def source(self,reconnect=False):
        theta=np.arange(6)*2*np.pi/6;centers=[[0,0,0],[-3,1.5,0],[-3,-1.5,0],[3,1.5,0],[3,-1.5,0]]
        v=np.concatenate([np.c_[np.full(6,c[0]),c[1]+.25*np.cos(theta),c[2]+.25*np.sin(theta)]for c in centers])*.001;faces=[]
        for a,b in [(0,1),(0,2),(0,3),(0,4)]+([(1,3)]if reconnect else []):
            first=len(faces)
            for i in range(6):
                j=(i+1)%6;faces.extend([[a*6+i,a*6+j,b*6+j],[a*6+i,b*6+j,b*6+i]])
            if centers[b][0]<centers[a][0]:faces[first:]=[f[::-1]for f in faces[first:]]
        return Mesh(v,np.asarray(faces),'fixture:metres')

    def test_four_source_branches_are_recovered_locally(self):
        m=self.source();r=explicit_ring_cells(m,arity=4);self.assertEqual(len(r['contracts']),1);c=r['contracts'][0];self.assertEqual(len(c['cells']),4)
        self.assertEqual({i for cell in c['cells']for i in cell['sourceFaceIds']},set(range(48)));self.assertFalse(c['clinicalAttachmentEstablished'])

    def test_remote_cycle_does_not_hide_a_local_branch(self):
        m=self.source(True);r=explicit_ring_cells(m,arity=4);self.assertEqual(len(r['contracts']),1);cells=r['contracts'][0]['cells'];self.assertEqual(len(cells),4)
        self.assertEqual({i for cell in cells for i in cell['sourceFaceIds']},set(range(48)),'Remote connector is retained, not swept into the local surgery')

    def test_fourway_carrier_keeps_four_ports_and_no_closure_disks(self):
        m=self.source();c=explicit_ring_cells(m,arity=4)['contracts'][0];q,p=union_coupled_collars(m,[x['sourceFaceIds']for x in c['cells']],[c['rootRawIds']]);self.assertTrue(p['carrierAudit']['topological_continuity_pass']);self.assertEqual(p['carrierAudit']['component_euler'],[-2]);self.assertEqual(len(p['sourcePortRawIds']),4)
        self.assertEqual(p['rootClosuresRetained'],0);self.assertEqual(p['emittedLumenDisks'],0);self.assertEqual(diagnostics(q)['self_intersection_pairs'],0);np.testing.assert_array_equal(q.vertices[:len(m.vertices)],m.vertices)

if __name__=='__main__':unittest.main()
