import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.curved_tissue_partition import split_graph_cells


def cube():
    p=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],float)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]],int)
    return p,f


class CurvedPartitionTests(unittest.TestCase):
    def test_curved_shared_surface_preserves_outer_carrier_and_closed_children(self):
        p,f=cube();height=lambda q:.1*np.sum(q*q,axis=1)
        r=split_graph_cells(p,{'whole':f},origin=np.zeros(3),basis=np.eye(3),height=height,
            plane_id='test-curve',maximum_interface_edge_m=.45,maximum_interface_faces=1000)
        self.assertTrue(np.array_equal(r.positions[:len(p)],p))
        for child in r.cells.values():
            self.assertTrue(audit(Mesh(r.positions,child,'test:metres'),require_closed=True,
                expected_components=1,expected_genus=[0])['passed'])
        self.assertLess(abs(r.evidence['cells']['whole']['volume_difference_m3']),1e-13)
        cap=r.interfaces['whole'];neg=r.cells[('whole',-1)];pos=r.cells[('whole',1)]
        self.assertTrue(np.array_equal(neg[-len(cap):],cap))
        self.assertTrue(np.array_equal(pos[-len(cap):],cap[:,[0,2,1]]))
        added=r.evidence['cells']['whole']['refinement']['inserted_controls']
        self.assertGreater(len(added),0)
        q=r.positions[added]
        np.testing.assert_allclose(q[:,2],height(q[:,:2]),atol=1e-14)

    def test_second_graph_splits_both_owners_without_cracks(self):
        p,f=cube();zero=lambda q:np.zeros(len(q))
        a=split_graph_cells(p,{'whole':f},origin=np.zeros(3),basis=np.eye(3),height=zero,
            plane_id='first',maximum_interface_edge_m=.6)
        basis=np.array([[0,1,0],[0,0,1],[1,0,0]],float)
        b=split_graph_cells(a.positions,a.cells,origin=np.zeros(3),basis=basis,height=zero,
            plane_id='second',maximum_interface_edge_m=.6)
        self.assertEqual(len(b.cells),4)
        for child in b.cells.values():
            self.assertTrue(audit(Mesh(b.positions,child,'test:metres'),require_closed=True,
                expected_components=1,expected_genus=[0])['passed'])
        # Every shared triangle in the full partition has reciprocal winding.
        from collections import defaultdict
        owners=defaultdict(list)
        for child in b.cells.values():
            for face in child:owners[tuple(sorted(face))].append(face)
        for faces in owners.values():
            if len(faces)==2:
                a1,a2=faces
                for x,y in zip(a1,np.roll(a1,-1)):
                    self.assertIn((y,x),list(zip(a2,np.roll(a2,-1))))

    def test_reflection_and_unbounded_refinement_are_rejected(self):
        p,f=cube();height=lambda q:.1*np.sum(q*q,axis=1)
        with self.assertRaisesRegex(ValueError,'proper rotation'):
            split_graph_cells(p,{'whole':f},origin=np.zeros(3),basis=np.diag([-1,1,1]),height=height,plane_id='bad')
        with self.assertRaisesRegex(ValueError,'budget'):
            split_graph_cells(p,{'whole':f},origin=np.zeros(3),basis=np.eye(3),height=height,
                plane_id='bounded',maximum_interface_edge_m=.01,maximum_interface_faces=4)


if __name__=='__main__':unittest.main()
