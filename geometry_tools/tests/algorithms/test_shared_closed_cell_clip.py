import unittest
import numpy as np
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank, split_closed_cell, ClosedCellClipError


def tube(n=20):
    angles=np.arange(n)*2*np.pi/n
    rings=[np.c_[r*np.cos(angles),np.full(n,y),r*np.sin(angles)]
           for r,y in [(2.,-1.),(2.,1.),(1.,-1.),(1.,1.)]]
    points=np.concatenate(rings);faces=[]
    for a,b,inward in [(0,n,False),(2*n,3*n,True)]:
        for i in range(n):
            j=(i+1)%n
            ts=[[a+i,b+i,a+j],[a+j,b+i,b+j]]
            faces.extend([t[::-1] if inward else t for t in ts])
    for outer,inner,reverse in [(0,2*n,False),(n,3*n,True)]:
        for i in range(n):
            j=(i+1)%n
            ts=[[outer+i,outer+j,inner+i],[outer+j,inner+j,inner+i]]
            faces.extend([t[::-1] if reverse else t for t in ts])
    return points,np.asarray(faces)


class ClosedTubeSplit(unittest.TestCase):
    def test_annular_cut_preserves_the_hole_and_volume(self):
        points,faces=tube();bank=SharedCutBank(points)
        cut=split_closed_cell(bank,faces,[0,1,0],.137,'section',allow_holes=True)
        self.assertEqual(len(cut.boundary_loops),2)
        self.assertTrue(cut.evidence['sourceControlsExact'])
        self.assertLess(abs(cut.evidence['volumeDifferenceMetersCubed']),1e-12)
        tri=np.asarray(bank.positions)[cut.interface_faces]
        # A genuine annular triangulation has no cap triangle containing
        # the centre of the lumen.
        for t in tri[:,:,[0,2]]:
            a,b=np.roll(t,-1,axis=0)-t,-t
            signs=a[:,0]*b[:,1]-a[:,1]*b[:,0]
            self.assertFalse(np.all(signs>=-1e-12) or np.all(signs<=1e-12))

    def test_unspecified_holes_are_rejected_without_mutating_bank(self):
        points,faces=tube();bank=SharedCutBank(points)
        with self.assertRaises(ClosedCellClipError):
            split_closed_cell(bank,faces,[0,1,0],.137,'section')
        self.assertEqual(len(bank.positions),len(points))
        self.assertFalse(bank.intersections)

    def test_opposed_source_faces_split_with_the_same_diagonals(self):
        points,faces=tube();bank=SharedCutBank(points)
        cut=split_closed_cell(bank,faces,[0,1,0],.137,'first',allow_holes=True)
        other=split_closed_cell(bank,faces[:,[1,2,0]],[0,1,0],.137,'first',allow_holes=True)
        canonical=lambda f: sorted(map(tuple,np.sort(f,axis=1)))
        self.assertEqual(canonical(cut.negative_faces),canonical(other.negative_faces))


if __name__=='__main__':unittest.main()
