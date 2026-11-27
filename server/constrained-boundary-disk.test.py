"""Nonplanar closure, blocked disks and explicit opposed-face cancellation."""
from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from geometry_tools.algorithms.constrained_boundary_disk import constrained_boundary_disk,cancel_opposed_faces,BoundaryDiskError
from ta98_triangle_intersections import sweep_hits


def tube_fixture(nonplanar):
    n=12;angle=np.arange(n)*2*np.pi/n;top=.0006*np.sin(2*angle)if nonplanar else np.zeros(n)
    def ring(radius,z):return np.column_stack((radius*np.cos(angle),radius*np.sin(angle),np.broadcast_to(z,n)))
    v=np.r_[ring(.010,-.020),ring(.010,top),ring(.008,-.020),ring(.008,top)];faces=[]
    for i in range(n):
        j=(i+1)%n
        faces.extend([[i,j,n+i],[j,n+j,n+i],[2*n+j,2*n+i,3*n+i],[2*n+j,3*n+i,3*n+j],
                      [n+i,n+j,3*n+i],[n+j,3*n+j,3*n+i],[j,i,2*n+i],[j,2*n+i,2*n+j]])
    return v,np.array(faces,int),np.arange(3*n,4*n)[::-1]


class BoundaryDiskTests(unittest.TestCase):
    def test_planar_and_nonplanar_rims_keep_literal_geometry(self):
        for nonplanar in (False,True):
            v,f,ring=tube_fixture(nonplanar);before=v.copy();cap,proof=constrained_boundary_disk(v,ring,f)
            self.assertEqual(len(cap),len(ring)-2);self.assertEqual(proof['newControls'],0)
            self.assertTrue(np.array_equal(v,before));self.assertTrue(set(cap.ravel())<=set(ring))
            self.assertEqual(sweep_hits(v[cap],v[f],False)['properOrCoplanarIntersections'],0)
            self.assertEqual(sweep_hits(v[cap],v[cap],True)['properOrCoplanarIntersections'],0)

    def test_same_oriented_obstacle_is_not_silently_reused(self):
        v=np.array([[0.,0.,0.],[.01,0.,0.],[0.,.01,0.]])
        with self.assertRaises(BoundaryDiskError):constrained_boundary_disk(v,np.array([2,1,0]),np.array([[0,1,2]]))

    def test_only_exact_opposed_boundary_faces_cancel(self):
        f=np.array([[0,1,2],[0,2,1],[2,3,4]])
        out,proof=cancel_opposed_faces(f)
        self.assertEqual(proof['opposedPairsCancelled'],1);self.assertTrue(np.array_equal(out,[[2,3,4]]))
        with self.assertRaises(BoundaryDiskError):cancel_opposed_faces(np.array([[0,1,2],[1,2,0]]))


if __name__=='__main__':unittest.main()
