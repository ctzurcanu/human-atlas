import itertools,unittest
import numpy as np
from scipy.spatial import ConvexHull
from geometry_tools.algorithms.annular_pipe_envelopes import recover_paired_sheets,PipeEnvelopeError,signed_volume,face_components
from geometry_tools.model import Mesh
from geometry_tools.audit import audit


def three_mouth_wall():
    phi=(1+np.sqrt(5))/2
    p=np.array([(0,a,b*phi)for a in [-1,1]for b in [-1,1]]+[(a,b*phi,0)for a in [-1,1]for b in [-1,1]]+[(b*phi,0,a)for a in [-1,1]for b in [-1,1]],float)
    f=ConvexHull(p).simplices.copy();t=p[f];flip=np.einsum('ij,ij->i',np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),t.mean(1))<0
    f[flip]=f[flip][:,::-1]
    holes=next(ids for ids in itertools.combinations(range(len(f)),3)if len(np.unique(f[list(ids)]))==9)
    keep=np.ones(len(f),bool);keep[list(holes)]=False;n=len(p);faces=[*f[keep],*(f[keep][:,::-1]+n)]
    for face in f[list(holes)]:
        for a,b in zip(face,np.roll(face,-1)):faces.extend([[a,b,b+n],[a,b+n,a+n]])
    return np.r_[p,p*.8],np.array(faces,int),np.r_[np.ones(n,int),np.zeros(n,int)]


class BranchedWallEnvelopeTests(unittest.TestCase):
    def test_three_ports_retain_native_wall_and_share_only_classifier_closures(self):
        p,f,labels=three_mouth_wall();r=recover_paired_sheets(p,f,labels,expected_mouths=3)
        self.assertEqual(r.report['mouths'],3)
        self.assertEqual(face_components(f[r.cap_face_ids])[0],3)
        # Source coordinates/triangles are only exactly welded and reindexed.
        original={tuple(sorted(map(tuple,t)))for t in p[f]}
        self.assertEqual(original,{tuple(sorted(map(tuple,t)))for t in r.positions[r.original_faces]})
        self.assertGreater(signed_volume(r.positions,r.outer_faces),signed_volume(r.positions,r.fluid_faces))
        for faces in [r.outer_faces,r.fluid_faces]:
            check=audit(Mesh(r.positions,faces,'fixture:metres'),require_closed=True,expected_components=1,expected_genus=[0])
            self.assertTrue(check['passed'],check['failures'])
        for disk in r.fluid_disk_faces:
            self.assertTrue(any(np.array_equal(disk,t)for t in r.outer_faces))

    def test_branch_count_cannot_be_inferred_or_silently_relaxed(self):
        p,f,labels=three_mouth_wall()
        with self.assertRaisesRegex(PipeEnvelopeError,'declared number'):
            recover_paired_sheets(p,f,labels)
        with self.assertRaisesRegex(PipeEnvelopeError,'declared number'):
            recover_paired_sheets(p,f,labels,expected_mouths=4)
        with self.assertRaisesRegex(PipeEnvelopeError,'bounded'):
            recover_paired_sheets(p,f,labels,expected_mouths=513)


if __name__=='__main__':unittest.main()
