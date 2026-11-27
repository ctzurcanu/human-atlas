import tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.coherent_candidate_context import context_snapshots


class CoherentCandidateContextTests(unittest.TestCase):
    def test_joint_translation_does_not_compare_a_new_vessel_with_its_old_peer(self):
        v=np.array([[0.,0,0],[0,.001,0],[0,0,.001]])
        f=np.array([[0,1,2]])
        primary=Mesh(v,f,'test:metres');peer=Mesh(v+[.0002,0,0],f,'test:metres')
        moved=Mesh(v+[.0002,0,0],f,'test:metres');newpeer=Mesh(v+[.0004,0,0],f,'test:metres')
        with tempfile.TemporaryDirectory()as directory:
            oldpath,newpath=Path(directory)/'old.npz',Path(directory)/'new.npz';peer.save(oldpath);newpeer.save(newpath)
            rows={'peer':dict(path=str(newpath),sha256=sha(newpath),sourcePath=str(oldpath),sourceSha256=sha(oldpath))}
            before,after=context_snapshots(peer,'peer',rows)
            for p,c in [(primary,before),(moved,after)]:
                contacts,proof=pair_contacts(p,c,0,max_tests=100,max_pairs=100)
                self.assertTrue(proof['complete_within_pair_budget']);self.assertEqual(len(contacts['pairs']),0)
            contacts,_=pair_contacts(moved,peer,0,max_tests=100,max_pairs=100)
            self.assertGreater(len(contacts['pairs']),0,'Old peer produces a false overlap after their common motion')
            with self.assertRaises(ValueError):context_snapshots(newpeer,'peer',rows)
            self.assertIs(context_snapshots(peer,'fixed',rows)[1],peer)


if __name__=='__main__':unittest.main()
