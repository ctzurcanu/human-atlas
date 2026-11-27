import unittest
import numpy as np
from geometry_tools.algorithms.source_annular_partition import annular_face_partition


class SourceAnnularPartitionTests(unittest.TestCase):
    def source(self):
        faces=[]
        for j in range(3):
            for k in range(6):
                a=j*6+k;b=j*6+(k+1)%6;c=a+6;d=b+6
                faces.extend([[a,b,c],[b,d,c]])
        return np.array(faces),np.arange(24)

    def test_literal_ring_partitions_every_face_once(self):
        f,a=self.source();x,y,p=annular_face_partition(f,a,np.arange(6,12),root_face=0)
        np.testing.assert_array_equal(np.sort(np.r_[x,y]),np.arange(len(f)))
        self.assertEqual(set(a[f[x]].ravel())&set(a[f[y]].ravel()),set(range(6,12)))
        self.assertEqual(p['sourceFacesDuplicated'],0);self.assertFalse(p['coordinatesMoved'])

    def test_existing_open_end_cannot_be_used_as_an_internal_cut(self):
        f,a=self.source()
        with self.assertRaises(ValueError):annular_face_partition(f,a,np.arange(6),root_face=0)

    def test_disconnected_source_is_rejected(self):
        f,a=self.source();f=np.r_[f,f+24];a=np.arange(48)
        with self.assertRaises(ValueError):annular_face_partition(f,a,np.arange(6,12),root_face=0)
