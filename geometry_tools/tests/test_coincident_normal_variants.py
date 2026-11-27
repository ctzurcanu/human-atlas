import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.coincident_normal_variants import resolve_normal_variants


class NormalVariantTests(unittest.TestCase):
    def fixture(self,*,opposed=False,inward=False,uv=False):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0]]*2)
        n=np.array([[.2,0,1]]*3+[[0,0,1]]*3)
        if inward:n[0]*=-1
        attrs={'normals':n,'normalSigns':np.ones(6)}
        if uv:attrs['uvs']=np.array([[0.,0],[1,0],[0,1],[0,0],[2,0],[0,1]])
        return Mesh(p,np.array([[0,1,2],[3,5,4] if opposed else [3,4,5]]),'fixture',attrs)

    def test_selects_existing_normal_variant_and_preserves_support(self):
        m=self.fixture();before=m.attributes['normals'].copy();keep,proof=resolve_normal_variants(m)
        self.assertEqual(keep.tolist(),[1]);self.assertEqual(proof['duplicateCopiesRemoved'],1)
        np.testing.assert_array_equal(m.attributes['normals'],before)
        self.assertFalse(proof['appearanceIdenticalClaim'])

    def test_uv_conflict_retained(self):
        m=self.fixture(uv=True)
        keep,proof=resolve_normal_variants(m);self.assertEqual(len(keep),2)
        self.assertEqual(proof['retainedConflicts'][0]['reason'],'Other corner field conflict retained')

    def test_opposed_sheet_retained(self):
        m=self.fixture(opposed=True)
        keep,proof=resolve_normal_variants(m);self.assertEqual(len(keep),2)
        self.assertFalse(proof['opposedSheetsMerged'])

    def test_inward_variant_not_silently_resolved(self):
        m=self.fixture(inward=True)
        keep,proof=resolve_normal_variants(m);self.assertEqual(len(keep),2)
        self.assertTrue(proof['retainedConflicts'])

    def test_zero_normal_rejected(self):
        m=self.fixture();n=m.attributes['normals'].copy();n[0]=0
        bad=Mesh(m.vertices,m.faces,m.frame,dict(m.attributes,normals=n))
        with self.assertRaises(ValueError):resolve_normal_variants(bad)
