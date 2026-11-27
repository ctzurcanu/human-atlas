import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.open_envelope_opposed_facets import cancel_open_envelope_internal_facets

def source(*,z=2.,outside=False,partition=False,same_winding=False):
    ring=np.array([[1.,0.],[0.,1.],[-1.,0.],[0.,-1.]])
    v=np.array([[x,y,height] for height in [0.,2.,4.] for x,y in ring])
    faces=[]
    for layer in range(2):
        for j in range(4):
            a=4*layer+j;b=4*layer+(j+1)%4;c=b+4;d=a+4
            faces += [[a,b,c],[a,c,d]]
    if partition:
        patch=[[4,5,6],[4,6,7]]
    else:
        x=1.7 if outside else 0.
        v=np.r_[v,[[x-.1,-.1,z],[x+.1,-.1,z],[x,.1,z]]]
        patch=[[12,13,14]]
    for f in patch:faces += [f,f if same_winding else f[::-1]]
    return Mesh(v,np.array(faces),"test:metres",{"colors":np.arange(len(v)*3,dtype=np.uint8).reshape(-1,3)})

class OpenEnvelopeOpposedFacetTests(unittest.TestCase):
    def test_interior_floating_sheet_cancelled_without_closing_source_openings(self):
        m=source();keep,p=cancel_open_envelope_internal_facets(m)
        np.testing.assert_array_equal(keep,np.arange(16))
        self.assertEqual(p['sourceOpenings'],2)
        self.assertEqual(p['computationCapsEmitted'],0)
        self.assertTrue(p['sourceBoundaryChainLiteral'])

    def test_shared_internal_partition_removed_with_original_tube_faces_literal(self):
        m=source(partition=True);keep,p=cancel_open_envelope_internal_facets(m)
        np.testing.assert_array_equal(m.faces[keep],m.faces[:16])
        self.assertEqual(len(p['removedSourceFaceIds']),4)
        self.assertEqual(p['retainedEnvelopeSelfPairs'],0)

    def test_external_opposed_faces_are_not_removed(self):
        with self.assertRaisesRegex(ValueError,'not inside'):
            cancel_open_envelope_internal_facets(source(outside=True))

    def test_mouth_adjacent_ambiguity_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'closure-independent'):
            cancel_open_envelope_internal_facets(source(z=.1))

    def test_same_winding_copies_do_not_qualify(self):
        with self.assertRaisesRegex(ValueError,'No exactly opposed'):
            cancel_open_envelope_internal_facets(source(same_winding=True))

if __name__=='__main__':unittest.main()
