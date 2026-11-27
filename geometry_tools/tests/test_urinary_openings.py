import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.urinary_openings import bind_opening_boundary, shared_lumen_rim


class OpeningWitnessTests(unittest.TestCase):
    def fixture(self):
        xy=np.array([[-2.,-2.],[2.,-2.],[2.,2.],[-2.,2.],[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])
        v=np.r_[np.c_[xy,np.full(8,.1)],np.c_[xy,np.full(8,-.1)]]
        f=[]
        for k in range(4):
            j=(k+1)%4;f.extend([[k,j,4+j],[k,4+j,4+k]])
        top=np.array(f);f.extend((top[:,[0,2,1]]+8).tolist())
        for k in range(4):
            j=(k+1)%4;f.extend([[k,k+8,j+8],[k,j+8,j],[k+4,j+4,j+12],[k+4,j+12,k+12]])
        mesh=Mesh(v,np.array(f),'test:metres')
        record=dict(tissueDrawGroups=[dict(start=0,count=24),dict(start=24,count=(len(f)-8)*3,neighbor='other')])
        return v[4:8],{'mucosa':(mesh,record)}

    def test_witnesses_follow_the_actual_annular_wall_without_a_lumen_disk(self):
        curve,supports=self.fixture();r=bind_opening_boundary('neck',curve,supports,'a'*64)
        self.assertEqual(r['representation'],'wall-rim');self.assertEqual(len(r['points']),16)
        self.assertLess(r['maximumSourceOffsetM'],1e-15)
        m=supports['mucosa'][0]
        q=np.einsum('ij,ijk->ik',r['barycentric'],m.triangles[r['faceIds']])
        np.testing.assert_allclose(q,r['points'],atol=1e-15)

    def test_a_generic_nearby_wall_cannot_establish_an_opening(self):
        curve,supports=self.fixture();curve[:,:2]*=.1
        with self.assertRaises(ValueError):bind_opening_boundary('neck',curve,supports,'a'*64)


class SharedLumenRimTests(unittest.TestCase):
    def pair(self):
        angles=np.arange(8)*np.pi/4
        def cell(y0,y1,pid,peer,seam_top):
            rings=[np.c_[r*np.cos(angles),np.full(8,y),r*np.sin(angles)]
                   for y,r in ((y0,1.),(y0,.5),(y1,1.),(y1,.5))]
            vertices=np.concatenate(rings);faces=[];groups=[]
            for name,first,second,flip in [('outer',0,16,False),('inner',8,24,True),('bottom',0,8,True),('top',16,24,False)]:
                start=len(faces)*3
                for k in range(8):
                    j=(k+1)%8;pair=[[first+k,first+j,second+j],[first+k,second+j,second+k]]
                    faces.extend([t[::-1] if flip else t for t in pair])
                group=dict(start=start,count=len(faces)*3-start)
                if name!='inner':group['neighbor']=peer if name==('top' if seam_top else 'bottom') else 'other'
                groups.append(group)
            return Mesh(vertices,np.array(faces),'fixture:metres'),dict(id=pid,tissueDrawGroups=groups)
        return cell(0.,1.,'first','second',True)+cell(1.,2.,'second','first',False)

    def test_literal_join_yields_only_the_actual_inner_lumen_loop(self):
        a,ra,b,rb=self.pair();r=shared_lumen_rim(a,ra,b,rb);p=np.array(r['points'])
        self.assertEqual(r['edge_count'],8);self.assertEqual(r['shared_tissue_triangles'],16)
        np.testing.assert_allclose(p[:,1],1.)
        np.testing.assert_allclose(np.linalg.norm(p[:,[0,2]],axis=1),.5)
        self.assertEqual(r['new_geometry_bytes'],0)
        self.assertFalse(r['opposed_normal_words_verified'])

    def test_a_micron_gap_is_not_silently_joined(self):
        a,ra,b,rb=self.pair();b=Mesh(b.vertices+np.array([0.,1e-6,0.]),b.faces,b.frame)
        with self.assertRaisesRegex(ValueError,'literal opposed'):shared_lumen_rim(a,ra,b,rb)

    def test_reversed_tissue_winding_is_rejected(self):
        a,ra,b,rb=self.pair();b.faces=b.faces[:,[0,2,1]]
        with self.assertRaisesRegex(ValueError,'literal opposed'):shared_lumen_rim(a,ra,b,rb)

    def test_mismatched_render_normal_words_are_rejected(self):
        a,ra,b,rb=self.pair()
        for mesh in (a,b):mesh.attributes['normals']=np.tile(np.array([0,32767,0],np.int16),(len(mesh.vertices),1))
        with self.assertRaisesRegex(ValueError,'opposed normal'):shared_lumen_rim(a,ra,b,rb)

    def test_outer_and_inner_free_rims_are_not_confused(self):
        a,ra,b,rb=self.pair()
        for r in (ra,rb):r['tissueDrawGroups'][0].pop('neighbor')
        with self.assertRaisesRegex(ValueError,'More than one'):shared_lumen_rim(a,ra,b,rb)


if __name__=='__main__':unittest.main()
