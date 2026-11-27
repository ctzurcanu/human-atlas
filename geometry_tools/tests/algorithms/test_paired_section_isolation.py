import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.algorithms.paired_section_isolation import isolate_paired_section


def returning_pipe():
    """Two vertical stems and a rounded return, with distinct wall sheets."""
    a=np.c_[np.zeros(5),np.zeros(5),np.linspace(0,2,5)]
    theta=np.linspace(0,np.pi,17)[1:]
    bend=np.c_[1-np.cos(theta),np.zeros(len(theta)),2+np.sin(theta)]
    b=np.c_[np.full(4,2.),np.zeros(4),np.linspace(1.5,0,4)]
    centers=np.r_[a,bend,b]
    tangents=np.r_[np.tile([0.,0.,1.],(5,1)),
        np.c_[np.sin(theta),np.zeros(len(theta)),np.cos(theta)],np.tile([0.,0.,-1.],(4,1))]
    right=np.c_[tangents[:,2],np.zeros(len(tangents)),-tangents[:,0]]
    turn=np.arange(8)*2*np.pi/8
    ring=np.cos(turn)[None,:,None]*right[:,None]+np.sin(turn)[None,:,None]*np.array([0.,1.,0.])
    v=np.r_[(centers[:,None]+.1*ring).reshape(-1,3),(centers[:,None]+.08*ring).reshape(-1,3)]
    total=len(centers)*8;outer=[]
    for k in range(len(centers)-1):
        for j in range(8):
            x=k*8+j;y=k*8+(j+1)%8
            outer.extend([[x,y,y+8],[x,y+8,x+8]])
    f=outer+[[x+total,z+total,y+total]for x,y,z in outer]
    cap_ids=[]
    for k,invert in [(0,False),(len(centers)-1,True)]:
        here=[]
        for j in range(8):
            x=k*8+j;y=k*8+(j+1)%8
            faces=[[x,x+total,y+total],[x,y+total,y]]
            if invert:faces=[q[::-1]for q in faces]
            here.extend(range(len(f),len(f)+2));f.extend(faces)
        cap_ids.append(np.array(here))
    return v,np.array(f),cap_ids


class PairedSectionIsolationTests(unittest.TestCase):
    def test_a_return_beyond_the_plane_stays_with_its_own_source_branch(self):
        v,f,caps=returning_pipe();r=isolate_paired_section(v,f,[0,0,.75],[0,0,1])
        labels=r['component_labels'];a=np.unique(labels[np.isin(r['source_face_ids'],caps[0])])
        b=np.unique(labels[np.isin(r['source_face_ids'],caps[1])])
        self.assertEqual(len(a),1);self.assertEqual(len(b),1);self.assertNotEqual(a[0],b[0])
        moving=r['faces'][labels==b[0]];used=np.unique(moving)
        self.assertTrue(np.any((r['positions'][used,0]>1.9)&(r['positions'][used,2]<.1)))
        self.assertEqual(r['other_plane_contours_reconnected'],2)
        self.assertTrue(topology(Mesh(r['positions'],r['faces'],'fixture:metres'))['closed_oriented_manifold'])
        np.testing.assert_array_equal(r['positions'][:len(v)],v)
        old=v[f];expected=np.linalg.norm(np.cross(old[:,1]-old[:,0],old[:,2]-old[:,0]),axis=1)/2
        p=r['positions'][r['faces']];areas=np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)/2
        actual=np.bincount(r['source_face_ids'],weights=areas,minlength=len(f))
        np.testing.assert_allclose(actual,expected,rtol=3e-14,atol=1e-17)
        self.assertEqual(r['anatomical_caps_added'],0)
    def test_a_query_outside_the_wall_pair_is_rejected(self):
        v,f,_=returning_pipe()
        with self.assertRaisesRegex(ValueError,'exactly two wall contours'):
            isolate_paired_section(v,f,[3,0,.75],[0,0,1])
    def test_unrelated_returning_faces_and_interfaces_retain_literal_triangles(self):
        v,f,caps=returning_pipe();r=isolate_paired_section(v,f,[0,0,.75],[0,0,1],restore_unrelated=True)
        keys={tuple(face)for face in r['faces']}
        returning=np.all(v[f][:,:,0]>1.8,axis=1)
        self.assertTrue(all(tuple(face)in keys for face in f[returning]))
        self.assertGreater(r['unrelated_subdivided_parents_restored'],0)
        self.assertTrue(topology(Mesh(r['positions'],r['faces'],'fixture:metres'))['closed_oriented_manifold'])
        self.assertEqual(set(r['component_labels']),{0,1})
        self.assertTrue(all(tuple(f[i])in keys for cap in caps for i in cap))


if __name__=='__main__':unittest.main()
