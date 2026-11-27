import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.profile_support_separation import split_contact_roles,support_rows,minimum_support_gap
from geometry_tools.distance import triangle_pairs


class ProfileSupportTests(unittest.TestCase):
    def test_artery_vein_crossings_are_never_conduit_coupling(self):
        p=[dict(system=x)for x in ['arterial','arterial','venous']];contacts=[dict(first=0,second=1,pairs=[[0,0]]),dict(first=0,second=2,pairs=[[0,0]])];kept,separate=split_contact_roles(p,contacts);self.assertEqual(kept,[contacts[0]]);self.assertEqual(separate,[contacts[1]])

    def test_support_inequality_separates_complete_triangles(self):
        a=np.array([[0.,0,0],[1,0,0],[0,1,0]]);b=np.array([[.25,.25,-1],[.25,.25,1],[.75,.25,0]]);meshes=[Mesh(v,np.array([[0,1,2]]),'fixture:metres')for v in [a,b]];groups=[[np.arange(3)]]*2;routes=[dict(artery=0,vein=1,direction=[1.,0,0],groupPairs=[[0,0]])];W,N,B,_=support_rows(meshes,groups,[np.array([0]),np.array([1])],routes,clearance_m=.01,columns=2);ds=np.array([[B[0],0,0],[0,0,0]]);self.assertGreaterEqual(float(np.sum(N*(W@ds))),B[0]-1e-12);moved=[Mesh(m.vertices+ds[i],m.faces,m.frame)for i,m in enumerate(meshes)];self.assertGreaterEqual(minimum_support_gap(moved,groups,routes),.01-1e-12);result=triangle_pairs(moved[0].triangles,moved[1].triangles);self.assertGreaterEqual(float(result['distance_m'][0]),.01-1e-12);self.assertNotIn(result['relation'][0],['crossing','coplanar-overlap'])


if __name__=='__main__':unittest.main()
