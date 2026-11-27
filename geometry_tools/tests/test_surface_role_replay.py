import unittest
import numpy as np
from geometry_tools.algorithms.surface_role_replay import split_by_source_role


class SurfaceRoleReplayTests(unittest.TestCase):
    def test_split_preserves_both_owners_and_literal_geometry(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[1,1,0]])
        f=np.array([[0,1,2],[1,3,2]])
        r=dict(positions=p,surfaces={'shared':f},cells={
            'a':[dict(surface='shared',side=1,neighbor='b')],
            'b':[dict(surface='shared',side=-1,neighbor='a')]},charts={
            'shared':dict(source_surfaces=['outer','inner'],source_faces=np.array([5,7]),
                normal_signs=np.array([1,-1]),barycentrics=np.repeat(np.eye(3)[None],2,axis=0))})
        out,roles=split_by_source_role(r,{'inner':'explicit-inner'})
        self.assertIs(out['positions'],p)
        self.assertEqual(set(roles.values()),{'','explicit-inner'})
        np.testing.assert_array_equal(np.concatenate(list(out['surfaces'].values())),f)
        for a,b in zip(out['cells']['a'],out['cells']['b']):
            self.assertEqual(a['surface'],b['surface']);self.assertEqual(a['side'],-b['side'])
        inner=next(k for k,v in roles.items() if v)
        self.assertEqual(out['charts'][inner]['source_faces'].tolist(),[7])
        self.assertEqual(out['charts'][inner]['normal_signs'].tolist(),[-1])


if __name__=='__main__':unittest.main()
