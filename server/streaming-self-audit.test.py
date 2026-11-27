import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts,self_contact_summary
from geometry_tools.spatial import BVH

class StreamingSelfAuditTest(unittest.TestCase):
    def source(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.2,.2,-.2],[.2,.2,.7],[.8,.2,.2]],float)
        return Mesh(v,np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3],[4,5,6]]),'test:metres')
    def test_stream_matches_full_witness_audit_and_budget_edges(self):
        m=self.source()
        for pairs,tests in [(1000,1000),(1,1000),(1000,1)]:
            a,s=pair_contacts(m,m,0.,selfcheck=True,max_pairs=pairs,max_tests=tests);t,_=self_contact_summary(m,max_pairs=pairs,max_tests=tests)
            for key in s:self.assertEqual(t[key],s[key],key)
            self.assertEqual(t['defect_pairs'],int(np.isin(a['relation'],['crossing','coplanar-overlap']).sum()))
    def test_disjoint_bvh_partition_is_complete(self):
        m=self.source();b=BVH(m,leaf_size=2);left,right=map(int,b.children[0]);whole,_=self_contact_summary(m,bvh=b,max_pairs=1000,max_tests=1000)
        pieces=[self_contact_summary(m,bvh=b,root_pair=p,max_pairs=1000,max_tests=1000)[0]for p in [(left,left),(left,right),(right,right)]]
        for key in ['candidate_triangle_pairs_tested','returned_triangle_pairs','crossing_face_pairs','coplanar_overlap_face_pairs','touching_face_pairs','defect_pairs']:
            self.assertEqual(sum(p[key]for p in pieces),whole[key],key)

if __name__=='__main__':unittest.main()
