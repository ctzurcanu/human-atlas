import unittest
from unittest.mock import patch
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.constraints import check, check_constraint
from geometry_tools.install import changed_constraint_meshes


class InstallAttributeConstraints(unittest.TestCase):
    def setUp(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        f=np.array([[0,1,2]])
        n=np.tile(np.array([0,0,32767],np.int16),(3,1))
        self.first=Mesh(p,f,'fixture:metres',{'normals':n})
        self.second=Mesh(p,f[:,::-1],'fixture:metres',{'normals':-n})
        self.row=dict(id='shared',kind='shared-interface',first_mesh_id='a',second_mesh_id='b',
                      parameters=dict(first_faces=[0],second_faces=[0]))
        self.registrations={pid:dict(id=pid,geometry_sha256=m.geometry_sha,
                                    attributes_sha256=m.attributes_sha)
                            for pid,m in [('a',self.first),('b',self.second)]}

    def run_check(self,candidates):
        def historical(db,mid):
            return None,{'a':self.first,'b':self.second}[mid]
        with patch('geometry_tools.constraints.list_constraints',return_value=[self.row]), \
             patch('geometry_tools.constraints.load_registered',side_effect=historical):
            return check(None,candidates=candidates)

    def changed_pair(self):
        words=np.tile(np.array([10,0,32767],np.int16),(3,1))
        a=Mesh(self.first.vertices,self.first.faces,self.first.frame,{'normals':words})
        b=Mesh(self.second.vertices,self.second.faces,self.second.frame,{'normals':-words})
        return a,b

    def test_corresponding_attribute_only_edit_uses_both_current_endpoints(self):
        a,b=self.changed_pair()
        self.assertEqual(a.geometry_sha,self.first.geometry_sha)
        self.assertFalse(self.run_check({'a':a})['passed'])
        candidates=changed_constraint_meshes(self.registrations,{'a':a,'b':b})
        self.assertEqual(set(candidates),{'a','b'})
        self.assertTrue(self.run_check(candidates)['passed'])

    def test_unmatched_new_normal_words_are_still_rejected(self):
        a,_=self.changed_pair()
        candidates=changed_constraint_meshes(self.registrations,{'a':a,'b':self.second})
        result=self.run_check(candidates)
        self.assertFalse(result['passed'])
        self.assertTrue(result['constraints'][0]['literal_opposed_triangles'])
        self.assertFalse(result['constraints'][0]['literal_opposed_normal_words'])

    def test_unchanged_endpoints_and_fixed_geometry_semantics_remain(self):
        self.assertEqual(changed_constraint_meshes(self.registrations,
                          {'a':self.first,'b':self.second}),{})
        a,_=self.changed_pair()
        fixed=dict(kind='fixed',parameters={'geometry_sha256':self.first.geometry_sha})
        self.assertTrue(check_constraint(fixed,a)['passed'])
        moved=Mesh(a.vertices+np.array([0.,0.,.001]),a.faces,a.frame,a.attributes)
        candidates=changed_constraint_meshes({'a':self.registrations['a']},{'a':moved})
        self.assertIn('a',candidates)
        self.assertFalse(check_constraint(fixed,candidates['a'])['passed'])


if __name__=='__main__':unittest.main()
