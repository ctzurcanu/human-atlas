import unittest
from unittest.mock import patch
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.interface_transfer import verify_transfer
from geometry_tools.constraints import check,_literal_cyclic_equal


def pair():
    p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
    f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
    a=Mesh(p,f,'fixture:metres',{'normals':np.tile([0,0,-32767],(4,1)).astype(np.int16)})
    q=p.copy();q[3,2]=-1
    b=Mesh(q,f[:,[0,2,1]],a.frame,{'normals':-a.attributes['normals']})
    return a,b,dict(kind='shared-interface',parameters={'first_faces':[0],'second_faces':[0]})


class InterfaceTransfer(unittest.TestCase):
    def test_literal_multisets_preserve_orientation_multiplicity_and_integer_bits(self):
        position=np.array([[[0,0,0],[1,0,0],[0,1,0]],[[0,0,0],[1,0,0],[0,1,0]]],float)
        words=np.full((2,3,3),2**60,dtype=np.int64);words[1,0,0]+=1
        fields=[position[:,:,k]for k in range(3)]+[words[:,:,k]for k in range(3)]
        reordered=[np.roll(f[::-1],1,axis=1)for f in fields]
        self.assertTrue(_literal_cyclic_equal(fields,reordered))
        changed=[f.copy()for f in reordered];changed[3][0,1]-=1
        self.assertFalse(_literal_cyclic_equal(fields,changed))
        reversed_fields=[f[:,[0,2,1]]for f in fields]
        self.assertFalse(_literal_cyclic_equal(fields,reversed_fields))
        self.assertFalse(_literal_cyclic_equal(fields,[f[:1]for f in fields]))

    def test_new_owners_can_reorder_only_the_literal_source_sheet(self):
        a,b,old=pair();na=Mesh(a.vertices,a.faces[::-1],a.frame,a.attributes)
        proof=verify_transfer(old,a,b,na,b,{'first_faces':[3],'second_faces':[0]})
        self.assertTrue(proof['passed']);self.assertTrue(proof['historical_constraint_unchanged'])

    def test_moved_seam_or_changed_normal_words_are_rejected(self):
        a,b,old=pair();q=a.vertices.copy();q[0,0]+=1e-12
        with self.assertRaisesRegex(ValueError,'every oriented source facet'):
            verify_transfer(old,a,b,Mesh(q,a.faces,a.frame,a.attributes),b,old['parameters'])
        words=a.attributes['normals'].copy();words[0,2]+=1
        with self.assertRaisesRegex(ValueError,'every oriented source facet'):
            verify_transfer(old,a,b,Mesh(a.vertices,a.faces,a.frame,{'normals':words}),b,old['parameters'])

    def test_unverified_or_duplicate_source_faces_are_rejected(self):
        a,b,old=pair()
        with self.assertRaisesRegex(ValueError,'verified shared tissue'):
            verify_transfer(dict(old,parameters=dict(old['parameters'],verification_pending=True)),a,b,a,b,old['parameters'])
        with self.assertRaisesRegex(ValueError,'nonduplicate'):
            verify_transfer(old,a,b,a,b,dict(first_faces=[0,0],second_faces=[0]))

    def test_fixed_constraints_cannot_be_transferred(self):
        a,b,_=pair();row=dict(id='fixed',kind='fixed',first_mesh_id='a',second_mesh_id=None,parameters={})
        with patch('geometry_tools.constraints.list_constraints',return_value=[row]):
            report=check(None,candidates={'a':a},endpoint_overrides={'fixed':(a,b)})
        self.assertFalse(report['passed'])

    def test_current_candidate_does_not_reload_the_discarded_old_buffer(self):
        a,_,_=pair();row=dict(id='fixed',kind='fixed',first_mesh_id='a',second_mesh_id=None,parameters={'geometry_sha256':a.geometry_sha})
        with patch('geometry_tools.constraints.list_constraints',return_value=[row]),patch('geometry_tools.constraints.load_registered',side_effect=AssertionError('Unnecessary historical reload')):
            report=check(None,candidates={'a':a})
        self.assertTrue(report['passed'])


if __name__=='__main__':unittest.main()
