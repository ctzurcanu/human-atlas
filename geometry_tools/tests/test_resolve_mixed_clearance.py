import unittest
import numpy as np
from geometry_tools.algorithms.rigid_section_clearance import fit_section_clearance
from geometry_tools.contact import self_contact_summary
from geometry_tools.model import Mesh

class MixedContactTests(unittest.TestCase):
    def test_existing_mixed_crossing_is_explicitly_resolved(self):
        v=np.array([[-.001,0,0],[.001,0,0],[0,.001,0],
                    [0,.0002,-.0003],[0,.0002,.0003],[0,-.0005,0]])
        f=np.array([[0,1,2],[3,4,5]])
        options=dict(reference_vertices=v,face_system=np.array([0,1]),gap_m=.00005,
                     step_m=.001,maximum_steps=2,wall_seconds=3.)
        old,proof=fit_section_clearance(v,f,[],[],**options)
        np.testing.assert_array_equal(old,v)
        self.assertEqual(proof['originalSourceContactPairsExplicitlyRetained'],1)
        fixed,proof=fit_section_clearance(v,f,[],[],resolve_mixed_contacts=True,**options)
        check,_=self_contact_summary(Mesh(fixed,f,'test:mixed'))
        self.assertEqual(check['defect_pairs'],0)
        self.assertEqual(proof['originalSourceContactPairsExplicitlyRetained'],0)
        self.assertTrue(proof['resolveOriginalArterialVenousContacts'])

    def test_coherent_escape_retains_profiles_and_fixed_vein(self):
        theta=np.arange(6)*np.pi/3
        def tube(x):
            v=np.concatenate([np.c_[x+.0003*np.cos(theta),.0003*np.sin(theta),np.full(6,z)]for z in [0,.004]])
            f=np.array([[k,(k+1)%6,(k+1)%6+6]for k in range(6)]+[[k,(k+1)%6+6,k+6]for k in range(6)])
            return v,f
        a,fa=tube(-.00015);b,fb=tube(.00015);v=np.r_[a,b];f=np.r_[fa,fb+len(a)]
        sections=[np.arange(k,k+6)for k in range(0,len(v),6)]
        q,proof=fit_section_clearance(v,f,sections,[],protected=np.arange(len(a),len(v)),reference_vertices=v,face_system=np.r_[np.zeros(len(fa),int),np.ones(len(fb),int)],gap_m=.00005,step_m=.0008,maximum_steps=2,wall_seconds=3.,resolve_mixed_contacts=True,coherent_mixed_sections=True)
        np.testing.assert_array_equal(q[len(a):],v[len(a):])
        check,_=self_contact_summary(Mesh(q,f,'test:coherent'))
        self.assertEqual(check['defect_pairs'],0)
        for ring in sections:
            np.testing.assert_allclose(q[ring]-q[ring].mean(0),v[ring]-v[ring].mean(0),atol=1e-15)

if __name__=='__main__':unittest.main()
