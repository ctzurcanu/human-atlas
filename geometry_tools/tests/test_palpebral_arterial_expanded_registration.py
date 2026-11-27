"""Arterial policy is separate, source-bound and cannot invent medial joins."""
import copy,importlib.util,unittest
from pathlib import Path
import numpy as np
from geometry_tools.tests.test_palpebral_install_verifier import source_tube,source_row,moved
spec=importlib.util.spec_from_file_location('palpebral_arterial_expanded_test',Path(__file__).resolve().parents[2]/'scripts/verify-palpebral-arterial-expanded-registration.py')
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)


def header():
 c=dict(complete=True,courseCase=h.CASE,sourceAtlasSha256=h.SOURCE_SHA,wholeAssembly=h.OWNERS,
  interfaces=[{}for _ in range(6)],externalInterfaces=[{}for _ in range(6)],pendingOffsetPairs=[],nonelasticRegistration=True,scopeVersion=2,**h.LIMITS)
 r=dict(complete=True,sourceAtlasSha256=h.SOURCE_SHA,rows=[dict(partId=p)for p in h.OWNERS])
 return c,r


class ExpandedArterialRegistrationTests(unittest.TestCase):
 def test_exact_ten_owner_source_and_predeclared_limits(self):
  c,r=header();self.assertEqual(len(h.OWNERS),10);h.verify_header(c,r)
  for key,value in [('scopeVersion',1),('sourceAtlasSha256','old-source'),('maximumPointMotion_m',.0151),
    ('maximumFacetTurnDegrees',66),('wholeAssembly',h.OWNERS[:-1]),
    ('pendingOffsetPairs',[dict(first='medial',second='arch')]),('donorElasticityClaim',True)]:
   bad=copy.deepcopy(c);bad[key]=value
   with self.subTest(key=key),self.assertRaises(AssertionError):h.verify_header(bad,r)

 def test_actual_15mm_motion_boundary_and_source_calibre(self):
  c,_=header();old=source_tube()
  h.verify_mesh(old,moved(old,old.vertices+[.014,0,0]),source_row(),c)
  with self.assertRaises(AssertionError):h.verify_mesh(old,moved(old,old.vertices+[.016,0,0]),source_row(),c)
  q=old.vertices.copy();q[:,:2]*=1.09
  with self.assertRaises(AssertionError):h.verify_mesh(old,moved(old,q),source_row(),c)

 def test_nonelastic_source_length_measurement_is_explicit_not_material_acceptance(self):
  c,_=header();old=source_tube();q=old.vertices.copy();q[:,2]*=2
  proof=h.verify_mesh(old,moved(old,q),source_row(),c)
  self.assertAlmostEqual(proof['metrics']['maximumSingularScale'],2.)
  self.assertTrue(proof['donorMaterialStrainDiagnosticOnly'])
  self.assertFalse(proof['finiteWallOrPhysiologicalPatencyAccepted'])

 def test_literal_source_owners_are_context_only_in_a_moving_subset(self):
  c,r=header();moving=[p for p in h.OWNERS if '_arch_ID'in p]
  c.update(installWrapper=True,physicsAssembly=h.OWNERS,movingOwners=moving,wholeAssembly=moving,
   fixedContextRows=[dict(partId=p)for p in h.OWNERS if p not in moving])
  r.update(rows=[dict(partId=p)for p in moving],passedGeometryGates=True,newIntrinsicPairs=0)
  h.verify_header(c,r)
  c['fixedContextRows']=c['fixedContextRows'][:-1]
  with self.assertRaises(AssertionError):h.verify_header(c,r)


if __name__=='__main__':unittest.main()
