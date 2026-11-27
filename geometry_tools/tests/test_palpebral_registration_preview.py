"""The registration case is explicit and cannot erase the strict material rejection."""
import copy,importlib.util,json,tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha
from geometry_tools.tests.test_palpebral_install_verifier import header,source_tube,source_row,moved

spec=importlib.util.spec_from_file_location('palpebral_registration_tests',
 Path(__file__).resolve().parents[2]/'scripts/verify-palpebral-registration-preview.py')
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)


def preview_header():
 c,r=header();c.update(courseCase=h.CASE,physicsAssembly=h.strict.EXPECTED_OWNERS,
  movingOwners=h.MOVING,literalRecipientOwners=h.FIXED,materialMetricAccepted=False,
  originalMaterialGateUnchangedAndRejected=True,maximumAnalyticSectionCorrection_m=.00025,
  registrationSourceReconciliation_m=0.)
 r.update(materialMetricAccepted=False,originalMaterialGateUnchangedAndRejected=True)
 return c,r


class RegistrationPreviewTests(unittest.TestCase):
 def test_explicit_separate_identity_retains_original_material_rejection(self):
  c,r=preview_header();h.verify_header(c,r)
  with self.assertRaises(AssertionError):h.strict.verify_header(c,r)
  for key,value in [('materialMetricAccepted',True),('originalMaterialGateUnchangedAndRejected',False),
    ('physicsAssembly',h.strict.EXPECTED_OWNERS[:-1]),('literalRecipientOwners',h.FIXED[:-1]),
    ('maximumFacetTurnDegrees',66.),('maximumAnalyticSectionCorrection_m',.0003)]:
   bad=copy.deepcopy(c);bad[key]=value
   with self.subTest(key=key),self.assertRaises(AssertionError):h.verify_header(bad,r)

 def test_axial_registration_measurement_does_not_claim_material_or_length_acceptance(self):
  old=source_tube();q=old.vertices.copy();q[:,2]*=2
  new=moved(old,q);c,_=preview_header()
  proof=h.verify_registration_mesh(old,new,source_row(),c)
  self.assertAlmostEqual(proof['metrics']['maximumSingularScale'],2.)
  self.assertTrue(proof['donorMaterialStrainDiagnosticOnly'])
  self.assertFalse(proof['donorLengthOrElasticityAcceptance'])
  with self.assertRaises(AssertionError):h.strict.verify_mesh_geometry(old,new,source_row(),c)

 def test_registration_preserves_calibre_normal_motion_and_fixed_stars(self):
  old=source_tube();c,_=preview_header()
  q=old.vertices.copy();q[:,:2]*=1.09
  with self.assertRaises(AssertionError):h.verify_registration_mesh(old,moved(old,q),source_row(),c)
  row=source_row();row['fixedStarRawIds']=[0]
  with self.assertRaises(AssertionError):h.verify_registration_mesh(old,moved(old,old.vertices+[1e-12,0,0]),row,c)
  with self.assertRaises(AssertionError):h.verify_registration_mesh(old,moved(old,old.vertices+[.013,0,0]),source_row(),c)

 def test_literal_recipients_reject_even_a_tiny_position_or_appearance_change(self):
  old=source_tube();self.assertTrue(h.verify_literal_recipient(old,old))
  with self.assertRaises(AssertionError):h.verify_literal_recipient(old,moved(old,old.vertices+[1e-13,0,0]))
  attributes={k:v.copy()for k,v in old.attributes.items()};attributes['uv'][0,0]+=1
  with self.assertRaises(AssertionError):h.verify_literal_recipient(old,moved(old,old.vertices,attributes=attributes))

 def test_four_packed_rows_require_six_separate_literal_context_rows(self):
  c,r=preview_header();c.update(installWrapper=True,wholeAssembly=h.MOVING,
   fixedContextRows=[dict(partId=p)for p in h.FIXED],registrationEvidencePending=False,registrationEvidenceVerified=True)
  r.update(rows=[dict(partId=p)for p in h.MOVING],passedGeometryGates=True,
   passedRegistrationPreviewCriteria=True,newIntrinsicPairs=0,registrationReviewPending=False)
  h.verify_header(c,r)
  c['fixedContextRows']=c['fixedContextRows'][:-1]
  with self.assertRaises(AssertionError):h.verify_header(c,r)

 def test_complete_hash_bound_independent_review_is_required(self):
  with tempfile.TemporaryDirectory()as folder:
   root=Path(folder);cp=root/'contract.json';rp=root/'report.json';ep=root/'repaired-strip-audit.json';mp=root/'registration-review.json'
   cp.write_text('{}');rp.write_text('{}');ep.write_text('{"actualAnnularJacobianPositive":true}')
   manifest=dict(complete=True,passed=True,sourceContractSha256=sha(cp),sourceReportSha256=sha(rp),
    evidence={str(ep):sha(ep)},**{key:True for key in h.REQUIRED_REVIEW})
   mp.write_text(json.dumps(manifest));h.verify_review_manifest(mp,cp,rp)
   for key in ['passed','actualAnnularGeometryPassed','completeContextPassed','noUnresolvedNewPhysicalCrossings','fullScopeRenderReviewed']:
    bad=copy.deepcopy(manifest);bad[key]=False;mp.write_text(json.dumps(bad))
    with self.subTest(key=key),self.assertRaises(AssertionError):h.verify_review_manifest(mp,cp,rp)
   mp.write_text(json.dumps(manifest));ep.write_text('{"actualAnnularJacobianPositive":false}')
   with self.assertRaises(AssertionError):h.verify_review_manifest(mp,cp,rp)

 def test_each_actual_repair_obeys_its_declared_direction_and_bound(self):
  old=source_tube().vertices;inferior=next(p for p in h.MOVING if '_inferior_'in p);superior=next(p for p in h.MOVING if '_superior_'in p)
  h.verify_actual_repair_delta(inferior,old,old+[.0002,0,0])
  h.verify_actual_repair_delta(superior,old,old+[0,0,.001])
  for pid,shift in [(inferior,[.000251,0,0]),(superior,[0,0,.001201]),
    (superior,[0,0,-1e-10]),(superior,[1e-12,0,.0001])]:
   with self.subTest(pid=pid,shift=shift),self.assertRaises(AssertionError):
    h.verify_actual_repair_delta(pid,old,old+shift)


if __name__=='__main__':unittest.main()
