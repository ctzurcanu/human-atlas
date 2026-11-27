"""Focused finite guards; no candidate packing, geometry edits or installation."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'server'))
from geometry_tools.model import Mesh
from geometry_tools.constraints import check_constraint
from tongue_basal_contract import (FRAME, MUCO, bind_features, corner_normal_mesh,
                                  covering_muscle_gate, retained_envelope_gate,selected_extrinsics_gate)


def pair_fixture():
    covers = {key:dict(path='/tmp/'+key, sha256=key) for key in ('mucosa','aponeurosis')}
    muscles = {str(i):dict(path='/tmp/muscle-'+str(i),sha256=str(i)) for i in range(10)}
    pairs = [dict(cover=a,muscle=b,cover_path=covers[a]['path'],cover_sha256=covers[a]['sha256'],
        muscle_path=muscles[b]['path'],muscle_sha256=muscles[b]['sha256'],
        witnesses='/tmp/witness.npz',witnesses_sha256='x',summary=dict(
            complete_within_pair_budget=True,truncated=False,crossing_face_pairs=0,
            coplanar_overlap_face_pairs=0,degenerate_pair_count=0)) for a in covers for b in muscles]
    return dict(schema='tongue-finite-covering-muscle-audit/v1',passed=True,pairs=pairs),covers,muscles


class ContractTests(unittest.TestCase):
    def test_complete_twenty_zero_pairs_pass(self):
        audit,covers,muscles=pair_fixture()
        self.assertEqual(covering_muscle_gate(audit,covers,muscles,check_files=False)['pairs'],20)

    def test_ignored_summary_pass_cannot_hide_crossings(self):
        audit,covers,muscles=pair_fixture()
        audit['pairs'][6]['summary']['crossing_face_pairs']=1
        with self.assertRaisesRegex(ValueError,'crossings block packing'):
            covering_muscle_gate(audit,covers,muscles,check_files=False)

    def test_duplicate_does_not_replace_missing_pair(self):
        audit,covers,muscles=pair_fixture()
        audit['pairs'][-1]=copy.deepcopy(audit['pairs'][0])
        with self.assertRaisesRegex(ValueError,'twenty identity-distinct'):
            covering_muscle_gate(audit,covers,muscles,check_files=False)

    def test_stale_repaired_gg_rejected(self):
        audit,covers,muscles=pair_fixture()
        audit['pairs'][0]['muscle_sha256']='old-native-gg'
        with self.assertRaisesRegex(ValueError,'another physical source'):
            covering_muscle_gate(audit,covers,muscles,check_files=False)

    def test_no_undeclared_coplanar_exception(self):
        audit,covers,muscles=pair_fixture()
        audit['pairs'][0]['summary']['coplanar_overlap_face_pairs']=1
        with self.assertRaisesRegex(ValueError,'coplanar occupied surfaces'):
            covering_muscle_gate(audit,covers,muscles,check_files=False)

    def test_literal_seam_requires_all_face_witnesses(self):
        audit,covers,muscles=pair_fixture()
        pair=audit['pairs'][0]
        pair['summary']['coplanar_overlap_face_pairs']=1
        seam={(pair['cover'],pair['muscle']):dict(cover_faces=[1],muscle_faces=[2],
            actual_check=dict(passed=True,literal_opposed_triangles=True,literal_opposed_normal_words=True))}
        with self.assertRaisesRegex(ValueError,'explicit source-face witnesses'):
            covering_muscle_gate(audit,covers,muscles,check_files=False,shared_interfaces=seam)
        pair['coplanar_face_pairs']=[[1,2]]
        self.assertTrue(covering_muscle_gate(audit,covers,muscles,check_files=False,shared_interfaces=seam)['passed'])
        pair['coplanar_face_pairs']=[[1,3]]
        with self.assertRaisesRegex(ValueError,'outside the literal declared seam'):
            covering_muscle_gate(audit,covers,muscles,check_files=False,shared_interfaces=seam)

    def test_actual_current_shell_failure_is_retained(self):
        path=ROOT/'.local-models/geometry-workbench/tongue-completion/corrected-covering-interface-audit-v1/summary.json'
        if not path.exists():
            self.skipTest('Local diagnostic is not available in this checkout')
        audit=json.loads(path.read_text())
        covers={r['cover']:dict(path=r['cover_path'],sha256=r['cover_sha256']) for r in audit['pairs']}
        muscles={r['muscle']:dict(path=r['muscle_path'],sha256=r['muscle_sha256']) for r in audit['pairs']}
        self.assertFalse(audit['passed'])
        with self.assertRaisesRegex(ValueError,'crossings block packing'):
            covering_muscle_gate(audit,covers,muscles)

    def test_actual_selected_corner_fields_do_not_change_geometry(self):
        directory=ROOT/'.local-models/geometry-workbench/tongue-completion/basal-macro-envelopes-v1-corner-charts'
        if not directory.exists():
            self.skipTest('Local corner chart candidate is not available')
        construction=json.loads((directory/'construction.json').read_text())
        for row in construction['meshes']:
            actual=corner_normal_mesh(row['mesh'],row['geometry_sha256'])
            original=Mesh.load(row['mesh'])
            self.assertTrue(np.array_equal(actual.triangles,original.triangles))
            self.assertEqual(actual.geometry_sha,original.geometry_sha)

    def test_partial_exposure_requires_features_to_be_regenerated(self):
        vertices=np.array([[0.,1.6,0.],[.01,1.6,0.],[0.,1.61,0.],[.01,1.61,0.]])
        source=Mesh(vertices,np.array([[0,1,2],[1,3,2]]),FRAME)
        carrier=Mesh(vertices,np.array([[0,1,2]]),FRAME)
        ancestry=dict(schema='tongue-covering-face-ancestry/v1',carrierPartId=MUCO,
            carrierBankSha256='a'*64,sourceGeometrySha256=source.geometry_sha,
            carrierGeometrySha256=carrier.geometry_sha,
            faces=[dict(sourceFace=0,carrierFace=0,cornerOrder=[0,1,2])])
        feature=dict(partIds=[MUCO],surface=dict(partId=MUCO,faces=[1]))
        features=dict(schema='tongue-construction-features/v1',hostPartId=MUCO,
            additionalSolidMeshes=0,terms=6,termGroups=[copy.deepcopy(feature) for _ in range(6)])
        with self.assertRaisesRegex(ValueError,'removed mucosal entry-boundary'):
            bind_features(features,ancestry,source,carrier,'a'*64)

    def test_builder_import_has_no_command_or_file_side_effects(self):
        path=ROOT/'server/build-tongue-coupled-candidate.py'
        spec=importlib.util.spec_from_file_location('tongue_coupled_test',path)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(len(module.REQUIRED),14)

    def test_job_preparation_import_is_inert(self):
        path=ROOT/'server/prepare-tongue-coupled-install.py'
        spec=importlib.util.spec_from_file_location('tongue_coupled_install_test',path)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(module.GG,{'ZA:Genioglossus muscle.l','ZA:Genioglossus muscle.r'})

    def test_cut_subdivision_requires_deliberate_carrier_index_features(self):
        source=Mesh(np.array([[0.,1.6,0.],[.02,1.6,0.],[0.,1.62,0.]]),np.array([[0,1,2]]),FRAME)
        bary=np.array([[1.,0.,0.],[.5,.5,0.],[.5,0.,.5]])
        carrier=Mesh(bary@source.triangles[0],np.array([[0,1,2]]),FRAME)
        ancestry=dict(schema='tongue-covering-face-ancestry/v1',carrierPartId=MUCO,carrierBankSha256='a'*64,
            sourceGeometrySha256=source.geometry_sha,carrierGeometrySha256=carrier.geometry_sha,
            faces=[dict(sourceFace=0,carrierFace=0,sourceBarycentric=bary.tolist())])
        group=dict(partIds=[MUCO],surface=dict(partId=MUCO,faces=[0]),
            anchor=dict(partId=MUCO,face=0,barycentric=[1/3]*3,point=carrier.triangles[0].mean(0).tolist()))
        features=dict(schema='tongue-construction-features/v1',hostPartId=MUCO,additionalSolidMeshes=0,
            terms=6,termGroups=[copy.deepcopy(group)for _ in range(6)])
        with self.assertRaisesRegex(ValueError,'explicitly regenerated carrier-index'):
            bind_features(features,ancestry,source,carrier,'a'*64)
        features.update(featureFaceSpace='carrier',exposureRegenerated=True,carrierGeometrySha256=carrier.geometry_sha)
        result=bind_features(features,ancestry,source,carrier,'a'*64)
        self.assertEqual(result[0]['surface']['sourceBankSha256'],'a'*64)
        self.assertEqual(result[0]['surface']['faces'],[0])

    def test_actual_regenerated_controls_bind_all_cut_facets(self):
        stage=ROOT/'.local-models/geometry-workbench/tongue-completion/surface-controls-exposed-v1'
        if not stage.exists():
            self.skipTest('Local exposed metadata is unavailable')
        features=json.loads((stage/'surface-controls.json').read_text())
        ancestry=json.loads((stage/'face-ancestry.json').read_text())
        source=Mesh.load(features['sourceWholeMesh'])
        actual=Mesh.load(features['carrierMesh'])
        result=bind_features(features,ancestry,source,actual,ancestry['carrierBankSha256'])
        self.assertEqual(len(result),6)
        self.assertTrue(any('sourceBarycentric'in r for r in ancestry['faces']))
        audit=json.loads((stage/'controls-audit.json').read_text())
        self.assertTrue(audit['apexLiteral'])
        self.assertEqual(audit['marginPoints'],92)
        inferior=next(r for r in audit['decisions']if r['key']=='inferior-surface')
        self.assertIn(114,inferior['excludedSourceFaces'])
        self.assertFalse(inferior['previousAnchorStillExposed'])

    def test_actual_three_cell_bank_has_two_literal_reciprocal_charts(self):
        stage=ROOT/'.local-models/geometry-workbench/tongue-completion/shared-coverings-exposed-apo-sl-v3'
        if not stage.exists():
            self.skipTest('Local shared bank is unavailable')
        proof=json.loads((stage/'construction.json').read_text())
        meshes={pid:Mesh.load(stage/'packed'/(pid.replace(':','_')+'.npz')) for pid in
                ('TA98TONGUE:male:mucosa','TA98TONGUE:male:aponeurosis',
                 'TA98TONGUE:male:intrinsic:superior-longitudinal')}
        self.assertEqual(len(proof['fixedBoneChecks']),24)
        self.assertTrue(proof['outerMucosaCoordinatesAndWordsLiteral'])
        self.assertLessEqual(proof['maximumInnerControlDisplacementM'],1e-8)
        self.assertEqual(proof['actualPositionReencodingErrorM'],0)
        self.assertEqual([len(r['parameters']['first_faces'])for r in proof['interfaces']],[62,122])
        for row in proof['interfaces']:
            actual=check_constraint(dict(kind='shared-interface',parameters=row['parameters']),
                                    meshes[row['first']],meshes[row['second']])
            self.assertTrue(actual['passed'])
            self.assertTrue(actual['literal_opposed_normal_words'])
        apo=meshes['TA98TONGUE:male:aponeurosis']
        original=Mesh.load(proof['APO_SLSourceSurfacesAndRawWordsLiteral']['aponeurosis']['sourceMesh'])
        self.assertEqual(apo.surface_sha,original.surface_sha)
        mapping=proof['APO_SLSourceSurfacesAndRawWordsLiteral']['aponeurosis']['sourceFaceToActualFace']
        self.assertEqual(len(set(mapping)),len(original.faces))

    def test_controls_rebound_to_actual_shared_bank(self):
        stage=ROOT/'.local-models/geometry-workbench/tongue-completion/surface-controls-shared-coverings-v1'
        if not stage.exists():
            self.skipTest('Local shared-bank feature metadata is unavailable')
        features=json.loads((stage/'surface-controls.json').read_text())
        ancestry=json.loads((stage/'face-ancestry.json').read_text())
        actual=Mesh.load(features['carrierMesh'])
        source=Mesh.load(features['sourceWholeMesh'])
        result=bind_features(features,ancestry,source,actual,ancestry['carrierBankSha256'])
        self.assertEqual(len(result),6)
        self.assertEqual(result[0]['surface']['sourceBankSha256'],ancestry['carrierBankSha256'])

    def test_no_implicit_native_hg_fallback(self):
        manifest=dict(schema='tongue-selected-extrinsic-candidates/v1',passed=True,
            sourceCatalogueSha256='source',fixedBonesModified=False,
            parts=[dict(originalPartId='ZA:Genioglossus muscle.l'),dict(originalPartId='ZA:Genioglossus muscle.r')])
        with self.assertRaisesRegex(ValueError,'old native HG is not an implicit fallback'):
            selected_extrinsics_gate(manifest,'source',ROOT/'.local-models')

    def test_failed_hg_producer_cannot_be_selected(self):
        manifest=dict(schema='tongue-selected-extrinsic-candidates/v1',passed=False,
            sourceCatalogueSha256='source',fixedBonesModified=False,parts=[])
        with self.assertRaisesRegex(ValueError,'actually passed'):
            selected_extrinsics_gate(manifest,'source',ROOT/'.local-models')

    def test_actual_selected_gg_and_hgv7_producer_bindings(self):
        path=ROOT/'.local-models/geometry-workbench/tongue-completion/hyoglossus-origin-v7-nonplanar-rim-topology/selected-extrinsics.json'
        if not path.exists():
            self.skipTest('Local selected extrinsics are unavailable')
        selected=json.loads(path.read_text())
        actual=selected_extrinsics_gate(selected,selected['sourceCatalogueSha256'],ROOT/'.local-models')
        self.assertEqual(set(actual),{'ZA:Genioglossus muscle.l','ZA:Genioglossus muscle.r',
                                     'ZA:Hyoglossus muscle.l','ZA:Hyoglossus muscle.r'})

    def test_actual_complete_twenty_and_twenty_four_audits(self):
        stage=ROOT/'.local-models/geometry-workbench/tongue-completion/selected-coupling-audit-v2'
        if not stage.exists():
            self.skipTest('Final local coupling audit is unavailable')
        cover=json.loads((stage/'covering-muscle/summary.json').read_text())
        sources={r['key']:r for r in cover['sources']}
        pending=json.loads((stage/'coupled-inputs.pending.json').read_text())
        row=pending['covering_muscle_shared_interfaces'][0]
        first=Mesh.load(sources[row['cover']]['path']);second=Mesh.load(sources[row['muscle']]['path'])
        actual=check_constraint(dict(kind='shared-interface',parameters=dict(first_faces=row['cover_faces'],
            second_faces=row['muscle_faces'])),first,second)
        shared={(row['cover'],row['muscle']):dict(row,actual_check=actual)}
        checked=covering_muscle_gate(cover,{k:sources[k]for k in ('mucosa','aponeurosis')},
            {k:v for k,v in sources.items()if k not in ('mucosa','aponeurosis')},shared_interfaces=shared)
        self.assertEqual(checked['pairs'],20)
        retained=json.loads((stage/'retained-context/summary.json').read_text())
        checked=retained_envelope_gate(retained,{r['key']:r for r in retained['intrinsic']},
                                      {r['key']:r for r in retained['context']})
        self.assertEqual(checked['pairs'],24)
        self.assertFalse(checked['tissue_ownership_claim'])

    def test_actual_shared_tv_and_all_four_producer_storage_contracts(self):
        work=ROOT/'.local-models/geometry-workbench/tongue-completion'
        manifest=work/'aponeurosis-sl-storage-combinatorial-v1/construction.json'
        selected_path=work/'hyoglossus-origin-v7-nonplanar-rim-topology/selected-extrinsics.json'
        if not manifest.exists()or not selected_path.exists():
            self.skipTest('Actual producer inputs are unavailable')
        spec=importlib.util.spec_from_file_location('tongue_actual_adapter_test',ROOT/'server/build-tongue-coupled-candidate.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        original=json.loads(manifest.read_text());before=json.dumps(original,sort_keys=True)
        adapted=module.qualified_intrinsics(original,manifest)
        shared=next(r for r in adapted['meshes']if len(r.get('terms',[]))==2)
        self.assertIn('one approximate',shared['coverageLimitation'].lower())
        self.assertEqual(shared['coverageQualificationSource']['representation_mode'],module.MODE)
        self.assertEqual(json.dumps(original,sort_keys=True),before,'Adapter may not mutate producer evidence')
        self.assertEqual(len(module.old.intrinsic_terms(shared)),2)
        selected=json.loads(selected_path.read_text())
        for row in selected['parts']:
            mesh=Mesh.load(row['mesh'])
            proof=module.stored_extrinsic_replay(mesh,row)
            self.assertEqual(proof['schema'],'local-float32-faces-attribute-normals')
            self.assertTrue(proof['actual_stored_replay_passed'])

    def test_actual_root_selected_preflight_without_packing(self):
        path=ROOT/'.local-models/geometry-workbench/tongue-completion/selected-coupling-audit-v2/coupled-inputs.selected.json'
        if not path.exists():
            self.skipTest('Root prepack selection has not been recorded')
        spec=importlib.util.spec_from_file_location('tongue_actual_preflight_test',ROOT/'server/build-tongue-coupled-candidate.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        actual=module.preflight(path)
        self.assertTrue(actual['cover_gate']['passed'])
        self.assertTrue(actual['retained_gate']['passed'])
        self.assertEqual(len(actual['extrinsics']),4)


if __name__=='__main__':
    unittest.main()
