"""Current thyroid movement ownership; no catalogue or geometry writes."""
import hashlib
import gzip
import importlib.util
import json
import sys
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('forest',ROOT/'scripts/build-spatial-dependency-forest.py')
forest_module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(forest_module)


class FrozenThyroidExteriorOwnership(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # The old exterior and future-metadata fixtures must stay anchored to
        # the preserved pre-installation catalogue. They must never overwrite
        # a real installed tissue marker in their synthetic assertions.
        cls.path=ROOT/'.local-models/ta98-review/thyroid-finite-v1/male/pre-installation/male-atlas.json.gz'
        cls.disk_raw=cls.path.read_bytes()
        cls.raw=gzip.decompress(cls.disk_raw)
        assert hashlib.sha256(cls.raw).hexdigest()=='aed625c3c7e785510cac1c98c6353c22f16a17f931997edaa57d870cdb016f3d'
        cls.atlas=json.loads(cls.raw)
        assert not cls.atlas.get('thyroidWallAssembly')
        cls.forest=forest_module.build_forest('male',cls.atlas)

    def test_complete_declared_source_regions_move_together(self):
        members=set(self.atlas['thyroidSubdivision']['parts'])
        self.assertEqual(len(members),4)
        self.assertEqual({self.forest['partToGroup'][pid] for pid in members},
                         {'assembly:thyroid-exterior'})
        group=next(g for r in self.forest['forest'] for g in r['children']
                   if g['id']=='assembly:thyroid-exterior')
        self.assertEqual(set(group['activePartIds']),members)
        self.assertTrue(group['moveScaleTogether'])
        self.assertEqual(group['basis'],'explicit-whole-reference-object-regional-partition')

    def test_alternative_native_gland_and_adjoining_organs_are_not_dragged(self):
        for pid in ['ZA:Thyroid gland','ZA:Trachea','ZA:Cricoid cartilage',
                    'ZA:Inferior parathyroid gland.l','ZA:Inferior thyroid artery.l']:
            self.assertNotEqual(self.forest['partToGroup'][pid],'assembly:thyroid-exterior')
        self.assertFalse(self.forest['anatomicallyAccepted'])
        self.assertFalse(self.forest['mechanicalAttachmentCompleteness'])

    def test_movement_plan_expands_only_its_actual_four_part_assembly(self):
        ids=[p['id'] for p in self.atlas['parts']]
        plan=forest_module.movement_plan(self.forest,['TA98PART:Thyroid:left-lobe'],
                                        ids,np.empty((0,2),dtype=np.int64))
        self.assertEqual(set(plan['movingPartIds']),set(self.atlas['thyroidSubdivision']['parts']))
        self.assertFalse(plan['geometryModified'])
        self.assertFalse(plan['cameraMovementAuthorized'])
        self.assertEqual(self.path.read_bytes(),self.disk_raw)

    def test_regional_source_identity_cannot_be_replaced_by_anatomical_alias(self):
        altered=json.loads(self.raw)
        part=next(p for p in altered['parts'] if p['id']=='TA98PART:Thyroid:left-lobe')
        part['sourceId']='Unrelated donor gland'
        with self.assertRaises(AssertionError):
            forest_module.build_forest('male',altered)

    def test_future_shared_tissues_and_lobule_alias_keep_one_owner_per_mesh(self):
        # Metadata fixture only; these records do not establish actual tissue
        # geometry. The runtime geometry/clinical gates remain separate.
        altered=json.loads(self.raw)
        original_ids=altered['thyroidSubdivision']['parts']
        first=next(p for p in altered['parts'] if p['id']==original_ids[0])
        members={role:[] for role in ['capsule','parenchyma']}
        for index,pid in enumerate(original_ids):
            original=next(p for p in altered['parts'] if p['id']==pid)
            for role in members:
                tissue=dict(original)
                tissue.update(id=f'TEST:thyroid:{index}:{role}',
                              pairedTissue=True,sectionSolid=True,
                              tissueWallAssembly='thyroid:male',tissueWallLayer=role,
                              geometryOrigin=[0,0,0])
                for key in ['chunk','positions','normals','vertexCount']:
                    tissue[key]=first[key]
                altered['parts'].append(tissue)
                members[role].append(tissue['id'])
            original['suppressed']=True
            original['replacementParts']=[members[role][-1] for role in members]
        all_ids=members['capsule']+members['parenchyma']
        whole=next(c for c in altered['concepts'] if c['id']=='ZA:Thyroid gland')
        whole.update(elements=all_ids,ta98Term='A11.3.00.001')
        altered['thyroidWallAssembly']={
            'schema':'ta98-thyroid-wall-v1','sex':'male','wholeConcept':whole['id'],
            'capsuleParts':members['capsule'],'parenchymaParts':members['parenchyma'],
            'lobuleParts':members['parenchyma'][:2]}
        built=forest_module.build_forest('male',altered)
        self.assertEqual({built['partToGroup'][pid] for pid in all_ids},
                         {'assembly:thyroid-tissues'})
        group=next(g for r in built['forest'] for g in r['children']
                   if g['id']=='assembly:thyroid-tissues')
        self.assertEqual(set(group['activePartIds']),set(all_ids))
        self.assertEqual(set(group['suppressedPartIds']),set(original_ids))
        self.assertNotIn('ZA:Thyroid gland',group['partIds'])
        self.assertEqual(len(built['partToGroup']),len(self.atlas['parts'])+8)
        altered['thyroidWallAssembly']['lobuleParts']=['ZA:Trachea']
        with self.assertRaises(AssertionError):
            forest_module.build_forest('male',altered)


class InstalledThyroidTissueOwnership(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
        cls.raw=cls.path.read_bytes()
        cls.atlas=json.loads(cls.raw)
        cls.marker=cls.atlas['thyroidWallAssembly']
        cls.members=set(cls.marker['capsuleParts']+cls.marker['parenchymaParts']+cls.marker.get('stromaParts',[]))
        cls.forest=forest_module.build_forest('male',cls.atlas)
        cls.group=next(g for r in cls.forest['forest'] for g in r['children']
                       if g['id']=='assembly:thyroid-tissues')

    def test_installed_tissue_bank_has_one_complete_movement_owner(self):
        self.assertEqual(len(self.marker['capsuleParts']),4)
        self.assertEqual(len(self.marker['parenchymaParts']),6)
        self.assertEqual(len(self.marker['stromaParts']),1)
        self.assertEqual(len(self.members),11)
        self.assertEqual({self.forest['partToGroup'][pid] for pid in self.members},
                         {'assembly:thyroid-tissues'})
        self.assertEqual(set(self.group['activePartIds']),self.members)
        self.assertEqual(set(self.group['suppressedPartIds']),
                         set(self.atlas['thyroidSubdivision']['parts']))
        self.assertTrue(self.group['moveScaleTogether'])
        self.assertEqual(self.group['basis'],'explicit-thyroidWallAssembly')
        whole=next(c for c in self.atlas['concepts'] if c['id']==self.marker['wholeConcept'])
        self.assertEqual(set(whole['elements']),self.members)

    def test_installed_stroma_and_illustrative_domain_union_keep_unique_cells(self):
        concepts={c['id']:c for c in self.atlas['concepts']}
        parts={p['id']:p for p in self.atlas['parts']}
        parenchyma=set(self.marker['parenchymaParts'])
        stroma=set(self.marker['stromaParts'])
        domains=set(self.marker['illustrativeDomainParts'])
        self.assertEqual(len(domains),2)
        self.assertTrue(domains < parenchyma)
        self.assertNotIn('lobuleParts',self.marker)
        self.assertNotIn('lobulesConcept',self.marker)
        self.assertFalse(stroma & parenchyma)
        self.assertEqual(set(concepts['THYROID:male:stroma']['elements']),stroma)
        domain_union=concepts[self.marker['illustrativeDomainConcept']]
        self.assertEqual(set(domain_union['elements']),domains)
        self.assertEqual(domain_union['partOf'],self.marker['parenchymaConcept'])
        self.assertEqual(domain_union['representation'],'illustrative-parenchymal-domain-union')
        self.assertNotIn('ta98Term',domain_union)
        self.assertFalse(any(c.get('ta98Term')=='A11.3.00.009'
                             for c in self.atlas['concepts']))
        self.assertFalse(any(p.get('ta98Term')=='A11.3.00.009' and not p.get('suppressed')
                             for p in self.atlas['parts']))
        self.assertEqual(set(concepts['THYROID:male:parenchyma']['elements']),parenchyma)
        leaves=[p['id'] for root in self.forest['forest']
                for group in root['children'] for p in group['children']]
        for pid in stroma | domains:
            self.assertEqual(leaves.count(pid),1)
            self.assertEqual(parts[pid]['conceptId'],pid)
            self.assertEqual(concepts[pid]['elements'],[pid])
            self.assertEqual(self.forest['partToGroup'][pid],'assembly:thyroid-tissues')
        # The illustrative union reuses two actual parenchymal owners, never a
        # second geometry record or a second movement tree for those cells.
        alias=next(c for c in self.forest['conceptReferenceSets']
                   if c['id']=='THYROID:male:lobules')
        self.assertEqual(set(alias['partIds']),domains)
        self.assertFalse(alias['coMove'])
        self.assertEqual(self.path.read_bytes(),self.raw)

    def test_installed_tissue_movement_expands_active_cells_without_dragging_context(self):
        ids=[p['id'] for p in self.atlas['parts']]
        selected=self.marker['capsuleParts'][0]
        plan=forest_module.movement_plan(self.forest,[selected],ids,
                                        np.empty((0,2),dtype=np.int64))
        self.assertEqual(set(plan['movingPartIds']),self.members)
        self.assertFalse(set(plan['movingPartIds'])&set(self.group['suppressedPartIds']))
        for pid in ['ZA:Thyroid gland','ZA:Trachea','ZA:Cricoid cartilage',
                    'ZA:Inferior parathyroid gland.l','ZA:Superior parathyroid gland.r',
                    'ZA:Inferior thyroid artery.l']:
            self.assertNotIn(pid,plan['movingPartIds'])
        self.assertFalse(plan['geometryModified'])
        self.assertFalse(plan['cameraMovementAuthorized'])
        self.assertEqual(self.path.read_bytes(),self.raw)

    def test_installed_tissue_bank_rejects_split_pose_and_foreign_lobule_ownership(self):
        altered=json.loads(self.raw)
        part=next(p for p in altered['parts'] if p['id']==self.marker['capsuleParts'][0])
        part['geometryOrigin']=[part['geometryOrigin'][0]+.001,*part['geometryOrigin'][1:]]
        with self.assertRaises(AssertionError):
            forest_module.build_forest('male',altered)
        altered=json.loads(self.raw)
        altered['thyroidWallAssembly']['lobuleParts']=['ZA:Trachea']
        with self.assertRaises(AssertionError):
            forest_module.build_forest('male',altered)
        self.assertEqual(self.path.read_bytes(),self.raw)



if __name__=='__main__':
    unittest.main()
