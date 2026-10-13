"""Movement closure separates independent mirrored anatomical instances."""
import importlib.util,sys,json,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('forest',ROOT/'scripts/build-spatial-dependency-forest.py');forest=importlib.util.module_from_spec(spec);spec.loader.exec_module(forest)
class SourceInstances(unittest.TestCase):
    def test_material_fragments_follow_only_the_selected_instance(self):
        parts=[dict(id=f'O3M:{side}:{material}',name='Root '+side,sourceId='Publisher right root',system='nervous',provenance={'label':'Publisher nerve','url':'https://example.invalid/nerve','mirrored':side=='left'}) for side in ['left','right'] for material in ['coat','core']]
        f=forest.build_forest('male',{'parts':parts});ids=[p['id'] for p in parts]
        plan=forest.movement_plan(f,['O3M:left:coat'],ids,np.empty((0,2),dtype=np.uint32))
        self.assertEqual(plan['movingPartIds'],['O3M:left:coat','O3M:left:core'])
        self.assertEqual(f['partToGroup']['O3M:right:coat'],f['partToGroup']['O3M:right:core'])
        self.assertNotEqual(f['partToGroup']['O3M:left:coat'],f['partToGroup']['O3M:right:coat'])
        frame=f['forest'][0]['id']
        whole=forest.movement_plan(f,[],ids,np.empty((0,2),dtype=np.uint32),[frame])
        self.assertEqual(set(whole['movingPartIds']),set(ids))
    def test_actual_left_root_scope_never_drags_right_root_or_cord(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text())
        f=forest.build_forest('male',a);ids=[p['id'] for p in a['parts']]
        selected=['ZA:Median nerve.l']+['O3M:'+name+'.l' for name in ['Lateral root of median nerve','Medial root of median nerve','Lateral cord of brachial plexus','Medial cord of brachial plexus']]
        plan=forest.movement_plan(f,selected,ids,np.empty((0,2),dtype=np.uint32))
        self.assertEqual(set(plan['movingPartIds']),set(selected))
    def test_installed_gallbladder_regions_and_coats_move_as_one_assembly(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text())
        f=forest.build_forest('male',a);ids=[p['id'] for p in a['parts']]
        expected=set(a['gallbladderWallIntegration']['partIds'])
        self.assertEqual(len(expected),16)
        if a.get('gallbladderFoldIntegration',{}).get('continuousNeck'):expected.add('FJ3080')
        if a.get('biliaryJunctionIntegration',{}).get('continuousConfluence'):expected.update(a['biliaryJunctionIntegration']['partIds'])
        # The installed cystic spiral fold shares the compound biliary carrier.
        self.assertEqual(a['cysticFoldIntegration']['foldPartId'],'CYSTICFOLD:male:A05.8.02.012')
        expected.add('CYSTICFOLD:male:A05.8.02.012')
        for selected in ['ZA:Gallbladder',sorted(expected)[0],'TA98PART:ZA:Gallbladder:A05.8.02.002']:
            plan=forest.movement_plan(f,[selected],ids,np.empty((0,2),dtype=np.uint32))
            self.assertEqual(set(plan['movingPartIds']),expected)
            self.assertIn('FJ3080',plan['movingPartIds'])
    def test_approved_oesophageal_wall_and_connector_move_together(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text())
        f=forest.build_forest('male',a);ids=[p['id'] for p in a['parts']]
        expected=set(a['oesophagealWallGeometry']['partIds'])
        for selected in ['ZA:Oesophagus','GASTRIC:male:cardiac-connection']+sorted(expected):
            plan=forest.movement_plan(f,[selected],ids,np.empty((0,2),dtype=np.uint32))
            self.assertEqual(set(plan['movingPartIds']),expected)
        stomach=forest.movement_plan(f,[a['gastricCompleteWall']['partIds'][0]],ids,np.empty((0,2),dtype=np.uint32))
        self.assertEqual(set(stomach['movingPartIds']),set(a['gastricCompleteWall']['partIds']))
        pylorus=forest.movement_plan(f,['GASTRIC:male:pyloric-connection'],ids,np.empty((0,2),dtype=np.uint32))
        self.assertEqual(pylorus['sourceFrameByGroup'][pylorus['mustMoveScaleTogetherGroups'][0]],'derived:gastric:male')
if __name__=='__main__':unittest.main()
